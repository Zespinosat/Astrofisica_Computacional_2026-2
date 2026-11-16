#!/usr/bin/env python3
"""
rt2d.py — la RT magnetizada 2D resuelta "a mano" (Proyecto 2, parte Python propio).

Reproduce el mismo problema que la corrida PLUTO 2D (config #07 del test problem):
RT monomodo con campo horizontal, ETA=2, GRAV=-0.1, CHI=0.2, malla 128x256
en x in [-0.5,0.5], y in [-1,1], hasta t=15.

El plan es el de la guia (seccion 8): volumenes finitos + solver de Riemann HLL
+ transporte restringido (CT, Evans & Hawley con promedio aritmetico) + RK2 +
limitador MC. Lo mismo que hace PLUTO, pero de segundo orden "casero" en vez de
PPM+RK3, y HLL en vez de Roe. La comparacion contra PLUTO vive en el notebook
comparacion_pluto_python.ipynb (carpeta de arriba).

Uso:
    python3 rt2d.py            # corrida completa (guarda py_rt2d.npz)
    python3 rt2d.py 1.0        # smoke test hasta t=1.0

Salida (py_rt2d.npz): t (Nt,), campos (Nt, NY, NX): rho, prs, vx, vy, bxc, byc,
mas xc, yc.
"""
import numpy as np
import time, sys, os

# ----------------------------------------------------------------------
# parametros (los del pluto.ini de la corrida 2D; si cambio el ini cambio esto)
# ----------------------------------------------------------------------
ETA, GRAV, CHI = 2.0, -0.1, 0.2          # densidad de arriba, gravedad, campo/critico
GAMMA = 5.0 / 3.0
NX, NY = 128, 256                        # misma malla que el preset 07
X1, X2 = -0.5, 0.5                       # rango en x (Lx = 1)
Y1, Y2 = -1.0, 1.0                       # rango en y
TSTOP = 15.0
DT_OUT = 0.5                             # un cuadro cada 0.5, como el dbl del ini
CFL = 0.35                               # conservador; el CT casero prefiere margen
P_FLOOR = 1.0e-8                         # presion minima (no deberia activarse)

dx = (X2 - X1) / NX
dy = (Y2 - Y1) / NY
xc = X1 + (np.arange(NX) + 0.5) * dx
yc = Y1 + (np.arange(NY) + 0.5) * dy
XC, YC = np.meshgrid(xc, yc)             # shape (NY, NX): filas=y, columnas=x

# campo inicial: B0 = chi*sqrt((eta-1)|g|)/sqrt(4pi), la misma formula del init.c
B0 = CHI * np.sqrt((ETA - 1.0) * abs(GRAV)) / np.sqrt(4.0 * np.pi)

# ----------------------------------------------------------------------
# estado: fluido en conservativas U = [rho, mx, my, E] (4, NY, NX)
#         campo en las CARAS: bxf (NX+1, NY), byf (NX, NY+1)  <- por el CT
# ----------------------------------------------------------------------
def caras_a_celdas(bxf, byf):
    "promedio de caras al centro (lo usan el HLL y las salidas)"
    bxc = 0.5 * (bxf[:-1, :] + bxf[1:, :]).T      # (NY, NX); bxf es (NX+1, NY)
    byc = 0.5 * (byf[:, :-1] + byf[:, 1:]).T      # (NY, NX)
    return bxc, byc


def condicion_inicial():
    rho = np.where(YC < 0.0, 1.0, ETA)            # pesado arriba (la "eyecta" 2D)
    prs = 1.0 / GAMMA + rho * GRAV * YC           # equilibrio hidrostatico del init.c
    vx = np.zeros_like(rho)
    vy = -1.0e-2 * (1.0 + np.cos(2.0 * np.pi * XC)) * np.exp(-YC**2 * 200.0)
    bxf = np.full((NX + 1, NY), B0)               # Bx constante -> divB = 0 exacto
    byf = np.zeros((NX, NY + 1))
    bxc, byc = caras_a_celdas(bxf, byf)
    E = (prs / (GAMMA - 1.0) + 0.5 * rho * (vx**2 + vy**2)
         + 0.5 * (bxc**2 + byc**2))
    return np.array([rho, rho * vx, rho * vy, E]), bxf, byf


def conservadas_a_primitivas(U, bxf, byf):
    rho, mx, my, E = U
    bxc, byc = caras_a_celdas(bxf, byf)
    rho = np.maximum(rho, 1.0e-12)
    vx, vy = mx / rho, my / rho
    prs = (E - 0.5 * rho * (vx**2 + vy**2) - 0.5 * (bxc**2 + byc**2)) * (GAMMA - 1.0)
    prs = np.maximum(prs, P_FLOOR)
    return rho, vx, vy, prs

# ----------------------------------------------------------------------
# reconstruccion LINEAR + limitador MC (sin limitar, 2do orden oscila: Godunov)
# ----------------------------------------------------------------------
def mc_slope(q, qm, qp):
    "limitador monotonized-central de van Leer / Woodward & Colella"
    dqm, dqp = q - qm, qp - q
    mag = np.minimum(np.minimum(np.abs(0.5 * (dqm + dqp)), 2 * np.abs(dqm)),
                     2 * np.abs(dqp))
    return np.where(dqm * dqp > 0.0, np.sign(dqp) * mag, 0.0)


def vecinos_y(q, par):
    "vecino de arriba/abajo en y con pared reflectiva: espejo con paridad"
    Q = np.concatenate([q[1::-1], q, q[-1:-3:-1]], axis=0)   # [q1,q0, ..., q-1,q-2]
    if par:                                                   # impar: vy, By
        Q[0] *= -1; Q[1] *= -1; Q[-1] *= -1; Q[-2] *= -1
    return Q[1:NY + 1], Q[3:NY + 3]                           # qm, qp de cada celda


def estados(q, eje, par=0):
    "devuelve (lado_menor, lado_mayor) de cada CELDA: q -/+ 0.5*pendiente"
    if eje == "x":                     # periodico: un roll y listo
        qm, qp = np.roll(q, 1, axis=1), np.roll(q, -1, axis=1)
    else:
        qm, qp = vecinos_y(q, par)
    s = mc_slope(q, qm, qp)
    return q - 0.5 * s, q + 0.5 * s

# ----------------------------------------------------------------------
# HLL (Harten-Lax-van Leer). U y F en el MISMO orden:
# [rho, m_normal, m_tangencial, E, B_tangencial]. El B_tangencial del flujo es
# la FEM que usa el CT (con su signo) — el campo NO se actualiza por aca.
# ----------------------------------------------------------------------
def cf2(rho, prs, bn, bt):
    "magnetosonica rapida al cuadrado (con el B de esa celda)"
    cs2 = GAMMA * prs / rho
    va2 = (bn**2 + bt**2) / rho
    a = np.maximum(0.25 * (cs2 + va2)**2 - cs2 * bn**2 / rho, 0.0)
    return 0.5 * (cs2 + va2) + np.sqrt(a)


def hll(rL, vnL, vtL, pL, bnL, btL, rR, vnR, vtR, pR, bnR, btR, flujo_f):
    "5 filas de flujo HLL por cara: [rho, mn, mt, E, Bt]"
    cfL, cfR = np.sqrt(cf2(rL, pL, bnL, btL)), np.sqrt(cf2(rR, pR, bnR, btR))
    SL = np.minimum(vnL - cfL, vnR - cfR)
    SR = np.maximum(vnL + cfL, vnR + cfR)
    FL = flujo_f(rL, vnL, vtL, pL, bnL, btL)
    FR = flujo_f(rR, vnR, vtR, pR, bnR, btR)

    def energ(rho, vn, vt, prs, bn, bt):
        return (prs / (GAMMA - 1.0) + 0.5 * rho * (vn**2 + vt**2)
                + 0.5 * (bn**2 + bt**2))

    UL = np.stack([rL, rL * vnL, rL * vtL, energ(rL, vnL, vtL, pL, bnL, btL), btL])
    UR = np.stack([rR, rR * vnR, rR * vtR, energ(rR, vnR, vtR, pR, bnR, btR), btR])
    with np.errstate(divide="ignore", invalid="ignore"):
        FH = (SR * FL - SL * FR + SL * SR * (UR - UL)) / (SR - SL)
    return np.where(SL >= 0.0, FL, np.where(SR <= 0.0, FR, FH))


def flujo_x_f(rho, vn, vt, prs, bn, bt):
    "direccion x: normal=vx, tang=vy, Bn=Bx, Bt=By"
    pT = prs + 0.5 * (bn**2 + bt**2)
    E = (prs / (GAMMA - 1.0) + 0.5 * rho * (vn**2 + vt**2) + 0.5 * (bn**2 + bt**2))
    return np.stack([rho * vn,
                     rho * vn**2 + pT - bn**2,
                     rho * vn * vt - bn * bt,
                     (E + pT) * vn - bn * (vn * bn + vt * bt),
                     vt * bn - vn * bt])            # F(By) = -Ez


def flujo_y_f(rho, vn, vt, prs, bn, bt):
    "direccion y: normal=vy, tang=vx, Bn=By, Bt=Bx"
    pT = prs + 0.5 * (bn**2 + bt**2)
    E = (prs / (GAMMA - 1.0) + 0.5 * rho * (vn**2 + vt**2) + 0.5 * (bn**2 + bt**2))
    return np.stack([rho * vn,
                     rho * vn**2 + pT - bn**2,
                     rho * vn * vt - bn * bt,
                     (E + pT) * vn - bn * (vn * bn + vt * bt),
                     vt * bn - vn * bt])            # F(Bx) = +Ez

# ----------------------------------------------------------------------
# rhs completo: dU/dt (4,NY,NX) y derivadas de los campos de caras
# ----------------------------------------------------------------------
def rhs(U, bxf, byf):
    rho, vx, vy, prs = conservadas_a_primitivas(U, bxf, byf)
    bxc, byc = caras_a_celdas(bxf, byf)

    # ---- barrido x (periodico): Bn = Bx del campo escalonado, NO se reconstruye
    rA, rB = estados(rho, "x"); vxA, vxB = estados(vx, "x")
    vyA, vyB = estados(vy, "x"); pA, pB = estados(prs, "x")
    byA, byB = estados(byc, "x")                    # By tangencial: si se reconstruye
    # cara derecha de la celda i: L=lado mayor de i, R=lado menor de i+1
    rR = np.roll(rA, -1, axis=1); vxR = np.roll(vxA, -1, axis=1)
    vyR = np.roll(vyA, -1, axis=1); pR = np.roll(pA, -1, axis=1)
    byR = np.roll(byA, -1, axis=1)
    bn = bxf[1:, :].T                               # Bx de la cara derecha de cada celda
    Fx = hll(rB, vxB, vyB, pB, bn, byB,
             rR, vxR, vyR, pR, bn, byR, flujo_x_f)  # (5, NY, NX) caras derechas
    # divergencia ciclica: cara derecha de i es Fx[...,i], la izquierda es la de
    # i-1 = roll(+1). Ojo con el signo del roll (este bug tuyo exploto en 50 pasos)
    dUx = -(Fx[:4] - np.roll(Fx[:4], 1, axis=2)) / dx
    # FEM en cada cara x: ojo, la fila 5 del flujo es F(By) = -Ez -> hay que
    # invertir el signo (me costo dos corridas explotadas encontrarlo)
    Ez_x_caras = -Fx[4]                             # (NY, NX) caras derechas
    EMFx_full = np.concatenate([Ez_x_caras[:, -1:], Ez_x_caras], axis=1).T  # (NX+1, NY)

    # ---- barrido y (reflectivo): Bn = By de las caras interiores
    rA, rB = estados(rho, "y"); vyA, vyB = estados(vy, "y", 1)
    vxA, vxB = estados(vx, "y"); pA, pB = estados(prs, "y")
    bxA, bxB = estados(bxc, "y")                    # Bx tangencial
    # caras INTERIORES (entre fila j y j+1): L=lado mayor de j, R=lado menor de j+1
    j = slice(0, NY - 1)
    rR = np.roll(rA, -1, axis=0)[j]; vyR = np.roll(vyA, -1, axis=0)[j]
    vxR = np.roll(vxA, -1, axis=0)[j]; pR = np.roll(pA, -1, axis=0)[j]
    bxR = np.roll(bxA, -1, axis=0)[j]
    bn = byf[:, 1:NY].T                             # By en las caras interiores (NY-1,NX)
    Fy = hll(rB[j], vyB[j], vxB[j], pB[j], bxB[j], bn,
             rR, vyR, vxR, pR, bxR, bn, flujo_y_f)  # (5, NY-1, NX)
    # caras de la PARED: flujo normal de masa/energia 0 y FEM 0, PERO el momento
    # normal lleva la presion del "fantasma" (si la pongo en 0 la fila de celdas
    # del borde se cae contra la pared: pared = espejo con vn=0, By=0)
    Fy_p = np.zeros((5, NY + 1, NX))
    Fy_p[:, 1:NY, :] = Fy
    p_pared = prs + 0.5 * bxc**2            # p + p_mag con By=0 en el espejo
    Fy_p[1, 0, :] = p_pared[0, :]
    Fy_p[1, -1, :] = p_pared[-1, :]
    dUy = -(Fy_p[:, 1:, :] - Fy_p[:, :-1, :]) / dy  # filas: [rho, d(my), d(mx), E, .]
    EMFy_full = Fy_p[4].T                           # Ez en caras y (NX, NY+1)

    # ---- ensamblo dU (recordando que flujo_y tiene my primero y mx despues)
    dU = dUx.copy()
    dU[0] += dUy[0]
    dU[1] += dUy[2]                                 # mx
    dU[2] += dUy[1]                                 # my
    dU[3] += dUy[3]
    dU[2] += rho * GRAV                             # fuente: gravedad en my
    dU[3] += rho * vy * GRAV                        # y su trabajo en E

    # ---- CT: Ez en las esquinas = promedio de las 4 FEM vecinas (Evans & Hawley)
    EzC = np.zeros((NX + 1, NY + 1))
    EzC[1:-1, 1:-1] = 0.25 * (EMFx_full[1:-1, 1:] + EMFx_full[1:-1, :-1]
                              + EMFy_full[1:, 1:-1] + EMFy_full[:-1, 1:-1])
    EzC[0, :] = EzC[-2, :]; EzC[-1, :] = EzC[1, :]  # periodicidad en x
    dbxf = -(EzC[:, 1:] - EzC[:, :-1]) / dy         # dBxf/dt = -dEz/dy
    dbyf = (EzC[1:, :] - EzC[:-1, :]) / dx          # dByf/dt = +dEz/dx
    return dU, dbxf, dbyf

# ----------------------------------------------------------------------
# RK2 (Heun) + CFL + diagnostico de divB
# ----------------------------------------------------------------------
def paso_temporal(U, bxf, byf, dt):
    dU1, dbx1, dby1 = rhs(U, bxf, byf)
    U1, bx1, by1 = U + dt * dU1, bxf + dt * dbx1, byf + dt * dby1
    dU2, dbx2, dby2 = rhs(U1, bx1, by1)
    return (U + 0.5 * dt * (dU1 + dU2),
            bxf + 0.5 * dt * (dbx1 + dbx2),
            byf + 0.5 * dt * (dby1 + dby2))


def paso_dt(rho, prs, bxc, byc, vx, vy):
    cfx = np.sqrt(cf2(rho, prs, bxc, byc))
    cfy = np.sqrt(cf2(rho, prs, byc, bxc))
    return CFL * min(dx / np.max(np.abs(vx) + cfx),
                     dy / np.max(np.abs(vy) + cfy))


def divB_max(bxf, byf):
    "divergencia discreta del campo de caras (deberia dar ~1e-16 siempre)"
    div = ((bxf[1:, :] - bxf[:-1, :]) / dx).T + ((byf[:, 1:] - byf[:, :-1]) / dy).T
    return np.max(np.abs(div)) * dx / B0

# ----------------------------------------------------------------------
# bucle principal
# ----------------------------------------------------------------------
def correr(tstop=TSTOP):
    t, nsteps = 0.0, 0
    U, bxf, byf = condicion_inicial()
    sal_t, campos = [], []

    def guardar(t):
        rho, vx, vy, prs = conservadas_a_primitivas(U, bxf, byf)
        bxc, byc = caras_a_celdas(bxf, byf)
        sal_t.append(t)
        campos.append([rho, prs, vx, vy, bxc, byc])

    guardar(0.0)
    t0 = time.time()
    while t < tstop - 1e-12:
        rho, vx, vy, prs = conservadas_a_primitivas(U, bxf, byf)
        bxc, byc = caras_a_celdas(bxf, byf)
        dt = min(paso_dt(rho, prs, bxc, byc, vx, vy), tstop - t)
        U, bxf, byf = paso_temporal(U, bxf, byf, dt)
        t += dt
        nsteps += 1
        if t - sal_t[-1] >= DT_OUT - 1e-9 or t >= tstop - 1e-12:
            guardar(t)
        if nsteps % 500 == 0:
            print(f"  t={t:6.2f}  dt={dt:.2e}  pasos={nsteps}  "
                  f"divB={divB_max(bxf, byf):.1e}  ({time.time()-t0:.0f}s)",
                  flush=True)

    print(f"listo: {nsteps} pasos en {time.time()-t0:.0f}s")
    return sal_t, campos


if __name__ == "__main__":
    tstop = float(sys.argv[1]) if len(sys.argv) > 1 else TSTOP
    sal_t, campos = correr(tstop)
    salida = {
        "t": np.array(sal_t),
        "rho": np.array([c[0] for c in campos]),
        "prs": np.array([c[1] for c in campos]),
        "vx": np.array([c[2] for c in campos]),
        "vy": np.array([c[3] for c in campos]),
        "bxc": np.array([c[4] for c in campos]),
        "byc": np.array([c[5] for c in campos]),
        "xc": xc, "yc": yc, "gamma": GAMMA, "B0": B0,
    }
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), "py_rt2d.npz")
    np.savez_compressed(fn, **salida)
    print("guardado:", fn)
