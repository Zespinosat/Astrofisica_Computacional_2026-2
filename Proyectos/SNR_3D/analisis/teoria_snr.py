# teoria_snr.py -- teoria lineal de la RT magnetizada en 3D para el caso de la corrida
#
# La idea: en 3D la perturbacion vive en el plano x-z de la interfaz y el campo
# B = Bx * xhat es horizontal. La relacion de dispersion de Chandrasekhar
# (unidades PLUTO, sin 4*pi) es
#
#     sigma^2 = A*g*k - 2*(k.B)^2/(rho1+rho2),   k = sqrt(kx^2 + kz^2)
#
# O sea: solo la componente kx (paralela al campo) siente la tension magnetica.
# Los modos "de intercambio" (kz != 0, kx ~ 0) crecen como en hidrodinamica pura.
# Eso es lo que hace la version 3D distinta de la 2D y es justo lo que pasa en
# la interfaz de un remanente de supernova: los dedos crecen alargados a lo
# largo del campo (Orlando et al. 2012).
#
# Caso de la corrida: ETA=4, GRAV=-0.01, CHI=0.2 (pluto.ini de pluto_run/)
# Uso:  python3 teoria_snr.py     (guarda figs/ y .csv en esta misma carpeta)

import numpy as np
import matplotlib
matplotlib.use("Agg")          # sin display, para correrlo en terminal
import matplotlib.pyplot as plt
import os

# ---------------- params de la corrida (a mano, ojo si se cambia el ini) -----
eta = 4.0        # rho2 (arriba)
rho1, rho2 = 1.0, eta
grav = 0.01      # |g| (en el ini va con signo: -0.01)
chi = 0.2        # campo en unidades del critico Bc
Lx = 1.0         # ancho horizontal -> unidad de longitud
gamma = 5.0/3.0

# ---------------- numeros derivados ------------------------------------------
A = (rho2 - rho1) / (rho2 + rho1)              # Atwood = 0.6
Bc = np.sqrt((rho2 - rho1) * grav * Lx)        # campo critico (def del init.c)
B = chi * Bc / np.sqrt(4.0 * np.pi)            # campo en unidades de codigo
cs = 1.0                                       # por diseno (p0 = 1/gamma, rho1=1)
vA1, vA2 = B/np.sqrt(rho1), B/np.sqrt(rho2)
p0 = 1.0/gamma                                 # presion en la interfaz
beta = 2.0*p0/(B**2)                           # beta = 2p/B^2 en unidades PLUTO

def sigma(kx, kz):
    """sigma(kx,kz), solo la parte kx acoplada al campo. Regresa np.nan donde
    sigma^2 < 0 (modos estabilizados, no crecen)."""
    k = np.sqrt(kx**2 + kz**2)
    s2 = A*grav*k - 2.0*(kx*B)**2/(rho1 + rho2)
    out = np.sqrt(np.where(s2 > 0, s2, np.nan))
    return out

# ---------------- mapa sigma(kx,kz) ------------------------------------------
# OJO: para k.perp.B la teoria incompresible no tiene ningun corte (sigma crece
# como sqrt(k) hasta el infinito). El que corta de verdad es la MALLA: no hay
# modos con k > kNyq = pi/dx. Con la malla de la corrida (48 celdas en 1.0):
dx = 1.0/48.0
kNyq = np.pi/dx                                # ~ 150.8
n = 400
kmax = 1.05*kNyq                               # hasta donde alcanza a ver la malla
kx = np.linspace(0.0, kmax, n)
kz = np.linspace(0.0, kmax, n)
KX, KZ = np.meshgrid(kx, kz, indexing="ij")
S = sigma(KX, KZ)

kc = 2.0*np.pi/(chi**2 * Lx)                    # corte paralelo al campo
lam_c = chi**2 * Lx

# valores puntuales que uso en el informe
s_diag  = sigma(2*np.pi/Lx, 2*np.pi/Lx)         # modo fundamental diagonal
s_inter = sigma(0.0, 2*np.pi/Lx)                # intercambio puro (k.B = 0)
s_paral = sigma(2*np.pi/Lx, 0.0)                # undular puro (k paralelo a B)
s_max   = float(np.nanmax(S))                   # cae en la esquina k.perp.B (grid scale)
# el maximo de los modos PARALELOS si sale de la teoria con campo:
kk_par = np.linspace(1.0, kc, 2000)
sig_par = sigma(kk_par, 0.0*kk_par)
i_best = int(np.nanargmax(sig_par))
s_max_par = sig_par[i_best]
k_best_par = kk_par[i_best]

print("=" * 62)
print("teoria lineal RT magnetizada 3D -- caso de la corrida")
print("=" * 62)
print(f"  A (Atwood)             = {A:.3f}")
print(f"  Bc = sqrt((eta-1)g Lx) = {Bc:.4f}")
print(f"  B = chi*Bc/sqrt(4pi)   = {B:.5f}   (vA = {vA1:.4f} abajo, {vA2:.4f} arriba)")
print(f"  beta (interfaz)        = {beta:.3e}   -> campo globalmente debil")
print(f"  lambda_c = chi^2 Lx    = {lam_c:.3f}   (modos paralelos a B con lambda<esta no crecen)")
print(f"  sigma max paralelo a B = {s_max_par:.4f} en k~{k_best_par:.1f} (lambda~{2*np.pi/k_best_par:.3f} = 2*lambda_c)")
print(f"  sigma max de TODO el mapa = {s_max:.4f}  <-- esquina k_perp B ~ escala de malla (kNyq={kNyq:.0f})")
print(f"  sigma fundamental diag (kx=kz=2pi) = {s_diag:.4f}")
print(f"  sigma intercambio puro (kx=0)      = {s_inter:.4f}   <- sin freno magnetico")
print(f"  sigma undular puro (kz=0)          = {s_paral:.4f}")
print(f"  e-foldings hasta tstop=30 (modo paralelo max) = {s_max_par*30.0:.1f}")
print("=" * 62)

# guardo el mapa en csv por si lo quiero mirar despues (parsear en el notebook)
np.savetxt("sigma_map.csv", S, delimiter=",",
           header=f"sigma(kx,kz); kx,kz de 0 a {kmax} en {n} pts (indexing ij); nan = estable")
np.savetxt("k_grid.csv", kx, delimiter=",", header="kx = kz grid")

# ---------------- figura 1: mapa ---------------------------------------------
fig, ax = plt.subplots(figsize=(6.2, 5.4))
im = ax.pcolormesh(KX, KZ, np.ma.masked_invalid(S), cmap="magma", shading="auto")
ax.contour(KX, KZ, S, levels=[0.0], colors="cyan", linewidths=1.2)  # borde del cono estable
ax.plot(k_best_par, 0.0, "c*", ms=13,
        label=f"máx paralelo: $\\sigma$={s_max_par:.3f}")
ax.plot([0, kmax], [0, kmax], "w--", lw=0.8, alpha=0.6)
ax.axhline(kNyq, color="w", lw=0.7, ls="--", alpha=0.5)
ax.annotate("k$_{Nyq}$ = $\\pi/\\Delta x$ (la malla no ve más allá)", (4, kNyq+5),
            color="k", fontsize=8, bbox=dict(boxstyle="round", fc="white", ec="none", alpha=0.8))
ax.set_xlabel("$k_x$  (paralelo a B)")
ax.set_ylabel("$k_z$  (perpendicular)")
ax.set_title("$\\sigma(k_x,k_z)$ — zona negra: estabilizada por tensión\n"
             "(caso de la corrida: $\\eta$=4, g=0.01, $\\chi$=0.2)")
fig.colorbar(im, ax=ax, label="$\\sigma$")
ax.legend(loc="lower left")
fig.tight_layout()
fig.savefig("figs/teoria_mapa_sigma.png", dpi=150)
plt.close(fig)

# ---------------- figura 2: cortes 1D ----------------------------------------
kk = np.linspace(1.0, kmax, 800)
fig, ax = plt.subplots(figsize=(6.8, 4.2))
ax.plot(kk, sigma(0*kk, kk), color="C0", lw=2, label="k $\\perp$ B (intercambio): $\\sigma=\\sqrt{Agk}$")
ax.plot(kk, sigma(kk, 0*kk), color="C3", lw=2, label="k $\\parallel$ B (undular): freno magnético")
ax.axvline(kc, color="gray", ls=":", lw=1)
ax.annotate(f"$k_c$=2$\\pi$/λ$_c$ (λ$_c$={lam_c:.2f})", (kc, 0.9*s_max), rotation=90, fontsize=9, color="gray")
ax.axhline(s_max, color="C2", ls="--", lw=1, label=f"$\\sigma_{{max}}$={s_max:.3f} (Chandrasekhar)")
ax.set_xlabel("$k$"); ax.set_ylabel("$\\sigma$")
ax.set_title("cortes del mapa: el campo solo frena los modos paralelos")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig("figs/teoria_cortes.png", dpi=150)
plt.close(fig)

print("listo: figs/teoria_mapa_sigma.png, figs/teoria_cortes.png, sigma_map.csv")
