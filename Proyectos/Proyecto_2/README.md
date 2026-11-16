# Proyecto 2 — Inestabilidad Rayleigh–Taylor magnetizada (MHD) y su extensión 3D a remanentes de supernova

Astrofísica Computacional 2026-II · Zeuz Espinosa

Proyecto del curso basado en el código **PLUTO** (volúmenes finitos, MHD ideal): el
problema [`Test_Problems/MHD/Rayleigh_Taylor`](https://pluto-code.github.io) — un fluido
pesado sobre uno ligero en gravedad efectiva con un campo magnético horizontal — con dos
partes:

1. **Parte 2D (guía de desarrollo)**: todo lo necesario para el problema base
   (config \#07: RK3 + PPM + transporte restringido) — teoría lineal verificada,
   receta PLUTO, plan del solver Python propio y análisis comparativo.
2. **Extensión 3D — SNR**: corrida real de la config \#09 (MHD 3D, RK2 + GLM) con
   perturbación aleatoria, como modelo local de la interfaz de un remanente de supernova
   joven (los "dedos" de eyecta de Tycho/Cas A frenados por el campo ambiental).

## Resultados clave de la corrida 3D

| Medición | Valor | Teoría |
|---|---|---|
| tasa efectiva de crecimiento de la capa de mezcla | σ_eff = 0.229 | 0.19 – 0.49 (rango lineal) |
| anisotropía espectral P(⊥B)/P(∥B) en k altos | ≈ 17 (17.3; ídem con 10 láminas) | > 1 (firma del RT magnético) |
| **control χ=0: misma razón, sin campo** | **1.0 exacto** | ≈ 1 (la anisotropía es magnética) |
| amplificación del campo en las puntas, max\|B\|/B₀ | 8.5 (p99: 4.4) | estiramiento de líneas congeladas |
| elongación de la mezcla ℓx/ℓz (envolventes 3D, t=30) | 1.50 | > 1 (modos ∥B crecen más) |
| conservación de la energía total | deriva ~7×10⁻¹¹ | ~ε de máquina |
| residuo ∇·B (GLM) | ~0.23 · B₀/Δx | ≠ 0 (precio de GLM vs CT) |
| **comparación 2D: σ (PLUTO #07 vs Python propio)** | 0.531 vs 0.507 (**4.5%**) | misma física; teoría fundamental 0.448 |
| **comparación 2D: L1(ρ) fase lineal / final** | 0.006 / 0.041 | difusión + desfasaje caótico |

## Estructura del repo

```
.
├── project_2.pdf               # enunciado del proyecto (curso)
├── Rayleigh_Taylor_MHD/        # PARTE 2D
│   ├── Guia_RT_Magnetizada.pdf # teoría completa, setup PLUTO, plan del código Python (+ .tex)
│   ├── comparacion_pluto_python.ipynb  # notebook EJECUTADO: PLUTO vs Python propio,
│   │                                    #   overlays de todas las variables, σ, L1/L2
│   ├── pluto_run_2d/           # corrida PLUTO 2D (preset #07 tal cual: RK3+PPM+CT,
│   │                           #   monomodo, 128×256, t=15) para la comparación
│   ├── python/
│   │   └── rt2d.py             # EL SOLVER PROPIO: FV + HLL + CT (Evans & Hawley) +
│   │                           #   RK2 + limitador MC (plan de la guía §8; 134 s/corrida)
│   └── figs/                   # figuras de la comparación (comp_*.png)
├── SNR_3D/                     # EXTENSIÓN 3D (remanente de supernova)
│   ├── README.md               # detalles de las corridas y cómo reproducirlas
│   ├── pluto_run/              # setup PLUTO de la corrida real
│   │   ├── init.c              #   init del test problem + parche RandomSeed (ver abajo)
│   │   ├── definitions.h       #   config #09 + USE_RANDOM_PERTURBATION YES
│   │   ├── pluto.ini           #   malla 48×96×48, t=30, dbl cada 2.0
│   │   └── corrida.log         #   log completo de la corrida (~22 min en serial)
│   ├── pluto_run_chi0/         # CORRIDA DE CONTROL: idéntica con CHI=0.0
│   │                           #   (regenerable ~22 min; aísla la anisotropía: 1.0 vs 17.3)
│   ├── analisis/
│   │   ├── teoria_snr.py       # teoría lineal 3D: mapa σ(kx,kz), λ_c, σ_max
│   │   ├── SNR_3D_analisis.ipynb   # notebook EJECUTADO: cortes, TODAS las variables,
│   │   │                           #   dedos, espectro por sectores, energías, divB,
│   │   │                           #   revisión (10 láminas, Bamp, vorticidad, control χ=0)
│   │   ├── SNR_3D_visualizacion.ipynb  # notebook EJECUTADO de lectura: envolventes 3D
│   │   │                                #   + GIFs (evolución y vista giratoria)
│   │   └── figs/               # figuras del análisis + viz3d_*.png, anim_*.gif
│   └── informe/
│       ├── Informe_SNR_3D.pdf  # EL INFORME (11 págs, + .tex): mapeo SNR, corrida,
│       │                       #   resultados, control χ=0, mapa de cumplimiento del
│       │                       #   enunciado y la comparación PLUTO vs Python
│       └── figs_obs/           # fotos NASA de Tycho y Cas A (con créditos)
```

## Requisitos

- **PLUTO 4.4** (compilado con gcc; aquí se usó serial) — `PLUTO_DIR` apuntando a la
  instalación.
- **Python ≥ 3.10** con `numpy`, `matplotlib`, `pyPLUTO` (aquí `py-pluto 1.2.6` de PyPI:
  `pip install --user "py-pluto"`), `jupyter`.
- **pdflatex** (dos pasadas) para los `.tex`.

## Cómo reproducir todo

```bash
# 0) la parte 2D: PLUTO #07 para la comparación (~2 min) + el solver Python (134 s)
cd Rayleigh_Taylor_MHD/pluto_run_2d
printf 'ARCH = Linux.gcc.defs\n' > makefile
python3 "$PLUTO_DIR/setup.py" --auto-update
echo 'CFLAGS += -D_DEFAULT_SOURCE' >> makefile
make && ./pluto
cd ../python && python3 rt2d.py
cd .. && jupyter nbconvert --to notebook --execute --inplace comparacion_pluto_python.ipynb

# 1) la corrida 3D (~22 min en un portátil; genera data.0000..0015.dbl)
cd SNR_3D/pluto_run
printf 'ARCH = Linux.gcc.defs\n' > makefile
python3 "$PLUTO_DIR/setup.py" --auto-update
echo 'CFLAGS += -D_DEFAULT_SOURCE' >> makefile   # fix srand48 en Fedora
make && ./pluto

# 2) la teoría lineal (mapa σ(kx,kz) + números + CSVs)
cd ../analisis
python3 teoria_snr.py

# 2b) la corrida de CONTROL χ=0 (~22 min más; CHI es runtime, no recompila)
mkdir ../pluto_run_chi0 && cd ../pluto_run_chi0
cp ../pluto_run/{init.c,definitions.h,pluto.ini,makefile,pluto} .
sed -i 's/^CHI  *0.2/CHI                         0.0/' pluto.ini
./pluto    # o: nohup ./pluto > corrida_chi0.log 2>&1 &
cd ../analisis

# 3) el análisis (regenera las figuras con los datos del paso 1; la sección 7
#    de revisión compara contra el control del paso 2b)
jupyter nbconvert --to notebook --execute --inplace SNR_3D_analisis.ipynb

# 4) las gráficas y animaciones 3D (envolventes + GIFs; tarda ~2 min)
jupyter nbconvert --to notebook --execute --inplace SNR_3D_visualizacion.ipynb

# 5) los PDFs (ya compilados en el repo; para rehacerlos)
cd ../informe && pdflatex Informe_SNR_3D.tex && pdflatex Informe_SNR_3D.tex
```

**Nota sobre el `init.c`**: el del test problem aborta con la perturbación aleatoria
activada (`maximum number of prns exceeded`); el de este repo lleva un parche de 5 líneas
que siembra el generador (`RandomSeed(26041988,0)` la primera vez), con semilla fija para
que la corrida sea reproducible.

## Qué NO está en el repo (y por qué)

- `pluto_run/data.*.dbl` (16 × 16 MB): salidas de la simulación, 100 % regenerables con
  el paso 1. `dbl.out` y `grid.out` sí van porque son diminutos y documentan la corrida.
- `pluto_run_2d/data.*.dbl` y `python/py_rt2d.npz` (42 MB): ídem, regenerables con el
  paso 0 en ~4 minutos entre ambos.
- `*.o`, binario `pluto`, `makefile`: artefactos de compilación (y el makefile lleva
  rutas absolutas de la máquina donde se generó).
- `analisis/*.csv`: derivados de `teoria_snr.py` en 2 segundos.
- Los documentos de estudio personales (preparación de la sustentación, glosario
  razonado, manual de la GUI) viven en una carpeta local fuera del repo.
- Todo lo demás (notebooks ejecutados, figuras, PDFs) sí está.

## Referencias principales

Chandrasekhar (1961), *Hydrodynamic and Hydromagnetic Stability* · Stone & Gardiner
(2007), Phys. Fluids 19, 094104 · Orlando et al. (2012), arXiv:1202.3593 · Mignone et al.
(2012), ApJS 198, 7 (PLUTO). La bibliografía comentada completa está en la guía 2D (§11)
y en el informe.
