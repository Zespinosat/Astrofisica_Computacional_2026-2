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
| anisotropía espectral P(⊥B)/P(∥B) en k altos | ≈ 17 | > 1 (firma del RT magnético) |
| conservación de la energía total | deriva ~7×10⁻¹¹ | ~ε de máquina |
| residuo ∇·B (GLM) | ~0.23 · B₀/Δx | ≠ 0 (precio de GLM vs CT) |

## Estructura del repo

```
.
├── project_2.pdf               # enunciado del proyecto (curso)
├── Rayleigh_Taylor_MHD/        # PARTE 2D
│   ├── Guia_RT_Magnetizada.pdf # guía completa: teoría, setup PLUTO, plan Python (+ .tex)
├── SNR_3D/                     # EXTENSIÓN 3D (remanente de supernova)
│   ├── README.md               # detalles de la corrida y cómo reproducirla
│   ├── pluto_run/              # setup PLUTO de la corrida real
│   │   ├── init.c              #   init del test problem + parche RandomSeed (ver abajo)
│   │   ├── definitions.h       #   config #09 + USE_RANDOM_PERTURBATION YES
│   │   ├── pluto.ini           #   malla 48×96×48, t=30, dbl cada 2.0
│   │   ├── corrida.log         #   log completo de la corrida (~22 min en serial)
│   │   └── (data.*.dbl NO va al repo: 245 MB, regenerables)
│   ├── analisis/
│   │   ├── teoria_snr.py       # teoría lineal 3D: mapa σ(kx,kz), λ_c, σ_max
│   │   ├── SNR_3D_analisis.ipynb  # notebook EJECUTADO: cortes, dedos, espectro,
│   │   │                          #   anisotropía por sectores, energías, divB
│   │   └── figs/               # las 9 figuras (usan también los PDFs)
│   └── informe/
│       ├── Informe_SNR_3D.pdf  # informe de la extensión 3D (9 págs, + .tex)
│       └── Glosario_Conceptos.pdf  # glosario razonado: conceptos, términos y
│                                    #   ecuaciones, con deducciones (+ .tex)
```

## Requisitos

- **PLUTO 4.4** (compilado con gcc; aquí se usó serial) — `PLUTO_DIR` apuntando a la
  instalación.
- **Python ≥ 3.10** con `numpy`, `matplotlib`, `pyPLUTO` (aquí `py-pluto 1.2.6` de PyPI:
  `pip install --user "py-pluto"`), `jupyter`.
- **pdflatex** (dos pasadas) para los `.tex`.

## Cómo reproducir todo

```bash
# 1) la corrida 3D (~22 min en un portátil; genera data.0000..0015.dbl)
cd SNR_3D/pluto_run
printf 'ARCH = Linux.gcc.defs\n' > makefile
python3 "$PLUTO_DIR/setup.py" --auto-update
echo 'CFLAGS += -D_DEFAULT_SOURCE' >> makefile   # fix srand48 en Fedora
make && ./pluto

# 2) la teoría lineal (mapa σ(kx,kz) + números + CSVs)
cd ../analisis
python3 teoria_snr.py

# 3) el análisis (regenera las figuras con los datos del paso 1)
jupyter nbconvert --to notebook --execute --inplace SNR_3D_analisis.ipynb

# 4) los PDFs (ya compilados en el repo; para rehacerlos)
cd ../informe && pdflatex Informe_SNR_3D.tex && pdflatex Informe_SNR_3D.tex
pdflatex Glosario_Conceptos.tex && pdflatex Glosario_Conceptos.tex
```

**Nota sobre el `init.c`**: el del test problem aborta con la perturbación aleatoria
activada (`maximum number of prns exceeded`); el de este repo lleva un parche de 5 líneas
que siembra el generador (`RandomSeed(26041988,0)` la primera vez), con semilla fija para
que la corrida sea reproducible.

## Qué NO está en el repo (y por qué)

- `pluto_run/data.*.dbl` (16 × 16 MB): salidas de la simulación, 100 % regenerables con
  el paso 1. `dbl.out` y `grid.out` sí van porque son diminutos y documentan la corrida.
- `*.o`, binario `pluto`, `makefile`: artefactos de compilación (y el makefile lleva
  rutas absolutas de la máquina donde se generó).
- `analisis/*.csv`: derivados de `teoria_snr.py` en 2 segundos.
- Todo lo demás (notebook ejecutado, figuras, PDFs) sí está.

## Referencias principales

Chandrasekhar (1961), *Hydrodynamic and Hydromagnetic Stability* · Stone & Gardiner
(2007), Phys. Fluids 19, 094104 · Orlando et al. (2012), arXiv:1202.3593 · Mignone et al.
(2012), ApJS 198, 7 (PLUTO). La bibliografía comentada completa está en la guía 2D (§11)
y en el glosario.
