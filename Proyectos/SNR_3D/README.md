# SNR_3D — RT magnetizada 3D como modelo local de la interfaz de un remanente de supernova

Carpeta de la **extensión 3D** del Proyecto 2 (MHD). El problema de base sigue siendo
`Test_Problems/MHD/Rayleigh_Taylor` de PLUTO, pero en su configuración **#09 (MHD 3D,
RK2 + LINEAR + GLM)** y con la perturbación **aleatoria** activada, para que en la interfaz
salgan "dedos" tipo Tycho/Cas A en vez de una sola burbuja.

## Contenido

| Carpeta | Qué hay |
|---|---|
| `pluto_run/` | La corrida PLUTO completa: `init.c` (parcheado, ver abajo), `definitions.h`, `pluto.ini` y `corrida.log`. **Las salidas reales** `data.0000..0015.dbl` (t=0 a 30, dbl cada 2.0) están en la máquina local pero NO van al repo de GitHub (245 MB; se regeneran con `make && ./pluto` en ~22 min, ver README raíz) |
| `analisis/` | `teoria_snr.py` (teoría lineal 3D: mapa σ(kx,kz), λ_c, números del caso) y `SNR_3D_analisis.ipynb` (notebook ejecutado con las figuras de la corrida). Figuras en `analisis/figs/` |
| `informe/` | `Informe_SNR_3D.tex` + PDF, y **`Glosario_Conceptos.tex/pdf`** (14 págs): todas las explicaciones de conceptos, términos y ecuaciones del proyecto, con deducciones paso a paso y tabla de símbolos |

## La corrida en una línea por parámetro

- Conf. #09 tal cual el test problem, solo cambié: `USE_RANDOM_PERTURBATION YES` y un
  parche de 5 líneas en `init.c` (había que llamar `RandomSeed()` una vez antes de
  `RandomNumber()`, si no PLUTO aborta con "maximum number of prns exceeded").
- Malla 48×96×48 en x∈[-0.5,0.5], y∈[-1.5,0.5], z∈[-0.5,0.5] (48³·8·... ≈ 221k celdas, 0.23 s/paso en serial, t=30 en ~22 min).
- `ETA=4.0` (ρ₂=4 arriba, eyecta), `GRAV=-0.01`, `CHI=0.2` → B = χ·Bc/√4π ≈ 9.8e-3.
- CFL 0.25, solver `roe`, dbl cada 2.0 (16 salidas).

## Física del mapeo SNR (resumen)

Durante la fase de eyecta dominada, el choque delantero desacelera → en el marco de la
interfaz de contacto hay una gravedad efectiva que apunta hacia el eyecta → el contacto es
Rayleigh–Taylor inestable → dedos de eyecta (los picos que caen en la simulación). Un campo
paralelo al contacto (aquí B = Bx x̂) frena los modos **a lo largo de B** más cortos que
λ_c = χ²Lx ≈ 0.04 y deja crecer los modos transversales (intercambio) → la mezcla queda
elongada a lo largo del campo, como en las sims 3D de Orlando et al. (2012).

## Cómo reproducir

```bash
cd pluto_run
make clean && make && ./pluto        # (~22 min) regenera data.0000..0015.dbl
cd ../analisis
python3 teoria_snr.py                # regenera figs de teoría + sigma_map.csv
jupyter nbconvert --to notebook --execute --inplace SNR_3D_analisis.ipynb
```

(Compilado con `setup.py --auto-update` + `CFLAGS += -D_DEFAULT_SOURCE`, la receta de
siempre; ver `corrida.log` para el log completo.)
