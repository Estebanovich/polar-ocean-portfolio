# Datos del Exp. 2 para el Exp. 3 (emulador del melt)

Fecha: 2026-11-10. Script: `analysis/export_exp3.py` (pruebas en `tests/test_export.py`). Lo corre también `reproduce.sh --quick`.

## Archivo

- **Ruta:** `/Volumes/LaCie/portfolio_runs/exp2/exp3_dataset.npz`, fuera del repositorio (como las corridas).
- **sha256:** `3d2eceb18be2d21c785b175e5a569dd2edbadd0b5437d9aba16ac71ae8d6200f` (1 344 225 bytes).
- **Registro en el repo:** `figures/exp3_dataset.json` (sha256, conteos, corridas excluidas y marcadas, forma y tipo de cada arreglo). Ese JSON es la fuente del README.
- **Reproducible byte a byte:** el .npz se escribe sin compresión, sin pickle y con la fecha de cada miembro del zip fija (1980-01-01). Dos exportaciones del mismo contenido dan el mismo sha256 (prueba `test_write_is_byte_reproducible_and_loads_without_pickle`; comprobado con dos `reproduce.sh --quick` seguidos).
- **Lectura:** `np.load(ruta, allow_pickle=False)`.

## Qué corridas entran

- Los miembros del ensamble de `cases.yaml` (`ensemble.registry_cases`): la referencia `ocean0_4km` y los 28 `ens_*`, todos a 4 km con ΓT_cal,4km = 0.029356 y directorios con sufijo `_r21`. Las corridas post hoc `hol_*` (R24) no son miembros del ensamble y no se exportan.
- **N3:** una corrida que `metrics.load_run` rechaza (no OK, incompleta, con NaN o con la ley de u\* anulada) o que no existe **se excluye**, no se imputa, y queda listada en `excluded_case` / `excluded_reason`. En esta exportación no hay ninguna: 29 corridas, 0 excluidas.
- **N1:** las corridas que fallan N1 se quedan y se marcan (`N1`, `N1_pass`), igual que en H1/H2 (aclaración R22). Son 5: `ens_T-0.5_zc250_s1.0`, `ens_T-0.5_zc550_s1.0`, `ens_T1.5_zc550_s1.0`, `ens_T-0.42_zc402_s0.8`, `ens_T-0.39_zc503_s1.2`.
- **Partición:** `split` = `test` para Tbot = 2.0 (los 3 casos de extrapolación), `train` para el resto (26, incluida la referencia).
- **Verificación cruzada:** N1 y el melt medio de la plataforma de las 29 corridas coinciden exactamente (diferencia máxima 0.0) con `figures/ensemble_table.json`.

## Esquema

n = 29 corridas; malla del modelo a 4 km con la fila de tierra: ny = 21, nx = 120; nz = 36.

| Arreglo | Forma | Tipo | Contenido |
|---|---|---|---|
| `case`, `group`, `split` | (n,) | str | nombre del caso; grupo (`reference`, `factorial`, `geometry`, `extrapolation`); `train` / `test` |
| `Tbot`, `z_c`, `d`, `s` | (n,) | float64 | parámetros del perfil (°C, m, m) y escala del calado; la referencia es la rampa con z_c = 360 m y d = 720 m |
| `N1`, `N1_pass` | (n,) | float64, bool | deriva N1 de pd0 y si es < 3 % |
| `mean_melt_shelf` | (n,) | float64 | melt medio de toda la plataforma, registros 19–24 (m a⁻¹; el ⟨m⟩ de H1) |
| `z` | (nz,) | float64 | altura del centro de cada nivel (RC), m, negativa hacia abajo |
| `T_profile`, `S_profile` | (n, nz) | float64 | perfil inicial (= de restauración) en `z` (°C, PSU), de `make_case._profile` con los parámetros de case.json |
| `x`, `y` | (nx,), (ny,) | float64 | centros de celda (m) |
| `melt` | (n, ny, nx) | float64 | melt medio de los últimos 6 meses (registros 19–24), m a⁻¹ de agua dulce, positivo = fusión; **NaN en celdas secas**; 0 en mar abierto |
| `draft` | (n, ny, nx) | float64 | base del hielo (`icetopo.bin`), m; NaN en celdas secas |
| `wet_mask`, `shelf_mask` | (n, ny, nx) | bool | columna mojada; plataforma = mojada y base < 0 (pd0). Cambian con s |
| `excluded_case`, `excluded_reason` | (k,) | str | corridas excluidas por N3 (k = 0 aquí) |

## Advertencias para el Exp. 3

- Todo es a 4 km. Las limitaciones de resolución del README (§4.4: ψ_bt +125 % frente a 2 km, r = 0.91 de los mapas de melt) se heredan.
- H1 falló (n > 2): la respuesta del melt a Tbot es más empinada que la ley de Holland 2008. El emulador no debe suponer una forma cuadrática.
- Con 26 corridas de entrenamiento, el conjunto es pequeño; la prueba de extrapolación (Tbot = 2.0) queda fuera del rango de entrenamiento en Tbot.
