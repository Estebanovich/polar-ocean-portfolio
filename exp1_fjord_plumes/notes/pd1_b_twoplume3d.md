# PD1 — Configuración (b) `b_twoplume3d`: dos plumas puntuales 3D contra la separación

Task 6, Steps 1–2. Material para el Punto de decisión 1 de (b). **PD1(b) aprobado por Esteban: dt = 2 s y los 14 casos, con animaciones GIF (fallo R29).** Preparación de la campaña en §6–§7; la campaña la lanza el controlador.

## 1. Tabla de parámetros

Todo lo que no aparece aquí es idéntico a (a) y está aprobado en PD1 (`pd1_a_line2d.md`, §1): EOS lineal, icefront, addMass, Sciascia 2013, no deslizamiento, E/S.

| Bloque | Parámetro | Valor | Fuente / razón |
|---|---|---|---|
| Malla | x, y, z | 100 × 120 × 30; dx = dy = dz = 10 m; 1.0 km × 1.2 km × 300 m; f = 0 | `cases.yaml` (b); brief del plan |
| | paredes laterales | filas de tierra j = 0 y Ny−1 (`land_j_rows: [0, -1]`); ancho mojado 1180 m; sin deslizamiento | Es una caja de fjord. Las filas también cierran el periodo en y |
| | última columna (i = 99) | tierra (`land_i_cols: [-1]`) | Igual que (a): cierra el periodo en x en la cara de hielo |
| Descomposición | SIZE.h | sNx = 100, sNy = 24, nPy = 5 (5 procesos MPI en y, los 5 núcleos P) | 120/5 = 24 exacto. Las filas de tierra quedan en los rangos 0 y 4 |
| Tiempo | `deltaT` | **2 s (aprobado en PD1(b) por Esteban)**; 1 s en las corridas cortas | §2 |
| | duración | 3 h por caso (`endTime` 10800 s) | plan, Step 2 |
| Fuentes | Q0 | 10 m³ s⁻¹ por fuente, en la celda mojada del fondo de la columna i = 0 (R8) | `place_sources`. Se usa la separación real `x0_actual` |
| | x0 = 0 | **una sola celda con 2·Q0** (las dos fuentes caen en la misma celda) | `make_case`. Prueba `test_b_x0_0_is_one_cell_of_2Q0` |
| | `_single` | una fuente de Q0 en j = 60: la referencia de M(x0)/(2·M_single) | prueba `test_b_single_is_one_source_of_Q0` |
| | x0 = 320 m | j = 44 y 76; a 435 y 425 m de las paredes | prueba `test_b_widest_separation_far_from_side_walls` |
| Estratificación | dos capas / Sermilik | interfaz (o centro de la picnoclina) a 100 m; capa inferior H2 = 200 m | Fallo R11 |
| Solvers | cg3d / cg2d | 500 iteraciones / 1e-8; 500 / 1e-13 | Igual que (a) |
| OBCS | este | `OB_Ieast = 120*-2` (columna Nx−1), Orlanski + esponja U/V de 10 celdas (100 m), sin `useOBCSbalance` | Igual que (a). El log confirma `OBCS_CHECK: set-up OK` y `check Inside Mask and OB locations: OK`, también con las filas de tierra |
| Salidas | diagnósticos | instantáneas de THETA, SALT, UVEL, VVEL, WVEL y ETAN cada 30 min; `icefrontDiag` promediado cada 15 min | Para `face_total_melt` y C4 |
| | pickups | `pChkptFreq` = 5400 s | |

## 2. Corridas cortas (30 min de modelo, `b_3d_twolayer_x0_40`, np = 5)

Directorios: `/Volumes/LaCie/portfolio_runs/exp1/b_3d_twolayer_x0_40_short30m{,_dt2}/`. Compilación MPI reutilizada (`/Volumes/ESTEBAN_MIT/portfolio_builds/exp1/b_twoplume3d/`).

| intento | dt | pasos | `run_status` | pared (s) | CFL u / v / w máx. | w máx. (m s⁻¹) | cg3d (iteraciones, en 600/1200/1800 s) | NaN | `pgrep` después |
|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 s | 1800 | OK | **192** | 0.044 / 0.042 / **0.128** | 1.28 | 203 / 190 / 199 (tope 500; residuo < 1e-8) | no | vacío |
| 2 | 2 s | 900 | OK | **106** | 0.085 / 0.085 / **0.256** | 1.28 | 237 / 215 / 224 | no | vacío |

- Con dt = 1 s el CFL máximo (0.128) quedó ≤ 0.15. Por eso se probó dt = 2 s, como pide la instrucción.
- La w máxima ya es estacionaria desde los 10 min: 1.267 → 1.276 → 1.277 m s⁻¹.
- Con dt = 2 s el CFL es 0.26, con margen hasta 0.5. dt = 3 s llevaría a ≈ 0.38, con poco margen si la turbulencia crece. Por eso no se propone.
- cg3d se registra solo cada `monitorFreq` (600 s; `solve_for_pressure.F:435-445`), así que hay 3 muestras por corrida: media ≈ 197 y máx. 203 (dt = 1 s); media ≈ 225 y máx. 237 (dt = 2 s). Siempre converge muy por debajo del tope de 500. Son menos iteraciones que en (a) (320–390 en la malla estirada).
- **dt = 1 s contra dt = 2 s a los 30 min:**
  - W: rms de la diferencia 7.8e-4 m s⁻¹ (rms del campo 0.020); máx. puntual 0.088.
  - T: rms 0.004 °C; máx. 0.36 °C en la intrusión (fase de remolinos).
  - Derretimiento total en la cara (ΣICFfwFlx·dy·dz, promedio de 15–30 min): −1344.0 contra −1344.8 kg s⁻¹, una diferencia de **0.06 %**.
- Avisos: ninguno en STDOUT. stderr solo trae `yaksa: 8 leaked handle pool objects` (mpich, inocuo).
- η medio sube 2.1 cm en 30 min, de forma casi lineal, por los 20 m³ s⁻¹ que entran en ~1.2 km² (≈ 3 cm en 30 min si todo quedara atrapado). Hay que vigilarlo en la campaña: con 3 h, si Orlanski no evacua el volumen, η llegaría a ~0.1 m. Es pequeño frente a H, pero conviene revisarlo como en (a).

## 3. Figura rápida y ruido de malla

`pd1_b_quicklook.png` (`pd1_b_quicklook.py`), con w a los 30 min (dt = 1 s), en tres paneles: la columna de la cara (y–z), la fila de una fuente (x–z) y un corte horizontal a −220 m.
- Las dos plumas (j = 58 y 62) suben pegadas a la cara con w ≈ 1.1–1.3 m s⁻¹. Cruzan la interfaz de 100 m y llegan a ~30–50 m de profundidad. La salida se separa de la cara entre ~50 y 150 m de profundidad y avanza unos 400 m. Entre las dos plumas hay una columna más débil (w ≈ 0.45–0.6 m s⁻¹): las plumas ya se tocan a 40 m de separación.
- Cada pluma tiene un núcleo de **1 celda** y una mitad de ancho de ~2 celdas (perfil en y a −150 m: 0.21, 0.62, **1.28**, 0.88, 0.59, …). Está **mal resuelta, pero no es ruidosa**. El campo es simétrico respecto al centro del par con una precisión de 1e-3 m s⁻¹.
- **Indicador de ruido de malla** G = rms(segunda diferencia)/(4 rms w): G = 1 para un tablero de ajedrez de 2Δ y G ≪ 1 para un campo suave.

  | región | G_x | G_y | G_z |
  |---|---|---|---|
  | ambiente (> 200 m de la cara) | 0.033 | 0.014 | 0.100 |
  | todo el dominio mojado | 0.137 | 0.243 | 0.078 |
  | columna de la cara (solo la pluma) | — | 0.268 | 0.035 |

  La fracción de puntos del ambiente en que w cambia de signo en celdas alternas a lo largo de y es 7e-5.
- Los valores altos de G en y vienen del ancho de 1–2 celdas de la pluma (resolución), no de oscilaciones: los perfiles son monótonos a cada lado del núcleo. En el ambiente, w rms = 3.4 mm s⁻¹ y |w| máx. = 0.045 m s⁻¹.
- **Conclusión:** a los 30 min no hay ruido de 2Δ pese a que el Re de malla es de ~20–40. Conviene repetir el chequeo al final de la campaña (3 h, flujo más turbulento).

## 4. Costo y extrapolación (14 casos × 3 h de modelo)

Costo medido por caso de 3 h = 6 × (pared de 30 min). Supone costo constante por paso; la incertidumbre principal es que cg3d suba con la turbulencia (se da también +50 %).

| opción | casos | dt | por caso | total | +50 % |
|---|---|---|---|---|---|
| A. completa, dt = 1 s | 14 | 1 s | 19.2 min | **4.5 h** | 6.7 h |
| **B. completa, dt = 2 s (recomendada)** | 14 | 2 s | 10.6 min | **2.5 h** | 3.7 h |
| C. sin x0 = 20 ni 320 (dt = 2 s) | 10 | 2 s | 10.6 min | 1.8 h | 2.7 h |
| D. solo dos capas (dt = 2 s) | 7 | 2 s | 10.6 min | 1.2 h | 1.9 h |
| E. dz = 15 m (Nr = 20; dt ≈ 3 s, mismo CFL_w), 14 casos | 14 | 3 s | ≈ 4.7 min | ≈ 1.1 h (estimado) | 1.6 h |
| F. caja de 0.6 × 0.8 km (60 × 80 × 30, sNy = 16), dt = 2 s, 14 casos | 14 | 2 s | ≈ 4.2 min | ≈ 1.0 h (estimado) | 1.5 h |

- E y F se escalan con celdas × pasos desde B, sin medirlos.
- E vuelve anisótropa la malla (10 × 10 × 15 m) en una pluma que ya tiene un núcleo de 1 celda.
- En F, x0 = 320 queda a ~230 m de las paredes, y la OB y la esponja se acercan a la salida (la intrusión llega a ~400 m en 30 min). No se recomiendan E ni F, porque el presupuesto alcanza con B.

**Presupuesto del Exp. 1 (24 h):**
- Ya gastado o comprometido: (a) 2.4 h, PD2 1 h + 2.6 h (conv_dx2p5) y (c) 8–12 h, es decir, 14–18 h.
- Con la opción B, 16.5–20.5 h; con +50 % en (b), 17.7–21.7 h. **Cabe.**
- Con la opción A, 18.5–22.5 h (hasta 24.7 h con +50 %). Queda al límite.
- Todas las opciones cumplen ≤ 8 h por corrida.

Disco: ~190 MB por directorio de 30 min en LaCie (exFAT, clústeres grandes); se estiman ≤ 0.6 GB por caso de 3 h. LaCie tiene 3.6 TB libres.

## 5. Puntos para que Esteban decida

1. **dt = 2 s** para los 14 casos (CFL 0.26, cg3d convergido, derretimiento total igual a 0.06 % del de dt = 1 s, 2.5 h en total). ¿Se acepta, o se mantiene dt = 1 s (4.5 h)?
2. Campaña completa de 14 casos (opción B), o una reducida (C/D) si (c) resulta cara.
3. **Resolución de la pluma 3D (clase (c)).** El núcleo es de 1 celda y con x0 = 20 m (2 celdas) las plumas se funden a escala de malla. El extremo de x0 pequeño de M(x0) está limitado numéricamente. Refinar a 5 m costaría ~16× por caso (8× celdas, 2× pasos) y no cabe.
4. **Para C4 (aviso, no se evalúa aún).** A los 30 min la pluma de Q0 = 10 m³ s⁻¹ atraviesa la interfaz de 100 m y llega a ~30–50 m de profundidad. La salida queda entre 50 y 150 m. La altura de despegue debe compararse con la teoría de pluma semicónica (`point_plume(geometry="half_cone")`) en la campaña; si la teoría predice despegue en la interfaz, la discrepancia aparecerá en el Step 5 (PD2).

## 6. Control de validez: retención de volumen (η medio a las 3 h)

Implementado en `analysis/b_checks.py::eta_retention(run_dir)` (pruebas en `tests/test_b_campaign.py`).
- addMass mete el volumen de la descarga dentro de la caja. El deshielo entra como flujo virtual de sal, sin volumen. Con superficie libre lineal, A dη̄/dt = Q_total − Q_OB.
- Si nada saliera por la OB, η̄ seguiría la recta η_atrapado(t) = Q_total t / A. Con A = 1.1682 km² (celdas mojadas), a las 3 h eso da **0.185 m** con dos fuentes y 0.092 m con `_single`.
- Se reportan dos fracciones:
  - `retained_cum` = η̄(3 h)/η_atrapado(3 h). Es informativa, porque incluye el arranque, cuando la salida todavía no llega a la OB.
  - `retained_last_hour` = [η̄(3 h) − η̄(2 h)] / (Q_total · 1 h / A).

**Umbral: se marca la corrida (flag) si `retained_last_hour` > X = 50 %**, es decir, si en la última hora la OB deja pasar menos de la mitad de la descarga. Justificación:
1. Lo que mide es la salud de la OB, no la física de la pluma. Aun con retención total, η̄ = 0.19 m es el 0.06 % de H. La desviación espacial de η es de 0.25 mm (casi uniforme), así que su efecto sobre la pluma es nulo. La señal es muy sensible: una velocidad media en la OB de solo 6e-5 m/s (20 m³/s en 1180 m × 300 m) decide el balance.
2. En (a) la retención acumulada fue 32 % a las 2 h y bajaba, con la OB funcionando.
3. Un factor 2 (la OB evacua menos de la mitad cuando ya debería estar en régimen) es una falla clara de la OB que no depende de ningún caso en particular. Por eso no se ajusta a una sola corrida.

Un caso marcado **no se descarta automáticamente**: se revisan su serie de η̄ (monitor cada 600 s, `dynstat_eta_mean`) y el transporte en la OB, y se compara su derretimiento con los casos vecinos. Cualquier diferencia se clasifica, probablemente como (c), elección numérica de la OB.

**Aviso de las corridas cortas (30 min).** Con dt = 1 s, η̄ = 0.96, 1.70 y 2.14 cm a los 10, 20 y 30 min frente a 1.03, 2.05 y 3.08 cm atrapados: retención acumulada del 69 %, y la OB empieza a evacuar. Con **dt = 2 s**, η̄ = 0.99, 1.90 y **3.40 cm**: 110 % a los 30 min. Entre los 20 y 30 min la OB *metió* agua (≈ +0.5 cm más que la descarga). El derretimiento total no cambió (0.06 %), y la salida baroclínica aún no llega a la OB (~400 m a los 30 min). Es transitorio de la OB de Orlanski, probablemente por su dependencia de dt (la velocidad de fase se recorta a CMAX·dx/dt, `CMAX = 0.45`). No se ha investigado más (no se permiten corridas mientras corre conv_dx2p5). **Sugerencia para el controlador:** al terminar el primer caso de la campaña (~11 min), correr `eta_retention` sobre él antes de dejar seguir la serie, o vigilar `dynstat_eta_mean` en `STDOUT.0000`. Si `retained_last_hour` > 50 % de forma sistemática, avisar a Esteban: el arreglo sería de la OB (por ejemplo, prescribir el transporte barotrópico de salida), no de la física.

## 7. Preparación de la campaña (R29)

- `cases.yaml` (b): `deltaT: 2.0` → `nTimeSteps` = 5400 por caso de 3 h. Se volvió a generar `b_3d_sermilik_x0_160` para comprobarlo: `deltaT=2.0`, `nTimeSteps=5400`, SIZE.h OK.
- **Diagnósticos de animación** (en `data.diagnostics`, listas 8 y 9; las listas 1–7 no cambian, lo verifica una prueba). Son solo salida, la física no cambia, y `code/` tampoco cambia, así que el binario es el mismo.
  - `animW`: WVEL 3D cada 120 s. pkg/diagnostics solo permite elegir niveles verticales (`levels(:,n)`, `diagnostics_readparms.F:49,326-331`), no un corte horizontal, y la sección y–z de la cara necesita todos los niveles. Son 2.88 MB por instantánea.
  - `animT`: THETA solo en los niveles 8 y 16 (z = −75 y −155 m), 0.19 MB por instantánea.
- **Disco extra:** 90 instantáneas por caso. `animW` ocupa 259 MB y `animT` 17 MB, es decir, **≈ 0.28 GB por caso de tamaño aparente y ≈ 0.35 GB en LaCie** (exFAT con clústeres de 256 KB, más los .meta). **Total para 14 casos: ≈ 3.9 GB aparentes, ≈ 4.8 GB en disco.** Con la salida previa, un caso ocupa ~0.7–0.8 GB y la campaña ~11 GB, sobre 3.6 TB libres.
- **`analysis/animate_b.py`** genera un GIF por caso en `figures/anim/<case>.gif`. Sigue las convenciones de `BayIW_Rectan/make_animation.py`: escalas fijas para todos los casos, cmocean (`pfcommon.style`), tierra gris y reloj.
  - Paneles: (a) w en la cara (y–z, fuentes marcadas); (b) T − T_init a −75 m; (c) T − T_init a −155 m; (d) barra de tiempo.
  - Se usa la anomalía porque a −155 m la pluma solo enfría ~0.2 °C, invisible en una escala de −1.6 a 4.1 °C.
  - Escalas fijas: w ±1.2 m/s; ΔT ±4 °C en (b) y ±0.2 °C en (c).
  - Probado con la corrida corta, que solo tiene 2 cuadros (0 y 30 min): ~25 KB por cuadro a 80 dpi, así que 90 cuadros dan ~2–3 MB, por debajo del límite de 8 MB. Si se pasa, el script duplica el paso entre cuadros. MP4 opcional con `--mp4` (ffmpeg).
- **Comando de la campaña** (lo lanza el controlador cuando termine conv_dx2p5 y `pgrep` esté vacío):

  ```
  cd /Users/estebancruz/Phd_project/portfolio && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /opt/anaconda3/envs/pyocean/bin/python exp1_fjord_plumes/scripts/run_campaign.py --config b_twoplume3d --cases 'b_3d_*' --np 5 > /Volumes/LaCie/portfolio_runs/exp1/campaign_b.log 2>&1
  ```

  Esperado: 14 casos en serie, ~10.6 min cada uno (≈ 2.5 h), todos `OK`; sale con código 1 si alguno no lo está.
- Después, por caso:
  - `python -c "from exp1_fjord_plumes.analysis.b_checks import eta_retention; print(eta_retention('<run_dir>'))"`
  - `python exp1_fjord_plumes/analysis/animate_b.py /Volumes/LaCie/portfolio_runs/exp1/<case>`
