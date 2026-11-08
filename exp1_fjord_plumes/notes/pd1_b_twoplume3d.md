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

**Umbral: se marca la corrida (flag) si |`retained_last_hour`| > X = 50 %**, es decir, si en la última hora la OB deja pasar menos de la mitad de la descarga (> +0.5) o drena más de 1.5 veces la descarga (< −0.5; corrección de la revisión). Justificación:
1. Lo que mide es la salud de la OB, no la física de la pluma. Aun con retención total, η̄ = 0.19 m es el 0.06 % de H. La desviación espacial de η es de 0.25 mm (casi uniforme), así que su efecto sobre la pluma es nulo. La señal es muy sensible: una velocidad media en la OB de solo 6e-5 m/s (20 m³/s en 1180 m × 300 m) decide el balance.
2. En (a) la retención acumulada fue 32 % a las 2 h y bajaba, con la OB funcionando.
3. Un factor 2 (la OB evacua menos de la mitad cuando ya debería estar en régimen) es una falla clara de la OB que no depende de ningún caso en particular. Por eso no se ajusta a una sola corrida.

Un caso marcado **no se descarta automáticamente**: se revisan su serie de η̄ (monitor cada 600 s, `dynstat_eta_mean`) y el transporte en la OB, y se compara su derretimiento con los casos vecinos. Cualquier diferencia se clasifica, probablemente como (c), elección numérica de la OB.

**Aviso de las corridas cortas (30 min).** Con dt = 1 s, η̄ = 0.96, 1.70 y 2.14 cm a los 10, 20 y 30 min frente a 1.03, 2.05 y 3.08 cm atrapados: retención acumulada del 69 %, y la OB empieza a evacuar. Con **dt = 2 s**, η̄ = 0.99, 1.90 y **3.40 cm**: 110 % a los 30 min. Entre los 20 y 30 min la OB *metió* agua (≈ +0.5 cm más que la descarga). El derretimiento total no cambió (0.06 %), y la salida baroclínica aún no llega a la OB (~400 m a los 30 min). **Hipótesis (sin verificar hasta leer el gate del primer caso):** es un transitorio de la OB de Orlanski que depende de dt (la velocidad de fase se recorta a CMAX·dx/dt, `CMAX = 0.45`). No se ha investigado más (no se permiten corridas mientras corre conv_dx2p5). Por eso la campaña lleva un gate después del primer caso (§7, fallo R30). Si el gate falla en (i), se avisa a Esteban: el arreglo sería de la OB (por ejemplo, prescribir el transporte barotrópico de salida), no de la física.

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
- **Procedimiento de la campaña en dos pasos (fallo R30).** Lo lanza el controlador cuando termine conv_dx2p5 y `pgrep -fl "mpirun|mitgcmuv"` esté vacío.

  **Paso 1: solo el primer caso** (~11 min):
  ```
  cd /Users/estebancruz/Phd_project/portfolio && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /opt/anaconda3/envs/pyocean/bin/python exp1_fjord_plumes/scripts/run_campaign.py --config b_twoplume3d --cases b_3d_twolayer_x0_40 --np 5 > /Volumes/LaCie/portfolio_runs/exp1/campaign_b_first.log 2>&1
  ```

  **Paso 2: gate** (sin MITgcm, segundos):
  ```
  /opt/anaconda3/envs/pyocean/bin/python exp1_fjord_plumes/scripts/b_gate.py /Volumes/LaCie/portfolio_runs/exp1/b_3d_twolayer_x0_40
  /opt/anaconda3/envs/pyocean/bin/python exp1_fjord_plumes/analysis/animate_b.py /Volumes/LaCie/portfolio_runs/exp1/b_3d_twolayer_x0_40 --check-scales
  ```

  **Criterio de PASS** (`b_gate.py` sale con 0):
  - (i) |`retained_last_hour`| ≤ 0.5 a las 3 h (§6);
  - (ii) en la iteración 900 (30 min), `icefrontDiag` y las instantáneas THETA, SALT, UVEL, VVEL y ETAN del caso son **idénticas byte a byte** a `icefrontDiag` y a los volcados T, S, U, V y Eta de la corrida corta de referencia `b_3d_twolayer_x0_40_short30m_dt2` (mismas entradas y binario, hecha antes de R29). Eso prueba que los diagnósticos de animación no cambiaron la física.

  Cualquier diferencia o archivo faltante es FAIL, con la diferencia máxima impresa. Si hay FAIL: **parar** y avisar a Esteban.

  `--check-scales` no decide PASS/FAIL. Si un panel tiene más de ~1 % de píxeles saturados en promedio, conviene ampliar `W_LIM`/`DT_LIM` (igual para todos los casos) antes de hacer los GIF.

  **Paso 3, solo si hay PASS: los 13 casos restantes** (el primero se salta porque ya está OK, igual que la compilación si está al día):
  ```
  cd /Users/estebancruz/Phd_project/portfolio && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 /opt/anaconda3/envs/pyocean/bin/python exp1_fjord_plumes/scripts/run_campaign.py --config b_twoplume3d --cases 'b_3d_*' --np 5 --skip-ok > /Volumes/LaCie/portfolio_runs/exp1/campaign_b.log 2>&1
  ```
  Esperado: 13 casos × ~10.6 min ≈ 2.3 h, todos `OK`.

- Después, por caso:
  - `python -c "from exp1_fjord_plumes.analysis.b_checks import eta_retention; print(eta_retention('<run_dir>'))"`
  - `python exp1_fjord_plumes/analysis/animate_b.py /Volumes/LaCie/portfolio_runs/exp1/<case>`

## Cambio de frontera a Stevens y expectativas de la puerta, fijadas antes de correr (2026-11-05)

**Resultado de la puerta con Orlanski:** FALLÓ. En `b_3d_twolayer_x0_40` (3 h) η medio llegó a 2.92 m, contra 0.18 m que explicaría la descarga, y entraban ~800 m³/s por la frontera desde 1–1.5 h. Orlanski deja libre la velocidad normal media en profundidad, y en 3D eso se volvió inestable.

**Cambio, aprobado por Esteban:** frontera de Stevens (1990), con u normal media en profundidad prescrita, u_ob = Q0_total / A_ob = 20 / 354000 = 5.65e-5 m/s (commits 5e76391 y siguientes). Se descartó Flather: con c = √(gH) ≈ 54 m/s, dx = 10 m y dt = 2 s, el número de Courant sería ≈ 11, inestable en forma explícita. MITgcm no lo trae.

**Puerta:** se vuelve a correr el mismo caso, 3 h con np = 5.
- **Criterio sin cambios:** |retained_last_hour| ≤ 0.5.
- **Expectativa física:** η medio ≈ 0. pkg/icefront no agrega volumen (solo modifica T y S), y Stevens saca exactamente Q0. Es decir, |η medio| ≪ 0.18 m.
- **Chequeo informativo**, que no decide la puerta: `obc_E_uVel_Int` en STDOUT ≈ +20 m³/s.
- **Antes de las 3 h:** prueba corta de 30 min. En STDOUT deben aparecer `ALLOW_OBCS_STEVENS`, `OBCS_FIELDS_LOAD: Reading initial data` y `useStevensEast = T`.

### Resultado de la puerta con Stevens (2026-11-05)

- **(i) Balance de volumen: PASA.**
  - η medio = 7.6e-15 m;
  - retained_last_hour = 9.6e-14, contra un límite de ≤ 0.5;
  - `obc_E_uVel_Int` = 20.000 m³/s.
- **(ii) Identidad bit a bit frente a la referencia de 30 min: "FAIL", pero por construcción, no por un defecto.** La referencia (`b_3d_twolayer_x0_40_short30m_dt2`) se corrió con Orlanski. Al cambiar el esquema de frontera, los campos no pueden ser idénticos.
- **Por qué igual se da por cumplido el fin de (ii):** (ii) verifica que la salida extra de animación (R29) no altera la física. Eso ya se demostró con el mismo código en el primer intento de la puerta (Orlanski, `b_3d_twolayer_x0_40_orlanski`): pasó (ii) y solo falló (i). La salida de diagnósticos no depende del esquema de frontera.
- **Ruling R32:** la puerta se da por PASADA. Se lanzan los 13 casos restantes.

## Definiciones operativas de C4 y C4′, fijadas antes de analizar (2026-11-05, 22:50)

Las 14 corridas (b) terminaron OK a las 22:33. Ninguna métrica de derretimiento ni de despegue se ha calculado todavía. Estas definiciones se fijan antes de calcular y no se cambian después de ver los resultados.

- **Ventana:** de 2 a 3 h, la última hora de cada corrida de 3 h. Como dato informativo, también se reporta la ventana de 1 a 2 h y el cambio relativo entre las dos, como indicador de estacionariedad. Ese indicador no entra en el veredicto.
- **`face_total_melt(run_dir, t_window_h)`:** es la suma, sobre las celdas de la cara (i = 0), de m_hielo · dy · dz · hFacC. Se promedia en la ventana con los promedios horarios de ICFfwFlx, con la convención de `plume_diag`. La unidad es m³ s⁻¹ de hielo.
- **Normalización:** M(x0) / (2 · M_single), con M_single tomado del caso `_single` de la misma estratificación. Se grafica contra x0_actual / x0c, con x0c = 57.14 m tomado de case.json.
  - En `sermilik`, x0c es nominal: no hay una interfaz de dos capas y se usa el mismo H2.
  - x0 = 0 es una sola celda con 2 · Q0, es decir, el límite de fuentes fusionadas.
- **C4 (±0.15 H = 45 m):** se usa la pluma única (`_single`).
  - Despegue del modelo: el mismo `detach_depth` (con la misma DETACH_FRAC que en (a)) aplicado al perfil de w de la columna de la cara. Ese perfil es el máximo en y del promedio de w en la ventana.
  - Despegue teórico: `point_plume(geometry="half_cone")` con Q0 = 10 m³ s⁻¹ y el ambiente con la misma convención que C2 en (a), es decir, T/S medios de la ventana lejos de la cara.
  - El despegue teórico con el perfil inicial se reporta como post-hoc.
  - Se evalúa en las dos estratificaciones. **C4 pasa solo si pasan ambas.**
- **C4′ ([0.4, 0.8]):** M(x0 = 0) / (2 · M_single) para cada estratificación. **Pasa solo si ambas están en el rango.**
- **Referencia de laboratorio:** los puntos de Cenedese & Gatto (2016a, Fig. 2) se digitalizan a mano desde el PDF de `literatura/` y se rotulan en la figura como "hand-digitised, approximate".
- **Modelo de plumas fusionadas** (opcional en el plan): no se implementa (regla R34).
- Si C4 o C4′ fallan, se pasa a **PD2**: se clasifica la diferencia como (a)–(d) y se presenta a Esteban. El criterio no se redefine.

## Hallazgo crítico: la frontera de Stevens genera velocidades espurias (2026-11-05, 23:15)

El análisis de la Task 6 (commit 2d86679) encontró que **4 corridas sermilik divergieron** aunque `run_status` dice OK:
- son x0_40, x0_80, x0_160 y x0_320;
- cg2d da NaN a partir de t ≈ 2.6–2.7 h;
- los volcados finales son NaN.

`require_ok` no lo detecta, porque solo revisa el estado y `last_iter`.

**Verificación del controlador.** Se tomó max|u| en las 3 columnas junto a la OB (i = 96–98) y en el interior (i = 10–89), en m/s:

| caso | 0.5 h | 1 h | 1.5 h | 2 h | 2.5 h |
|---|---|---|---|---|---|
| twolayer x0_0 | 0.01 / 0.26 | 0.01 / 0.25 | 0.10 / 0.26 | **1.24** / 0.28 | |
| twolayer x0_40 | 0.02 / 0.33 | 0.29 / 0.32 | **2.11** / 0.65 | 1.42 / 0.43 | |
| twolayer x0_320 | 0.01 / 0.24 | 0.11 / 0.23 | **1.75** / 0.41 | 1.50 / 0.45 | |
| twolayer single | 0.01 / 0.23 | 0.02 / 0.20 | 0.55 / 0.21 | **2.13** / 0.62 | |
| sermilik x0_0 | 0.02 / 0.24 | 0.01 / 0.24 | 0.02 / 0.28 | 0.05 / 0.25 | |
| sermilik x0_20 | 0.02 / 0.35 | 0.02 / 0.30 | 0.04 / 0.36 | 0.14 / 0.31 | |
| sermilik single | 0.01 / 0.08 | 0.02 / 0.08 | 0.04 / 0.08 | 0.30 / 0.08 | |
| sermilik x0_40 | 0.02 / 0.28 | 0.07 / 0.26 | 0.62 / 0.26 | 1.49 / 0.39 | NaN a las 2.67 h |

(Las horas están en los mismos tiempos que los volcados UVEL; "/" separa la OB del interior.)

- La inestabilidad empieza en la columna u de la OB (i = 98), a k ≈ 9–10 (~100 m), es decir, en la capa de salida de la pluma.
- Crece cuando la intrusión llega a la frontera, a ~1.5–2 h, y se propaga hacia el interior.
- **Afecta a casi todas las corridas (b)**, no solo a las que dieron NaN.
- **(c) muestra la misma firma, más débil.** En la prueba de 1 día, a las 24 h, max|u| es 0.55 m/s en las columnas 196–198 frente a 0.10 m/s en el interior. No hay NaN.
- El gate de (b) (R32) y la prueba de (c) verificaron η y el transporte de la OB, que son correctos, pero no las velocidades junto a la frontera. Ese es un hueco del diseño del gate.

**Consecuencia.**
- La ventana pre-registrada de 2–3 h cae en corridas contaminadas. **Los veredictos C4 (FAIL) y C4′ (PASS) de 2d86679 quedan SUSPENDIDOS**: no se reportan como resultado. Tampoco se cambia la ventana post-hoc para salvarlos.
- Hay que arreglar la frontera y volver a correr. Esto es una elección numérica nuestra (clase (a)/(c)) y no física.

## Diagnóstico de la OB: variantes y criterios, fijados antes de correr (2026-11-06; plan aprobado por Esteban)

**Caso de prueba:** `b_3d_twolayer_x0_40`, el peor caso: 2.11 m/s junto a la OB a 1.5 h. Se corre 3 h con np = 5, igual que en la campaña. Las variantes se prueban en orden y se paran en la primera que pasa:
- **V1:** `useStevensPhaseVel=.FALSE.` (`--set stevens_phasevel=0`). Stevens queda solo con advección (obcs_calc_stevens.F:141). Sufijo `_stv_nophase`.
- **V2:** zona RBCS de T/S en (b), como en (c). Solo se prueba si V1 falla.

**Criterios (todos deben pasar):**
1. `require_ok` pasa. Desde c5561f8 rechaza NaN o Infinity en el log.
2. En cada volcado de UVEL (cada 0.5 h): max|u| en las columnas i = 96–98 ≤ 1.5 × max|u| en el interior, i = 10–89.
3. |η medio| < 1e-3 m y `obc_E_uVel_Int` = 20 m³/s, con error relativo < 1e-6, en todas las salidas del monitor.

**Si pasa:**
- Se aplica a (c) y se verifica con `c_fjord_shelf` a 30 d. Los criterios son 1 y 3 con Q0 = 0 (transporte 0), y el criterio 2 en las columnas 196–198 frente a 10–189, en cada volcado de U.
- Después se vuelven a correr los 14 casos de (b) y los 4 de (c).

### Resultado de V1 y nueva evidencia (2026-11-06)

**V1 (`useStevensPhaseVel=.FALSE.`) FALLA el criterio 2.** La razón OB/interior llega a 3.51 (sin V1 era 3.31). Los criterios 1 y 3 pasan. La velocidad de fase no es la causa.

**Evidencia nueva:**
- **La corrida Orlanski del mismo caso tiene una OB limpia.** Es `b_3d_twolayer_x0_40_orlanski`, con la esponja U/V y sin esponja T/S, y su razón máxima es 0.41. Su único problema era la deriva de η. Como tampoco tenía relajación de T/S, la falta de esa relajación no explica la celda espuria, y V2 pierde el sustento.
- **En el caso de referencia, a las 2 h, hay una celda de vuelco espuria pegada a la OB:** u = −0.55 m/s en k = 6–8 y +0.55 m/s en k = 12, con una anomalía de T de +1.3 °C en k = 10, en i = 98.
- **Mecanismo candidato (fuente MITgcm).**
  - Con Stevens, la parte baroclínica de la u de la OB es la propia velocidad de la OB en el paso anterior: `OBEuStevens`, guardada en `obcs_save_uv_n.F:99` y usada en `obcs_calc_stevens.F:200-219`. Es decir, se integra con su propia ecuación de momento y no se radía.
  - La esponja U/V relaja hacia `OBEu` (`obcs_sponge.F:141-147`), que bajo Stevens es esa velocidad autoactualizada y no la prescrita.
  - Resultado: la esponja arrastra el interior hacia el perfil espurio de la frontera en lugar de amortiguarlo. Es una posible retroalimentación, la misma clase de incompatibilidad que la de la esponja T/S (R33).

**Nuevo orden de variantes (regla R36, dentro del plan de diagnóstico aprobado; fijado antes de correr). Los criterios 1–3 no cambian:**
- **V4:** esponja U/V apagada (`--set ob_uv_sponge=0`), sin recompilar. Si pasa, la esponja era el amplificador.
- **V3:** `#define OBCS_STEVENS_USE_INTERIOR_VELOCITY` (`obcs_calc_stevens.F:39-42`, copia la u interior). Se probaría en un build aparte con la configuración de (b).
- **V5:** reemplazar la esponja OBCS U/V por una relajación RBCS de u/v hacia la u prescrita.
- **V2:** queda al final.

### Resultado de V4, y criterios de V3 y V6 (fijados antes de correr, 2026-11-06)

**V4 (esponja U/V apagada) FALLA el criterio 2.** La razón llega a 3.09: la velocidad junto a la OB es de 1.18 m/s, frente a 2.1 m/s con la esponja. La esponja amplifica la celda espuria, pero no la causa. La causa está en la persistencia baroclínica de Stevens: la u baroclínica de la OB es su propia velocidad del paso anterior.

**Literatura** (revisión de la biblioteca de la tesis; informe de trabajo local, no incluido en el repositorio):
- La opción mejor respaldada es fijar el transporte medio con una restricción global de volumen y radiar el resto (Marchesiello et al. 2001, §4.2, ecs. 16–20; Marsaleix et al. 2006; Palma & Matano 1998).
- Marchesiello et al. (2001, §3.3 y §6) reportan "sin éxito" el uso de la física reducida tipo Stevens en una malla C.
- La OB de la tesis BayIW era Orlanski con CMAX 0.45 y cvelTimeScale 250 s, una esponja U/V de 1200 → 60 s sin esponja T/S, y `useOBCSbalance`. Archivo: `/Volumes/ESTEBAN_MIT/MITgcm/verification/BayIW_Rectan/run_expand/data.obcs`.

**Builds de diagnóstico.** Son las configuraciones `b_obdiag` y `b_obdiag_orl`, con los mismos binarios de entrada que la producción según `test_obdiag_cases_reuse_the_b_generator`.
- **V3:** `b_obdiag_iv_twolayer_x0_40`. Es Stevens con `OBCS_STEVENS_USE_INTERIOR_VELOCITY`, es decir, la u baroclínica de la OB se copia de la columna interior. Usa los mismos criterios 1–3.
- **V6:** `b_obdiag_orl_twolayer_x0_40`. Es Orlanski más `useOBCSbalance`, el esquema de BayIW, y equivale a la plantilla anterior a Stevens más el balance.
  - El balance no cuenta `addMass` (`obcs_balance_flow.F:304-345`), así que el transporte neto por la OB es ≈ 0 y el dominio se llena a una tasa Q0·t/A.
  - Por eso sus criterios son:
    1. `require_ok` pasa.
    2. Igual que el criterio 2 anterior (razón ≤ 1.5).
    3. |`obc_E_uVel_Int`| ≤ 0.2 m³/s (1 % de Q0) en todas las salidas, y |η̄(t) − Q0·t/A| ≤ 0.02 m en todas las salidas, con Q0 = 20 m³/s y A = área mojada de hFacC en superficie.
  - **Adoptar V6 cambia la física:** el agua de la pluma no sale en neto. En (b) es despreciable (η sube 0.19 m en 3 h sobre 300 m). En (c) serían ~1.2 m en 30 d y un sesgo de 100 m³/s en el intercambio. Adoptarlo en (c) requiere la aprobación de Esteban.

### Resultados de V3 y V6 (2026-11-06)

| variante | C1 | C2 (razón máxima) | C3 |
|---|---|---|---|
| V3: Stevens + u baroclínica interior | PASA | **0.50**, PASA | \|η̄\| 1.7e-16 m; error de transporte 1e-12. PASA |
| V6: Orlanski + balance (BayIW) | PASA | **0.46**, PASA | \|η̄ − Q0·t/A\| máx 1.9e-3 m (A = 1.1682 km²); \|transporte\| máx 1.8e-12. PASA |

(El evaluador genérico marcó C3 de V6 como FALLA porque supone η = 0. El criterio pre-registrado de V6 es η = Q0·t/A, y ese pasa.)

- **Diagnóstico:** la celda espuria venía de la persistencia baroclínica de Stevens (`OBEuStevens`). Copiar la u interior (V3) o radiar con Orlanski (V6) la eliminan.
- **Elección (regla R37): V3.** Conserva la física: el agua de la pluma sale exactamente a Q0. Además, se queda dentro de la decisión ya aprobada de usar Stevens. V6 hace que se acumule agua.
  - El encabezado de MITgcm advierte que V3 "can give unexpected results" (`obcs_calc_stevens.F:39-42`). Por eso, antes de adoptarlo en las campañas, se verifica en (c) con `c_fjord_shelf` a 30 d, con los criterios pre-registrados.
- V3 sigue escribiendo `pickup_stevensE`, así que el reinicio de dos etapas de (c) no cambia.

### Re-corrida de las 14 corridas (b) con V3 (2026-11-06)

Todas terminan OK y **todas pasan los criterios 1–3 de la OB**. La razón OB/interior máxima por caso es:

| | x0_0 | x0_20 | x0_40 | x0_80 | x0_160 | x0_320 | single |
|---|---|---|---|---|---|---|---|
| twolayer | 0.28 | 0.36 | 0.50 | 0.76 | 0.88 | 0.84 | 0.52 |
| sermilik | 0.10 | 0.11 | 0.44 | 0.32 | 0.57 | 0.66 | 0.36 |

Las 4 sermilik que divergían ahora son finitas. Las corridas anteriores quedan en `*_stv0`. C4 y C4′ se evalúan con las definiciones pre-registradas sobre estas corridas.

## PD2(b): resultado de C4 y C4′ sobre las corridas V3 (2026-11-06)

Los números vienen de la re-ejecución del revisor de la Task 6 sobre las 14 corridas V3, todas finitas. La ventana es la pre-registrada, 2–3 h. El ambiente teórico sigue la redacción vinculante de (b): media de la ventana lejos de la cara (regla R39). La lectura alternativa con el perfil inicial da el mismo veredicto.

| criterio | twolayer | sermilik | veredicto |
|---|---|---|---|
| **C4** (\|Δz_despegue\|/H ≤ 0.15) | 0.214 (modelo −1.9 m; teoría −66.2 m) | 0.187 (−1.5 m; −57.5 m) | **FALLA** |
| **C4′** (M(0)/(2·M_single) ∈ [0.4, 0.8]) | 0.520 | 0.525 | **PASA** |

**C4′.** La interacción de fuentes juntas reduce el derretimiento total a ~0.52 de dos plumas independientes, en las dos estratificaciones. Coincide con el ≈ 0.5 del laboratorio (Cenedese & Gatto 2016a). La escala de Slater/Cowton daría ≈ 0.66; esa cifra se usa solo como verificación de plausibilidad.

**C4 FALLA: la pluma del modelo llega a la superficie y la teoría semicónica despega a ~60 m.** Se revisó el perfil de w en la cara (máximo en y, media de la ventana) de `_single`:
- w ≈ 0.45 m/s todavía a 10 m de profundidad, con un máximo de ~1.3 m/s a 150 m.
- T en la cara ≈ 3.8 °C en toda la columna, contra −1.5 °C en el ambiente superior. El agua profunda sube casi sin diluirse.
- El modelo arrastra mucho menos agua ambiente que la teoría.

**Clasificación (hipótesis, sin verificar):**
- **(b)/(d), física de pluma de pared frente a pluma semicónica.**
  - Kimura et al. (2014, resumen, l. 19–20): "The presence of an ice face reduces the entrainment of seawater into the meltwater plumes, so the plumes remain attached to the ice front, in contrast to previous simple models".
  - Ezhova et al. (2018, l. 29): el flujo de volumen de una pluma de pared es la mitad del de una pluma cónica con el doble de flujo de flotabilidad.
  - La teoría semicónica puede sobreestimar el arrastre de una pluma pegada a la pared.
- **(c), resolución y cierre.**
  - A dx = 10 m, la pluma mide pocas celdas de ancho (b ≈ 6α z/5 ≈ 24 m a 200 m sobre la fuente), así que el arrastre lo fija sobre todo el cierre de submalla.
  - En (a), la PD2 mostró una sensibilidad fuerte al cierre (`sens_xu2012`: R_Q de 0.547 a 0.720).
- No es error nuestro de análisis (clase a): el revisor verificó el estimador. La profundidad de despegue del modelo está acotada en ~−2 m por w = 0 en la superficie.

**Pruebas propuestas (para decidir con Esteban; se pre-registran antes de correr):**
- **D1 (diagnóstico post-hoc, barato, sin correr el modelo):** comparar el flujo de volumen Q(z) de la pluma del modelo junto a la cara con el Q(z) teórico. Si el modelo arrastra mucho menos, eso confirma el mecanismo, pero no distingue entre (b)/(d) y (c).
- **D2 (convergencia):** `b_3d_twolayer_single` a dx = dy = dz = 5 m. Son ~8× celdas y ~2× pasos, unas 3 h. Si Δz_despegue no se mueve con la resolución, se descarta que (c) sea la causa dominante.

C4 se reporta como FALLA bajo su definición. No se redefine.

### Corrección: C4′ pasa por construcción y no informa sobre la interacción (2026-11-06)

C4′ **PASA bajo su definición pre-registrada** (0.520 y 0.525), y eso no se cambia. Pero el pase **no es evidencia de interacción entre plumas**:
- M_single vale 1.38 y 1.36 m³/s (twolayer y sermilik), mientras que M con dos fuentes vale 1.43–1.48 m³/s. Agregar una segunda pluma sube el derretimiento total de la cara solo un 4–8 %.
- La razón M(x0)/(2·M_single) vale ≈ 0.52–0.54 **para todas las separaciones**, incluso x0 = 320 m = 5.6·x0c. Es casi 0.5 por construcción.
- Perfil del derretimiento a lo largo de la cara (columna integrada en z, ventana 2–3 h): **~90 % del derretimiento total es de fondo** (derretimiento ambiente sobre los 1180 m de cara). La pluma única agrega solo ~0.10–0.12 m³/s sobre ~1.3 m³/s de fondo.
- **Falla de diseño de la métrica (clase a, error nuestro).** El total sobre una cara ancha mide sobre todo el derretimiento ambiente. En el laboratorio de Cenedese & Gatto la cara era estrecha y el derretimiento por pluma dominaba. La comparación con su ≈ 0.5 no es válida.

**Métrica post-hoc, etiquetada como tal (no es veredicto).** Se definió el exceso de la pluma como el total menos el fondo, con fondo = mediana de la columna × número de columnas mojadas. Este estimador se eligió después de ver los datos.

| x0 | 0 | 20 | 40 | 80 | 160 | 320 |
|---|---|---|---|---|---|---|
| exceso / (2·exceso_single), twolayer | 0.75 | 0.90 | 0.91 | 0.92 | 0.94 | 0.85 |
| exceso / (2·exceso_single), sermilik | 0.78 | 0.90 | 0.89 | 0.91 | 0.94 | 0.83 |

- Con fuentes juntas, cada pluma derrite ~25 % menos que una pluma independiente, y el valor se recupera con la separación hasta ~0.94. Es cualitativamente como en el laboratorio, pero mucho más débil que su ≈ 0.5.
- La caída en x0 = 320 m probablemente es un artefacto del estimador de fondo, porque las plumas ensanchadas contaminan la mediana. Está sin verificar.
- Es sensible al método de fondo. Antes de usarla en el README hay que definir un estimador robusto, por ejemplo una corrida sin descarga para el fondo, y fijarlo por escrito.

### Métrica de exceso de la pluma con fondo limpio: definición fijada antes de calcular (2026-11-06, aprobada por Esteban)

**Casos de fondo.** `b_bg_twolayer_nodischarge` y `b_bg_sermilik_nodischarge`: (b) con Q0 = 0, el mismo ambiente, la misma OB y 3 h. Las entradas son idénticas a las de `_single` salvo la descarga, según `test_b_background_case_has_no_discharge_and_same_ambient`.

**Definición:**
- E(caso) = M(caso) − M(fondo de la misma estratificación). M es `face_total_melt` en la ventana de 2–3 h.
- R_E(x0) = E(x0) / (2·E_single), graficada contra x0_actual/x0c.
- También se reportan M_fondo/M_single y el ruido del fondo, ε = |M_fondo(2–3 h) − M_fondo(1–2 h)|.
- Si E_single < 10·ε, R_E se declara no resoluble.

**Estatus.** Es una métrica **posterior a C4′, informativa**. No reemplaza a C4′, que sigue PASA bajo su definición, con la advertencia anterior. No tiene umbral de pase. Se compara cualitativamente con el ≈ 0.5 del laboratorio, sin forzar el acuerdo.

**Supuesto declarado.** El fondo es aditivo: la pluma no cambia el derretimiento ambiente fuera de su huella. Si la pluma altera la circulación de toda la cara, E mezcla los dos efectos.

### D1 y D2 para el fallo de C4: fijados antes de calcular o correr (2026-11-06, aprobados por Esteban)

**D1 (post-hoc, sin modelo).** Q_m(z) es el flujo de volumen ascendente junto a la cara: ∫ w⁺ dx dy sobre i < 10 (100 m), con la media de w en la ventana de 2–3 h, en `b_3d_{twolayer,sermilik}_single`. Se compara con Q_t(z) de `point_plume(half_cone)` con el mismo ambiente que C4.
- **Lectura:** si Q_m/Q_t < 0.6 en la mitad inferior de la columna, el modelo arrastra mucho menos que la teoría, y eso explica que la pluma llegue a la superficie.
- D1 no separa la causa física (b/d) de la numérica (c).

**D2 (convergencia).** `b_dx5_twolayer_single`: dx = dy = dz = 5 m, nx × ny × nz = 200 × 240 × 60, dt = 1 s, 3 h. La física queda fija: viscosidades y difusividades iguales, mismo perfil y Q0 = 10 m³/s, según `test_b_dx5_single_matches_10m_case_physics`.
- **Diferencias declaradas:**
  - la fuente es una celda de 5 × 5 m (antes 10 × 10 m), así que la velocidad de entrada es 4× mayor;
  - la esponja U/V tiene 10 celdas = 50 m (antes 100 m);
  - el ancho mojado es 1190 m (antes 1180 m).
- **Antes de la corrida larga** se hace una corrida cronometrada de 10 min de modelo. Si la extrapolación supera 8 h, no se corre y se reporta.
- **Criterios de interpretación (C4 con el mismo código y la misma ventana):**
  - |Δ(diff/H)| < 0.03 entre 5 y 10 m → la resolución no domina C4, lo que favorece la clase (b)/(d).
  - diff/H baja en más de 0.05 → la resolución es importante (clase c).
  - En medio → no concluyente.
  - Los criterios OB 1–3 deben pasar; si no, D2 no se interpreta.

### Resultado de D1 (2026-11-06). Script: `analysis/d1_volume_flux.py`

**D1 tal como se pre-registró:** la mediana de Q_m/Q_t en la mitad inferior es **11.6 (twolayer) y 8.8 (sermilik)**. El modelo lleva mucho **más** flujo ascendente que la teoría, no menos. La hipótesis de sub-arrastre **no se sostiene con esta métrica**.

**Lectura de la métrica (post-hoc).** La franja i < 10 abarca los 1180 m de cara, así que mezcla la pluma con un ascenso que ocupa toda la cara. Separando ambas cosas, en los snapshots de 2 y 2.5 h:
- **Lejos de la pluma, bajo la interfaz, hay ascenso en toda la cara:** w medio de +1–2 cm/s en twolayer a z = −150 a −250 m (sermilik: +1 cm/s a −200/−250 m). Encima de la interfaz es casi 0. Probablemente es la capa límite convectiva del derretimiento ambiente, la misma que da el ~90 % de fondo en el derretimiento. Se verifica con los casos sin descarga.
- **Cerca de la pluma** (|y − y_src| ≤ 200 m, i < 10), Q va de 435 a 1445 m³/s en twolayer, frente a Q_t de 500 a 884 m³/s:
  - bajo la interfaz es mayor que la teoría, y ahí el ascenso de toda la cara aporta buena parte del flujo;
  - encima de la interfaz (z ≥ −110 m) es comparable, de 0.5 a 1.5 × Q_t.

**Conclusión de D1.** El fallo de C4 **no se explica por sub-arrastre** de la pluma. Queda abierto. Candidatos sin verificar:
- el ascenso de toda la cara cambia el ambiente efectivo de la pluma;
- la dinámica de la capa límite de pared en la cara (calor y sal de pkg/icefront) mantiene agua de la cara ascendiendo hasta la superficie.

D2 (resolución) y los casos sin descarga (¿existe el ascenso de toda la cara sin pluma?) son los siguientes discriminantes. C4 sigue en FALLA.

### Resultado: exceso de la pluma con fondo limpio y ascenso sin pluma (2026-11-06)

Script: `analysis/excess_clean_bg.py` → `figures/fig2_excess_clean.json`. Ventana 2–3 h, `face_total_melt`. x0c = 57.14 m (nominal en sermilik).

**Ruido y resolubilidad (pre-registrado):**

| | M_fondo (m³/s) | ε | E_single | 10·ε | ¿Resoluble? | M_fondo/M_single |
|---|---|---|---|---|---|---|
| twolayer | 1.2801 | 1.8e-4 | 0.1036 | 1.8e-3 | sí | 0.925 |
| sermilik | 1.2445 | 7.2e-4 | 0.1182 | 7.2e-3 | sí | 0.913 |

**R_E(x0) = E(x0)/(2·E_single), informativa, sin umbral:**

| x0 (m) | 0 | 20 | 40 | 80 | 160 | 320 |
|---|---|---|---|---|---|---|
| x0/x0c | 0 | 0.35 | 0.70 | 1.40 | 2.80 | 5.60 |
| twolayer | 0.762 | 0.953 | 0.918 | 0.927 | 0.952 | 0.983 |
| sermilik | 0.783 | 0.910 | 0.890 | 0.908 | 0.955 | 0.994 |

**Lectura (cualitativa, como se pre-registró):**
- Para x0 ≥ 20 m (x0/x0c ≥ 0.35), R_E ≈ 0.89–0.99: en el modelo, dos plumas separadas derriten casi como dos plumas independientes. No aparece la caída hacia ≈ 0.5 que muestra el laboratorio para x0 < x0c.
- Solo la fusión total (x0 = 0, una fuente de 2Q0) baja a 0.76–0.78.
- Con el fondo limpio desaparece la caída en x0 = 320 m del estimador anterior (0.85 → 0.98): confirma que era un artefacto del estimador de mediana.
- Sigue en pie el supuesto declarado de fondo aditivo.
- Cuidado al comparar con el laboratorio: con dx = 10 m, x0 = 20 m son solo dos celdas entre fuentes. Si la interacción depende de la resolución, esto no se puede saber sin D2 o equivalente.

**Ascenso en toda la cara sin pluma (verificación de la hipótesis de D1, post-hoc):** mediana a lo largo de la cara de w medio (2–3 h) en i = 0, en las corridas sin descarga:
- twolayer: ≈ 0 sobre −90 m; +0.8 cm/s a −120 m; +2.2 a +3.0 cm/s entre −150 y −240 m; +1.9 cm/s a −270 m.
- sermilik: ≈ 0 sobre −150 m; +1.1 a +2.3 cm/s entre −180 y −270 m.

Confirmado: el ascenso a lo largo de toda la cara existe sin pluma. Es la capa límite convectiva del derretimiento ambiente bajo la interfaz (agua más cálida). Explica el Q_m/Q_t ≈ 9–12 de D1, que medía sobre toda la cara. Esto no explica todavía el fallo de C4, que sigue abierto hasta D2.

### Resultado de D2: convergencia de C4 con dx = 5 m (2026-11-06)

Corrida: `/Volumes/LaCie/portfolio_runs/exp1/b_dx5_twolayer_single`. Llegó a 3 h (10 800 pasos) en 18 240 s de reloj con 5 rangos; la estimación previa era 23 328 s.

**Criterios OB 1–3 (requisito para interpretar): PASAN.**
- 1. `require_ok`: PASA.
- 2. Razón OB/interior máx 0.48, con columnas 196–198 frente a 20–179 (los mismos 100–900 m que i = 10–89 a 10 m). Con 10–189 da 0.41. A 10 m daba 0.52.
- 3. |η̄| máx 1.2e-16 m; error del transporte de la OB 1.9e-12 m³/s.

**C4, mismo código (`c4_for`) y misma ventana de 2–3 h:**

| | z_detach modelo | z_detach teoría | diff/H |
|---|---|---|---|
| 10 m (`b_3d_twolayer_single`) | −1.9 m | −66.1 m | 0.214 |
| 5 m (`b_dx5_twolayer_single`) | −55.7 m | −66.4 m | 0.036 |

**Veredicto pre-registrado: Δ(diff/H) = −0.179, una caída de más de 0.05, así que la resolución es importante (clase c).** El fallo de C4 a 10 m se atribuye a la resolución, no a física ausente ni a la teoría.

**Perfil de w (máximo en y, i = 0, media 2–3 h), para entender el cambio:**
- A 10 m, la pluma sigue a 0.35 m/s en z = −15 m y llega a la superficie.
- A 5 m, w cae de 0.57 m/s (−97.5 m) a 0.11 (−62.5) y a 0.02 (−52.5), y es ~0 encima. Se desprende bajo la interfaz, como predice la teoría.

**Lo que esto NO dice:**
- C4 sigue siendo **FALLA** a la resolución de producción (10 m); no se reetiqueta. El valor de 5 m está bajo el umbral, pero C4 se definió para la campaña.
- Solo dos resoluciones y una estratificación (twolayer). No demuestra convergencia: 5 m puede seguir sin converger, aunque la diferencia con la teoría ya es del 3.6 % de H.
- Hay diferencias declaradas: la fuente de 5 × 5 m tiene una velocidad de entrada 4× mayor, y la esponja es de 50 m.

**Consecuencia para (b):** todas las corridas de la campaña a 10 m tienen la pluma llegando a la superficie por falta de resolución. R_E y la Figura 2 a 10 m deben llevar esta advertencia o rehacerse a 5 m. Una corrida de 5 m cuesta ~5 h; la campaña completa (14 corridas, ~70 h) excede el presupuesto de 24 h. Decisión de Esteban.

### Subconjunto de (b) a 5 m: fijado antes de correr (2026-11-06, opción 1 aprobada por Esteban)

**Autorización.** Esteban aprobó ~20 h de cómputo adicional, por encima del presupuesto de 24 h de Exp. 1.

**Casos** (config `b_twoplume3d_dx5`, twolayer, misma física que a 10 m según `test_b_dx5_subset_matches_10m_cases`):
- `b_dx5_bg_twolayer_nodischarge` (Q0 = 0);
- `b_dx5_twolayer_x0_0`, `b_dx5_twolayer_x0_20`, `b_dx5_twolayer_x0_80`;
- `b_dx5_twolayer_single` ya existe (D2).

Orden: fondo, x0 = 20, x0 = 0, x0 = 80. Cada uno tarda ~5 h con 5 rangos.

**Requisitos por corrida:** criterios OB 1–3 (razón en columnas 196–198 frente a 20–179). Si una corrida falla, no se interpreta.

**Métrica (la misma que a 10 m, informativa, sin umbral):**
- R_E(x0) = E(x0)/(2·E_single), con E = M − M_fondo a 5 m, ventana 2–3 h.
- Se mantiene la regla ε: resoluble solo si E_single ≥ 10·ε.
- Se reporta junto a R_E a 10 m y se compara cualitativamente con el ≈ 0.5 del laboratorio.
- También se reporta C4 de cada corrida de dos plumas (z_detach del máximo en y de w), como contexto.

**Lectura fijada ahora:**
- Si R_E a 5 m difiere de 10 m en más de 0.1 para algún x0, la interacción entre plumas depende de la resolución, y la Figura 2 se presenta a 5 m.
- Si no, la conclusión "las plumas casi no interactúan" es robusta a la resolución dentro de este rango.

### Resultado del subconjunto de (b) a 5 m (2026-11-08)

Corridas en `/Volumes/LaCie/portfolio_runs/exp1/b_dx5_*`:
- fondo: 40 767 s de reloj;
- x0 = 20: 22 838 s, aceptada por R41 (falló `MPI_Finalize` tras la última escritura; ver `pd2_convergence.md`);
- x0 = 0: 19 129 s;
- x0 = 80: 20 979 s.

Script: `analysis/excess_dx5.py` → `figures/fig2_excess_dx5.json`.

**Criterios OB 1–3 (requisito): PASAN las cuatro corridas** (columnas 196–198 frente a 20–179):

| corrida | 1 `require_ok` | 2 razón OB/interior máx | 3 \|η̄\| máx; error de transporte |
|---|---|---|---|
| fondo (Q = 0) | PASA | 0.10 | 1.5e-18 m; 1.3e-13 m³/s |
| x0 = 20 | PASA | 0.52 | 7.4e-16 m; 2.0e-12 m³/s |
| x0 = 0 | PASA | 0.51 | 4.5e-16 m; 2.0e-12 m³/s |
| x0 = 80 | PASA | 0.88 | 6.4e-16 m; 3.0e-12 m³/s |

**Ruido (regla ε):** M_fondo = 1.2893 m³/s; ε = 2.4e-4; E_single = 0.0639 ≥ 10·ε = 2.4e-3, así que es **resoluble**. M_fondo/M_single = 0.953. A 10 m, E_single valía 0.1036: el exceso de la pluma única baja 38 % con la resolución.

**R_E(x0), ventana 2–3 h:**

| x0 (m) | 0 | 20 | 80 |
|---|---|---|---|
| R_E a 5 m | **0.977** | **1.204** | **0.947** |
| R_E a 10 m | 0.762 | 0.953 | 0.927 |
| diferencia | +0.215 | +0.251 | +0.020 |
| \|dif\| > 0.1 | sí | sí | no |

**Lectura fijada de antemano: se cumple la condición de dependencia con la resolución** (x0 = 0 y x0 = 20). La interacción entre plumas depende de la resolución, y la Figura 2 se presenta a 5 m. La conclusión de 10 m («las plumas casi no interactúan; solo la fusión total baja a 0.76») **no es robusta a la resolución**.

**C4 de cada corrida de dos plumas (solo contexto, sin umbral).** z_despegue del máximo en y de w en la cara. La teoría usa el ambiente lejano medio de la ventana:

| x0 | modelo | teoría Q0 | teoría 2Q0 |
|---|---|---|---|
| 0 | −9.3 m | −66.3 m | −37.8 m |
| 20 | −17.1 m | −66.5 m | −37.8 m |
| 80 | −55.1 m | −66.4 m | −37.7 m |
| `_single` (D2) | −55.7 m | −66.4 m | — |

Con x0 ≤ 20 m (≤ 4 celdas a 5 m), las plumas juntas suben más que una sola y casi llegan a la superficie. Superan incluso la teoría para 2Q0 en una fuente puntual. Con x0 = 80 m, cada pluma despega como la única.

**Sensibilidad a la ventana (post-hoc).** R_E en las mitades 2–2.5 h y 2.5–3 h:
- x0 = 0: 0.972 / 0.982;
- x0 = 20: 1.203 / 1.206;
- x0 = 80: 0.943 / 0.951.

Es estable dentro de 0.01.

**Interpretación (cautelosa):**
- A 5 m, R_E ≈ 1 o mayor para x0 ≤ 20 m: dos plumas cercanas derriten igual o más que dos independientes. Esto es lo contrario de la caída hacia ≈ 0.5 del laboratorio para x0 < x0c.
- Un mecanismo posible, sin verificar: la pluma fusionada o casi fusionada tiene más flujo de flotación, sube por encima de la interfaz y pone en contacto con la cara un flujo más rápido a lo largo de más profundidad. El x0 = 0 a 10 m no lo mostraba porque ahí todas las plumas llegaban a la superficie (C4 FALLA a 10 m).
- R_E > 1 en x0 = 20 y no en x0 = 0 queda **sin explicar**. Hay que verlo en el campo de derretimiento por profundidad antes de afirmar nada.
- Limitaciones:
  - solo tres x0 y una estratificación (twolayer);
  - dos resoluciones: no demuestra convergencia a 5 m;
  - se mantiene el supuesto de fondo aditivo;
  - la fuente de 5 × 5 m tiene una velocidad de entrada 4× mayor que a 10 m (diferencia declarada en D2).
