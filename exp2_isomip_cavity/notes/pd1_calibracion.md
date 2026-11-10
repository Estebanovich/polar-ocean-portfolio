# PD1a — Cronometraje, paso de tiempo y calibración de ΓT (Exp. 2, Task 4 fase A)

Fecha: 2026-11-08. Criterios vinculantes: `pd0_preregistro.md` (K1, V1, N1–N3). Presupuesto B1: ≤ 32 h de cómputo para todo el Exp. 2 y ≤ 8 h por etapa de corrida. Corridas en `/Volumes/LaCie/portfolio_runs/exp2/`, logs en `…/exp2/logs/`.

## 1. Corridas de cronometraje

Todas con hilos = 1, `pgrep -fl mpirun` vacío antes y después, y ΓT = 0.03 (el valor no afecta el costo). El costo por paso es el "Wall clock time" de la sección `MAIN LOOP` del propio modelo, dividido por el número de pasos.

| Corrida | Malla | np × simultáneas | dt (s) | Duración | s/paso |
|---|---|---|---|---|---|
| `cal2_gT0.03_smoke30d_walls` (Task 2) | 2 km | 5 × 1 | 200 | 30 d | 0.01955 |
| `cal2_gT0.03_tim30d_dt400` | 2 km | 5 × 1 | 400 | 30 d | 0.02040 |
| `tim4m_gT0.03_tim10d` | 4 km | 5 × 1 | 400 | 10 d | 0.00580 |
| `tim4m_gT0.03_tim30d_dt400` | 4 km | 5 × 1 | 400 | 30 d | 0.00660 |
| `tim4m_gT0.03_tim30d_dt800` | 4 km | 5 × 1 | 800 | 30 d | 0.00619 |
| `ens_T*_zc400_s1.0_tim10d` (4 de 5, ver §5) | 4 km | 1 × 4 | 400 | 10 d | 0.02125 c/u |
| `ens_T*_s0.8_tim10d` | 4 km | 1 × 5 | 400 | 10 d | 0.02385 c/u |

- La caja np = 5 a 4 km usa `config/ocean0_4km_mpi` (caso nuevo `tim4m_gT0.03`, misma física que `ocean0_4km`).
- El costo por paso de los primeros 10 días subestima el de 30 días (0.00580 contra 0.00660 a np = 5: factor 1.14, más iteraciones de cg2d). Para extrapolar se aplica ese factor a las corridas serie de 10 días.
- **Rendimiento a 4 km** (dt = 800, por corrida de 2 años, con el factor 1.14): np = 5 en serie, 0.136 h por corrida; np = 1 × 5 simultáneas, 0.594 h por tanda de 5, es decir **0.119 h por corrida**. El arreglo de mayor rendimiento es np = 1 con 5 simultáneas (13 % mejor).

## 2. Paso de tiempo por CFL

**Regla** (escrita después de ver los CFL a dt = 200/400 s y antes de correr los dt dobles): dt = el mayor de la sucesión de duplicaciones desde el valor inicial (200 s a 2 km, 400 s a 4 km) cuyo CFL máximo en 30 días sea ≤ 0.5 en advección (u, v, w con hFac) y en ondas internas. El barotrópico es implícito. Monitor cada hora en las corridas de dt doble y diario en las de Task 2.

| Malla, dt | CFL u | CFL v | CFL w (= W_hf) | CFL ondas internas | f·dt |
|---|---|---|---|---|---|
| 2 km, 200 s (30 d) | 0.084 | 0.057 | 0.134 (día 22) | 0.065 | 0.028 |
| 2 km, **400 s** (30 d) | 0.170 | 0.114 | **0.271** (día 22) | 0.131 | 0.056 |
| 4 km, 400 s (30 d) | 0.094 | 0.053 | 0.159 (día 26) | 0.065 | 0.056 |
| 4 km, **800 s** (30 d) | 0.186 | 0.106 | **0.317** (día 30) | 0.131 | 0.112 |

- **CFL de ondas internas:** c₁·dt/Δx, con c₁ = velocidad del primer modo baroclínico, calculada del perfil discreto de 36 niveles con la EOS lineal (problema de autovalores N²w = c²·(−d²w/dz²)). Ocean0: c₁ = 0.654 m s⁻¹ (H = 720 m, igual a NH/π). Ensamble (rampa, d = 300 m): el máximo es c₁ = 0.91 m s⁻¹ (Tbot = −0.5, z_c = 400), que da un CFL de 0.091 (4 km, 400 s) y 0.18 (4 km, 800 s).
- **Duplicar otra vez** (800 s a 2 km, 1600 s a 4 km) daría un CFL w ≈ 0.54 y 0.63, por encima de 0.5. Ese dt queda descartado sin correrlo, porque el CFL escala linealmente con dt.
- **Elección: dt = 400 s a 2 km y 800 s a 4 km.** dt/Δx es el mismo en las dos mallas. Sustituye los valores provisionales de la Task 2 (200 y 400 s) en `cases.yaml`.
- **Sensibilidad al dt** (30 días, ΓT = 0.03, `pd1_smoke_check.py`):

  | | dt | melt plataforma (m a⁻¹) | melt profundo (m a⁻¹) | S1 / S2 / S4 / S5 | S3 |
  |---|---|---|---|---|---|
  | 2 km | 200 | 19.827 | 57.781 | PASS (S2 2.31e-6) | FAIL |
  | 2 km | 400 | 19.836 (+0.05 %) | 57.811 (+0.05 %) | PASS (S2 2.30e-6) | FAIL (T mín −2.507) |
  | 4 km | 400 | 20.498 | 61.327 | PASS (S2 1.77e-6) | FAIL |
  | 4 km | 800 | 20.488 (−0.05 %) | 61.265 (−0.10 %) | PASS (S2 1.75e-6) | FAIL (T mín −2.449) |

  S3 falla igual que en la Task 2 (R10: el criterio no se redefine). Las demás condiciones del smoke se mantienen con el dt nuevo.
- **Limitación:** el CFL w a 4 km y 800 s alcanza su máximo el día 30, el último de la corrida, así que puede crecer más adelante. Se vigila en las corridas de 2 años (§4).

## 3. Costo del Exp. 2 contra B1

### 3.1 Estimación con el s/paso de 30 días

Estimación previa a la cola (2026-11-08). **Supuesto:** el costo a 2 km usa el s/paso de la corrida de 30 días (0.02040 s/paso a 400 s, es decir 0.894 h por corrida de 2 años). Para 4 km ver §1.

| Bloque | Corridas | Arreglo | Horas |
|---|---|---|---|
| Calibración 2 km | 4 × 2 años (3 + verificación) | np = 5, una a la vez | 3.58 |
| Calibración 4 km | 3 × 2 años simultáneas + verificación | np = 1 (×3, luego ×1) | ≤ 1.06 |
| Ocean1 2 km | 20 años = 4 etapas × 5 años (2.24 h c/u ≤ 8 h) + prueba de reinicio | np = 5 | 8.94 + ≤ 0.25 |
| Ensamble 4 km | 28 × 2 años, tandas 5,5,5,5,5,3 | np = 1 × 5 | 3.56 |
| Ya gastado | smokes de la Task 2 y cronometraje | — | ≈ 0.3 |
| **Total** | | | **≈ 17.7 h ≤ 32 h (B1 se cumple)** |

- Con una segunda iteración de K1 en las dos mallas (+0.89 + 0.53 h), el total sería ≈ 19.1 h.
- **Para comparar:** con los dt provisionales (200 y 400 s), el mismo plan costaría ≈ 33 h (6.85 + 2.1 + 17.1 + 7.0), por encima de B1.
- La cota de "calibración 4 km" usa el costo con 4 simultáneas (0.02125 s/paso × 1.14), que es una cota superior para 3 simultáneas y para una corrida serie sola.

**Arreglo de la calibración a 4 km:** np = 1, 3 simultáneas, y la verificación `ocean0_4km` también serie. np = 1 es el arreglo de mayor rendimiento medido (§1). Con solo 3 corridas, np = 5 una tras otra terminaría ≈ 0.1 h antes (3 × 0.136 h = 0.41 h contra ≈ 0.53 h), y la verificación sola ≈ 0.3 h antes. Se prefiere np = 1 porque así la calibración, la verificación y el ensamble usan el mismo binario. La verificación `ocean0_4km` es además la referencia del ensamble. La diferencia (≤ 0.5 h) no afecta a B1.

### 3.2 Corrección con las corridas de 2 años (ruling de la revisión de código)

Tiempos de pared medidos de las primeras corridas de calibración (`MAIN LOOP`, todas OK):

| Corrida | Arreglo | Pared (s) | s/paso |
|---|---|---|---|
| `cal2_gT0.015` | np = 5 | 3181.2 | 0.02017 |
| `cal2_gT0.03` | np = 5 | 3180.2 | 0.02017 |
| `cal2_gT0.06` | np = 5 | 3198.9 | 0.02029 |
| `cal4_gT0.015`, `0.03`, `0.06` | np = 1 × 3 simultáneas | 1530.1 / 1532.1 / 1531.5 | 0.01942 |

- **2 km:** 0.885 h por corrida de 2 años, 1 % menos que lo estimado con 30 días (0.894 h).
- **4 km:** la tanda de 3 tarda 0.425 h. La cota de §3.1 (0.53 h) era conservadora.

**Total corregido**, ya con la enmienda de §4.2 (2 corridas más de 2 años en cada malla):

| Bloque | Horas |
|---|---|
| Gastado: smokes de la Task 2 y cronometraje | 0.30 |
| Gastado: tanda `cal4` (3 corridas) | 0.43 |
| Gastado: `cal2` (3 × 0.885) | 2.66 |
| `cal4_gT0.0075` + verificación `ocean0_4km` (serie, solas, ≤ 0.425 c/u) | ≤ 0.85 |
| `cal2_gT0.0075` + verificación `ocean0_2km` (2 × 0.885) | 1.77 |
| Ocean1, 2 km (10 × 0.885; 4 etapas de 2.21 h) + prueba de reinicio | 8.85 + ≤ 0.25 |
| Ensamble 4 km (estimación de §3.1, conservadora) | 3.56 |
| **Total** | **≈ 18.7 h ≤ 32 h (B1 se cumple)** |

Con la iteración adicional de K1 en las dos mallas (+0.885 + 0.425 h), el total sería ≈ 20.0 h, también dentro de B1.

## 4. Calibración K1

### 4.1 Primera cola (3 puntos): NOT_BRACKETED → PD2

La cola `scripts/queue_calibration.sh` corrió del 2026-11-08 20:40 al 23:46 CST, en el commit 3c3a568. Logs en `…/exp2/logs/queue1/`. Las 6 corridas terminaron OK, y `calibrate.py` las leyó con `require_ok`: N3 pasa en todas.

| Malla | ΓT | ⟨m⟩ profundo, meses 19–24 (m a⁻¹) | ⟨m⟩ plataforma (m a⁻¹) | N1 |
|---|---|---|---|---|
| 4 km | 0.015 | 41.216 | 13.595 | 0.0014 |
| 4 km | 0.03 | 60.187 | 18.327 | 0.0021 |
| 4 km | 0.06 | 78.899 | 22.806 | 0.0026 |
| 2 km | 0.015 | 39.910 | 13.208 | 0.0056 |
| 2 km | 0.03 | 57.574 | 17.465 | 0.0097 |
| 2 km | 0.06 | 74.680 | 21.761 | 0.0025 |

- **Resultado mecánico de K1:** el objetivo de 30 m a⁻¹ **no queda acotado** a ninguna resolución, porque el melt mínimo muestreado es 41.2 (4 km) y 39.9 (2 km). pd0 prohíbe extrapolar, así que `calibrate.py` devolvió NOT_BRACKETED (código 3) y la cola no lanzó verificaciones. Esto fue a PD2.
- **Registro:** `figures/calibration_{4km,2km}_3pt.json`.
- **Extrapolación log–log informativa** con los dos primeros puntos. No es la calibración ni sustituye a ninguna corrida:
  - 4 km: pendiente 0.546, ΓT ≈ 0.0084;
  - 2 km: pendiente 0.529, ΓT ≈ 0.0087.

  Las dos caen por debajo del límite inferior de V1 (0.011).

### 4.2 Enmienda del muestreo (aprobada por Esteban, 2026-11-09, PD2)

- **Muestreo:** se añade ΓT = 0.0075, con ΓS = ΓT/35, a las dos resoluciones (casos `cal2_gT0.0075` y `cal4_gT0.0075` en `cases.yaml`). El conjunto queda en {0.0075, 0.015, 0.03, 0.06}.
- **Criterios:** K1 (30 ± 2 m a⁻¹, melt profundo de los meses 19–24) y V1 (0.011 ≤ ΓT_cal,2km ≤ 0.2) **no se redefinen**. `pd0_preregistro.md` no se edita.
- **Procedimiento:** por lo demás es el de pd0.
  - Interpolación de ln⟨m⟩ contra ln ΓT, lineal por tramos, en el tramo que acota 30 m a⁻¹, sin extrapolar.
  - Una corrida de verificación.
  - Como máximo, una iteración más, con todos los puntos (los 4 muestreos y la verificación).
  - Si sigue fallando, se vuelve a PD2.
  - Si 0.0075 tampoco acota los 30 m a⁻¹, el resultado vuelve a ser NOT_BRACKETED y PD2.
- **Si ΓT_cal,2km < 0.011:** V1 se reporta como **FAIL** y se clasifica como **(c), elección numérica**.
  - Motivo: el modelo es de nivel z con celdas superiores de 20 m (capa límite de Losch de una celda).
  - Gwyther et al. (2020, *Ocean Modelling* 147:101569; `literatura/Gwyther_2020 - Vertical Processes and Resolution Impact Ice Shelf Basal Melting.md`, L30–31) concluyen: "Z-coordinate models, which generally have coarser vertical resolution in ice shelf cavities, may simulate higher melt rates compared to terrain-following coordinate models". En L34 añaden que los modelos coinciden "when the effective vertical resolution adjacent to the ice shelf base is similar".
  - Un melt mayor a igual ΓT exige un ΓT_cal menor.
  - V1 no se redefine. Esta clasificación solo explica el fallo.
- **Cola:** `scripts/queue_calibration.sh`, ampliada y relanzada (ver `task-4-report`). Corre las dos corridas de 0.0075, `calibrate.py` con los 4 puntos de cada malla y las verificaciones `ocean0_4km` / `ocean0_2km` con `--set gammaT=<ΓT_cal>`. Luego evalúa K1 de la verificación (`calibrate.py --verify`, 28 ≤ ⟨m⟩ ≤ 32). La iteración adicional (`calibration_<res>_iter2.json`, corrida `<caso>_iter2`) solo se lanza si K1 falla, como permite pd0.
- **Resultados:** §4.3.

### 4.3 Resultado y evaluación (segunda cola, 2026-11-09 00:41–03:13 CST, commit f7fdd41)

La evaluación es mecánica: `analysis/pd1_calibration_eval.py` escribe `figures/pd1_calibration_eval.json`, y los veredictos se copian en `figures/calibration_{4km,2km}.json` (clave `evaluation`) y en `figures/verify_{4km,2km}.json`. Todas las corridas se leyeron con `metrics.load_run`.

**Calibración (4 puntos por malla).** En las dos mallas, el tramo que acota los 30 m a⁻¹ es [0.0075, 0.015].

| Malla | ΓT_cal (bruto) | ΓT_cal (usado, 5 cifras) | Pendiente local log–log |
|---|---|---|---|
| 4 km | 0.0094686201 | **0.0094686** | 0.690 |
| 2 km | 0.0099627405 | **0.0099627** | 0.698 |

La pendiente local es ln(m₀.₀₁₅/m₀.₀₀₇₅)/ln 2.

**Corridas** (melt de los meses 19–24, en m a⁻¹; CFL máximo de los 2 años con el monitor diario, 731 líneas):

| Corrida | ΓT | Melt profundo | Melt plataforma | N1 | N2 | N3 | CFL u / v / w (R17) |
|---|---|---|---|---|---|---|---|
| cal4_gT0.0075 | 0.0075 | 25.541 | 9.288 | 0.0144 ✓ | 4.1e-8 ✓ | ✓ | 0.095 / 0.070 / 0.189 |
| cal4_gT0.015 | 0.015 | 41.216 | 13.595 | 0.0014 ✓ | 4.4e-8 ✓ | ✓ | 0.131 / 0.083 / 0.262 |
| cal4_gT0.03 | 0.03 | 60.187 | 18.327 | 0.0021 ✓ | 9.8e-8 ✓ | ✓ | 0.176 / 0.105 / 0.334 |
| cal4_gT0.06 | 0.06 | 78.899 | 22.806 | 0.0026 ✓ | 1.5e-7 ✓ | ✓ | 0.211 / 0.135 / 0.431 |
| **ocean0_4km** (verificación) | 0.0094686 | **30.410** | 10.661 | 0.0129 ✓ | 4.8e-8 ✓ | ✓ | 0.104 / 0.075 / 0.215 |
| cal2_gT0.0075 | 0.0075 | 24.610 | 9.104 | 0.0100 ✓ | 1.1e-6 ✓ | ✓ | 0.116 / 0.091 / 0.192 |
| cal2_gT0.015 | 0.015 | 39.910 | 13.208 | 0.0056 ✓ | 1.4e-6 ✓ | ✓ | 0.145 / 0.096 / 0.288 |
| cal2_gT0.03 | 0.03 | 57.574 | 17.465 | 0.0097 ✓ | 4.0e-6 ✓ | ✓ | 0.168 / 0.128 / 0.362 |
| cal2_gT0.06 | 0.06 | 74.680 | 21.761 | 0.0025 ✓ | 4.6e-6 ✓ | ✓ | 0.186 / 0.158 / 0.414 |
| **ocean0_2km** (verificación) | 0.0099627 | **30.408** | 10.711 | 0.0098 ✓ | 1.3e-6 ✓ | ✓ | 0.128 / 0.095 / 0.228 |

En todas las corridas, W_hf coincide con w. **Ningún CFL supera 0.5** (R17); el mayor es 0.431, en w de `cal4_gT0.06`.

**Veredictos:**
- **K1, 4 km: PASA.** Melt profundo 30.41 m a⁻¹, dentro de [28, 32].
- **K1, 2 km: PASA.** Melt profundo 30.41 m a⁻¹, dentro de [28, 32]. No hizo falta la iteración adicional (pd0 solo la permite si K1 falla).
- **V1: FALLA.** ΓT_cal,2km = 0.0099627 < 0.011. Clase **(c), elección numérica**, según la enmienda (§4.2): modelo de nivel z con celdas superiores de 20 m (Gwyther 2020 L30–31). V1 no se redefine.
- **N1–N3: PASAN en las 10 corridas.** N1 máximo 0.0144 (< 0.03); N2 máximo 4.6e-6 (< 0.01).
- ΓT_cal,2km = 0.0099627 queda fijo para Ocean0 y Ocean1 a 2 km, y ΓT_cal,4km = 0.0094686 para el ensamble. No se vuelven a tocar (pd0 K1).

**Sobreenfriamiento (R15)** en las corridas de verificación. ΔT = T − T_f(S, p), con el liquidus de AD16 Eq. 25 y p = ρ₀g|z_c| en el centro de la celda. Se dan el estado final (instantáneo) y, entre paréntesis, la media de los meses 19–24.

| | ocean0_4km | ocean0_2km |
|---|---|---|
| Celdas con ΔT < 0 / celdas mojadas | 225 / 41 249 (209) | 1015 / 164 922 (1041) |
| Mínimo de ΔT (°C) | −0.217 (−0.171) | −0.339 (−0.305) |
| Celdas con ΔT < −0.1 | 56 (49) | 585 (596) |
| Filas de pared (y = 0 y y = 80 km) | 15, mín −0.012 | 0 (3), mín −0.0004 |
| Interior | 210, mín −0.217 | 1015, mín −0.339 |
| Vecino lateral seco, fuera de las filas de pared | 78, mín −0.217 | 189, mín −0.304 |
| Sin vecino seco | 132, mín −0.168 | 826, mín −0.339 |
| Celda superior delgada (hFacC < 0.5) | 123, mín −0.217 | 483, mín −0.339 |
| Otras celdas | 102, mín −0.161 | 532, mín −0.299 |
| Celda superior de la columna, cualquier hFac | **225 de 225** | **1000 de 1015** |
| Celda más fría | k = 8 (z = −170 m), hFacC = 0.20, T = −2.184, S = 33.54 | k = 10 (z = −210 m), hFacC = 0.24, T = −2.335, S = 33.50 |

- **Dónde está:** casi todo el sobreenfriamiento está en la celda mojada superior, justo bajo el hielo, y es más intenso en las celdas delgadas. No se concentra en las paredes de y = 0 / 80 km, a diferencia del congelamiento en las paredes laterales de Y26 (L1187). A 2 km hay más celdas y el mínimo es más frío que a 4 km.
- **Clasificación:** limitación de tipo (b)/(c), como prevé R15: capa límite de Losch de una celda, sin frazil. Es menor que en el smoke de 30 días (1592 celdas, −0.49 °C, con ΓT = 0.03), lo que concuerda con un ΓT menor.

**Prueba de resolución (pd0, descriptiva).** Ocean0 a 2 km contra 4 km, cada uno con su ΓT_cal, meses 19–24:

| Métrica | 2 km | 4 km | (4 − 2)/2 |
|---|---|---|---|
| Melt medio de la plataforma (m a⁻¹) | 10.711 | 10.661 | −0.5 % |
| Máx. overturning (Sv) | 0.1936 | 0.1952 | +0.9 % |
| Máx. \|ψ_bt\| en la cavidad (Sv) | 0.2654 | 0.6219 | **+134 %** |
| Correlación espacial del melt (mapa de 2 km promediado 2×2 a 4 km; 707 celdas con plataforma en alguno de los dos) | — | — | r = 0.968 |

- La diferencia de melt de la plataforma (0.5 %) está por debajo del 25 %, así que no hay bandera de pd0. Esta coincidencia es, en gran parte, el objetivo de la calibración (ΓT calibrado por separado a cada resolución al mismo melt profundo, K1), no evidencia independiente de convergencia.
- **Observación descriptiva** (sin umbral; no se interpreta aquí): el máximo de |ψ_bt| en la cavidad (y > 50 km, x < 640 km) es 2.3 veces mayor a 4 km. El melt y la overturning, en cambio, casi coinciden. Hechos (revisión de código, comprobados con `metrics.barotropic_streamfunction`):
  - el máximo a 4 km está en y = 52 km, x = 608 km, dentro de un giro amplio 2–3 veces más intenso que a 2 km; a 2 km, el máximo está en y = 64 km, x = 592 km;
  - el revisor señala un zigzag de escala de malla cerca de y ≈ 52–60 km. En el perfil de ψ_bt a x = 608 km se ve además un escalón no monótono en y ≈ 68–72 km (0.142 → 0.191 Sv).
  - **Antes de usar ψ_bt del ensamble** (4 km), hay que revisar el ruido de escala de malla en U a 4 km. No se atribuye causa aquí.

## 5. Hallazgo aparte (no afecta a la calibración)

`ens_T0.0_zc400_s1.0_tim10d` no arrancó: `INI_THETA: found 10142 wet grid-pts with theta=0 identically` (STOP en `ini_theta.F`, `checkIniTemp`).
- Con Tbot = 0.0, el perfil rampa da T = 0 exactamente bajo z_c + d/2, y MITgcm lo toma como un campo sin inicializar.
- Afecta a los 3 casos del ensamble con Tbot = 0.0 (z_c = 250, 400, 550).
- No afecta a Ocean0 ni a Ocean1.
- Resuelto por el ruling R18 (§6).
- Esta corrida fallida no entra en el cronometraje: la tanda de 5 se repitió con los casos s = 0.8.

## 6. Ruling R18: `checkIniTemp=.FALSE.`

Desde el 2026-11-09, antes de cualquier corrida nueva, la plantilla `config/ocean0/input_templates/data` (PARM05, la única con ese grupo) lleva `checkIniTemp=.FALSE.`. Fuente: `/Volumes/ESTEBAN_MIT/MITgcm`, solo lectura.

**Por qué no cambia la solución:**
- En `model/src/ini_theta.F:110-119`, `checkIniTemp` solo decide qué pasa cuando hay puntos mojados con θ = 0 exactamente. Con `.TRUE.` detiene la corrida (`STOP 'ABNORMAL END: S/R INI_THETA'`); con `.FALSE.` escribe un aviso (`:120-125`). θ no se modifica en ningún caso.
- El otro uso está en `ini_theta.F:130`: `IF ( checkIniTemp .AND. allowFreezing )`, que recorta θ a −1.9 °C. Ese bloque solo actúa con `allowFreezing = .TRUE.`. Su valor por omisión es `.FALSE.` (`model/src/set_defaults.F:220`), y ninguna plantilla lo cambia; lo comprueba la prueba `test_r18_check_ini_temp_off_without_freezing`.
- Por eso, en esta configuración, el valor de `checkIniTemp` no altera el estado inicial ni la integración.

**Corridas anteriores:** las 6 primeras de calibración (`cal4_gT{0.015,0.03,0.06}`, `cal2_gT{0.015,0.03,0.06}`), los smokes y el cronometraje corrieron con el valor por omisión (`.TRUE.`). Ninguna tenía θ = 0 (salvo `ens_T0.0_zc400_s1.0_tim10d`, que no arrancó, §5), así que el cambio no tiene efecto físico sobre ellas. Se reutilizan con `--skip-ok`, que compara deltaT, endTime y gammaT (anulado por R21: todas esas corridas usaron la ley errónea de u\* y ya no se reutilizan).

## Error (a): forma de u\* (R21)

Ruling R21, aprobado por Esteban el 2026-11-09. Clase **(a), error nuestro**. Detalle en `pd1_config.md`, sección «Error (a): forma de u\* (R21)».

**Qué estaba mal.** AD16 Eq. 27 (GMD 9, p. 2485 = p. 15 del PDF, verificada en la página impresa) es u\*² = C_D,top·(u_w² + u_tidal²), con piso u\* = √C_D·u_tidal = 5e-4 m s⁻¹. El modelo calculaba u\*² = Cd·u_w² + u_tidal², con piso 0.01 m s⁻¹. **Origen:** la fórmula del plan se copió de una extracción de texto del PDF que perdió los paréntesis.

**Evidencia.** u\* medio de la región profunda (base del hielo < −300 m) en Ocean0: MITgcm-JPL (`data_external/osf_jpl/Ocean0_yosi_COM_MITgcm_JPL_new.nc`, `frictionVelocity`, meses 7–12 de su año) **0.0040 m s⁻¹**; el nuestro, **0.0109 m s⁻¹** (corrida anulada `ocean0_2km`, ΓT = 0.0099627, media ponderada por área de `SHIuStar` mensual en los meses 19–24; a 4 km, `ocean0_4km`, 0.0109 también). Un u\* inflado ~2.7 veces exige un ΓT menor para dar 30 m a⁻¹; es una explicación candidata del ΓT_cal bajo de §4.3 (0.0095–0.0100) y del fallo de V1, que solo la nueva calibración puede confirmar o descartar.

**Corridas ANULADAS** (se conservan en disco, en `/Volumes/LaCie/portfolio_runs/exp2/`, solo como registro histórico; ningún resultado suyo se reutiliza):

- los smoke: `cal2_gT0.03_smoke30d`, `cal2_gT0.03_smoke30d_run1_mixedsnap`, `cal2_gT0.03_smoke30d_walls` (y `cal2_gT0.03_tim30d_dt400`);
- toda la calibración: `cal2_gT{0.0075,0.015,0.03,0.06}` y `cal4_gT{0.0075,0.015,0.03,0.06}` (todos los ΓT);
- las verificaciones `ocean0_2km` y `ocean0_4km`;
- el cronometraje `tim4m_*` (los costos medidos siguen valiendo como costo por paso: el cambio no altera el número de operaciones);
- la serie diagnóstica R20, `diag4_*` (`diag4_adv33`, `diag4_hfac025`, `diag4_ustartopdz`, `diag4_convadj`).

Los cronometrajes cortos del ensamble (`ens_*_tim10d`) también usaron la ley errónea; solo servían para medir costo y tampoco se reutilizan.

**Condicionado a la ley errónea:**
- la **enmienda de PD2** (§4.2, ΓT = 0.0075 en las dos mallas) y todo §4.1–§4.3 (ΓT_cal,4km = 0.0094686, ΓT_cal,2km = 0.0099627, K1 PASA, V1 FALLA con clase (c)). Se conservan como registro y **no se reutilizan**: la enmienda no se aplica de forma automática a la nueva calibración;
- la **serie diagnóstica R20** (`pd1_diagnostico_capa_limite.md`), anulada completa.

**Lo que no cambia:** K1, V1, N1–N3 y todos los criterios de pd0 siguen tal como están pre-registrados (`pd0_preregistro.md` no se edita). El paso de tiempo de §2 (dt = 400 s a 2 km, 800 s a 4 km) y los costos por paso de §1 y §3 se mantienen: el cambio no altera el costo de un paso; el CFL de las nuevas corridas se vuelve a revisar (R17).

### Relanzamiento de la calibración (desde cero, procedimiento de pd0 K1)

- **Muestra:** solo ΓT ∈ {0.015, 0.03, 0.06}, ΓS = ΓT/35, a 2 km y a 4 km. Luego `calibrate.py` (interpolación log–log, sin extrapolación), una corrida de verificación por malla y, solo si K1 falla, la única iteración adicional con los 4 puntos.
- **Si 30 m a⁻¹ no queda acotado** en una malla, la cola se detiene para esa malla (sin verificación, sin extrapolar) y se va a PD2. Las corridas de muestra de la otra malla siguen, porque son las de pd0 de todos modos.
- **Manejo de las corridas anuladas:** las nuevas corridas llevan el sufijo de directorio **`_r21`** (`cal4_gT0.015_r21`, …, `ocean0_4km_r21`, `ocean0_2km_r21`; iteración: `<caso>_r21_iter2`). Las anuladas se quedan donde están, sin tocar: `--skip-ok` solo mira el directorio con sufijo, así que nunca las recoge, y `calibrate.py` recibe explícitamente los directorios `_r21`. No se movieron a `void_r21/`: el sufijo basta para aislarlas y evita copiar y borrar ~14 GB en el exFAT (copia y verificación sin beneficio). Los resultados nuevos van a `figures/calibration_<res>_r21.json` y `figures/verify_<res>_r21.json` (los JSON anulados se conservan); los logs, a `logs/r21/`. `pd1_calibration_eval.py` evalúa por omisión los directorios `_r21` y la muestra de pd0.
- **Cola:** `scripts/queue_calibration.sh`, lanzada con nohup; PID, log y resultados se anotan en la sección siguiente cuando termine.
- **Candado contra corridas anuladas (después de la revisión).** `cases.yaml` lleva `melt_law: AD16_eq27_r21` en la configuración base, que `make_case` escribe en `case.json`. `run_campaign.py --skip-ok` compara además `melt_law`: una corrida OK sin esa clave (ley anulada) no se salta ni se sobrescribe; la campaña se niega a arrancar. `metrics.load_run` rechaza cualquier corrida cuyo `melt_law` no sea el vigente, así que ningún análisis (calibración, ensamble, Ocean1) puede leer una corrida anulada. El sufijo por omisión de `run_campaign.py` es ahora `_r21`, y `calibrate.py` escribe por omisión `figures/*_<res>_r21.json`. A las 8 corridas `_r21` (hechas en el commit 4c5a655 justo antes del candado) se les añadió `melt_law` en `case.json` y `input/case.json` después de comprobar que su commit es 4c5a655 y que su ejecutable es byte a byte el binario R21 (sha256 `f69062d2…` a 2 km, `4ddad51a…` a 4 km); la clave `melt_law_backfill` lo registra. Las dos verificaciones (`ocean0_4km_r21`, `ocean0_2km_r21`) tienen `git_dirty: True` en `case.json`: cuando se generaron sus entradas, el árbol tenía sin rastrear los JSON de calibración que la propia cola acababa de escribir (`figures/calibration_*_r21.json`); el código era el del commit 4c5a655, y su ejecutable es el mismo binario R21.

## Resultados R21

Cola `scripts/queue_calibration.sh`, PID 13965, 2026-11-09 08:11–12:28 CST, commit 4c5a655. Log en `…/exp2/logs/r21/queue_calibration_r21.log`. Evaluación mecánica contra pd0: `analysis/pd1_calibration_eval.py` (por omisión, corridas `_r21`) → `figures/pd1_calibration_eval_r21.json`; calibración en `figures/calibration_{4km,2km}_r21.json` y verificación en `figures/verify_{4km,2km}_r21.json`. Todas las corridas se leyeron con `metrics.load_run` (`require_ok` y el candado de `melt_law`).

**Calibración (3 puntos de pd0 por malla).** 30 m a⁻¹ queda acotado en las dos mallas, sin enmienda.

| Malla | Tramo que acota | ΓT_cal (bruto) | ΓT_cal (usado, 5 cifras) |
|---|---|---|---|
| 4 km | [0.015, 0.03] | 0.029355867 | **0.029356** |
| 2 km | [0.03, 0.06] | 0.030959375 | **0.030959** |

**Corridas** (melt de los meses 19–24, m a⁻¹; CFL máximo de los 2 años con el monitor diario, 731 líneas):

| Corrida | ΓT | Melt profundo | Melt plataforma | N1 | N2 | N3 | CFL u / v / w (R17) |
|---|---|---|---|---|---|---|---|
| cal4_gT0.015_r21 | 0.015 | 14.160 | 5.424 | 0.0076 ✓ | 6.3e-8 ✓ | ✓ | 0.072 / 0.059 / 0.127 |
| cal4_gT0.03_r21 | 0.03 | 30.737 | 10.409 | 0.0045 ✓ | 4.1e-8 ✓ | ✓ | 0.107 / 0.077 / 0.212 |
| cal4_gT0.06_r21 | 0.06 | 53.922 | 16.580 | 0.0032 ✓ | 8.5e-8 ✓ | ✓ | 0.161 / 0.100 / 0.304 |
| **ocean0_4km_r21** (verificación) | 0.029356 | **30.023** | 10.186 | 0.0002 ✓ | 4.3e-8 ✓ | ✓ | 0.106 / 0.075 / 0.208 |
| cal2_gT0.015_r21 | 0.015 | 12.988 | 5.033 | 0.0112 ✓ | 3.1e-7 ✓ | ✓ | 0.098 / 0.081 / 0.140 |
| cal2_gT0.03_r21 | 0.03 | 29.239 | 9.960 | 0.0094 ✓ | 8.5e-7 ✓ | ✓ | 0.135 / 0.098 / 0.216 |
| cal2_gT0.06_r21 | 0.06 | 51.502 | 15.846 | 0.0022 ✓ | 3.9e-6 ✓ | ✓ | 0.209 / 0.125 / 0.348 |
| **ocean0_2km_r21** (verificación) | 0.030959 | **30.002** | 10.191 | 0.0083 ✓ | 1.2e-6 ✓ | ✓ | 0.137 / 0.098 / 0.220 |

W_hf coincide con w en todas. **Ningún CFL supera 0.5** (R17); el mayor es 0.348, en w de `cal2_gT0.06_r21`.

**Veredictos (pd0, sin cambios en los criterios):**
- **K1, 4 km: PASA.** Melt profundo 30.02 m a⁻¹, dentro de [28, 32].
- **K1, 2 km: PASA.** Melt profundo 30.00 m a⁻¹, dentro de [28, 32]. No hizo falta la iteración adicional.
- **V1: PASA.** ΓT_cal,2km = 0.030959, dentro de [0.011, 0.2]. Informativo: queda dentro del rango de los modelos z de Y26 (0.0135–0.0325) y cerca del de MITgcm-JPL (0.0325).
- **N1–N3: PASAN en las 8 corridas.** N1 máximo 0.0112 (< 0.03); N2 máximo 3.9e-6 (< 0.01); N3 sin fallos.
- ΓT_cal,2km = 0.030959 queda fijo para Ocean0 y Ocean1 a 2 km, y ΓT_cal,4km = 0.029356 para el ensamble (pd0 K1). Sustituyen a los valores anulados de §4.3.

**Sobreenfriamiento (R15, estado final de las verificaciones; no es criterio).**
- 4 km: 109 de 41 249 celdas mojadas bajo T_f(S, p) local; mínimo −0.149 °C; 7 por más de 0.1 °C. Ninguna en las filas de pared; las 109 son celdas superiores de columna (80 delgadas, hFacC < 0.5, entre ellas las 7 de < −0.1 °C); 48 tienen vecino lateral seco. La más fría está en k = 11 (z = −230 m), hFacC = 0.2.
- 2 km: 184 de 164 922 celdas; mínimo −0.115 °C; 5 por más de 0.1 °C. Ninguna en las filas de pared; todas son celdas superiores (127 delgadas); 48 con vecino lateral seco. La más fría, en k = 9 (z = −190 m), hFacC = 0.49.

**Prueba de resolución (pd0, descriptiva, sin umbral).** 4 km contra 2 km, cada una con su ΓT_cal, meses 19–24:
- melt medio de la plataforma: 10.186 contra 10.191 m a⁻¹ (−0.05 %); no hay bandera de pd0 (> 25 %). **Esta coincidencia es, en gran parte, el objetivo de la calibración:** ΓT se calibró por separado a cada resolución para dar el mismo melt profundo de Ocean0 (K1, 30 m a⁻¹). No es evidencia independiente de convergencia. Las métricas de resolución que no se calibraron son ψ_bt (+125 %) y la correlación de los mapas de melt (r = 0.915), abajo;
- máximo de la overturning: 0.1845 contra 0.1851 Sv (−0.3 %);
- máximo de |ψ_bt| en la cavidad: 0.584 contra 0.260 Sv (+125 %, 2.2 veces). Es la misma observación de §4.3 (allí 2.3 veces) y sigue pendiente la revisión del ruido de escala de malla en U a 4 km antes de usar ψ_bt del ensamble;
- correlación espacial de los mapas de melt (2 km promediado 2 × 2 sobre la malla de 4 km): r = 0.915, en 707 celdas.


## Producción (PD1a aprobada por Esteban, 2026-11-09)

- ΓT_cal escrito en `cases.yaml` (ya no depende de `--set`): 0.029356 en la configuración `ocean0_4km` (los 28 `ens_*` y `ocean0_4km`), 0.030959 en `ocean0_2km` (`ocean0_2km`, `ocean1_2km`). ΓS = ΓT/35 lo deriva `make_case`. Los casos `cal*` y `diag4_*` conservan sus valores propios. Prueba: `test_calibrated_gammas_in_registry`.
- Entradas generadas en seco para los 31 casos (28 `ens_*`, `ocean0_4km`, `ocean0_2km`, `ocean1_2km`): ΓT y ΓS correctos en `case.json` y en `data.shelfice`; `melt_law` presente; `checkIniTemp=.FALSE.` (R18). Ocean1: estado inicial COLD (AD16 Tabla 5, comprobado contra `profiles.cold`), restauración WARM (`T.restore.bin` = `profiles.ocean0`), 20 años = 1 576 800 pasos de 400 s, en 4 etapas de 394 200 pasos = 5 años; 5 años es múltiplo de `pChkptFreq` (1 año), de los diagnósticos mensuales (60) y de la serie de 5 días (365).
- Procedencia: desde ahora `run_campaign.py` escribe en `input/case.json` el sha256 del ejecutable con que corre cada caso (`binary_sha256`, `binary_path`, `binary_recorded_at`).
- Cola `scripts/queue_production.sh`: primero el ensamble (28 casos, 4 km, np = 1, 5 a la vez; `ocean0_4km_r21` se reutiliza como referencia con `--skip-ok`), luego N1–N3 sobre las 29 corridas (`analysis/ensemble_n_check.py`; un fallo de N3 queda como hueco y no se vuelve a correr), luego `ocean1_2km` (np = 5, 4 etapas). Logs en `…/exp2/logs/production/`.
- **Tbot = 0.0 (R18), comprobado:** los 3 casos `ens_T0.0_*` arrancan. Prueba de 1 día (108 pasos, `--end-time 86400`, sufijo `_r21_start1d`): OK, con el aviso `INI_THETA: found 10142 wet grid-pts with theta=0 identically` en lugar del STOP. Un primer intento de 10 pasos (sufijo `_r21_start10`) terminó normalmente en el modelo, pero quedó FAILED en `run_status`: el runner lee la última iteración del monitor diario, y 8000 s no llegan a una línea de monitor. Es un artefacto de la prueba, no del modelo; esos directorios quedan como registro. El `binary_sha256` registrado en la prueba coincide con el binario R21 de 4 km (`4ddad51a…`).
