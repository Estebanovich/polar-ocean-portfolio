# PD1 — Configuración (a) `a_line2d`: pluma lineal 2D no hidrostática

Task 4, Steps 1–3. Material para el Punto de decisión 1. **No se ha lanzado la campaña**: se espera el visto bueno de Esteban.

## 1. Tabla de parámetros

| Bloque | Parámetro | Valor | Fuente / razón |
|---|---|---|---|
| Malla | x | `nx = 193`; dx = 5 m en la cara, crece ×1.03 hasta 200 m; L = 20.14 km | `cases.yaml` (Task 3, `grids.stretched_x`); dx0 = 5 m está entre Xu 2012 (20 m, y prueba a 2 m) y Sciascia 2013 (10 m) |
| | última columna (i = 192) | tierra (`land_columns: [[0, -1]]`) | Cierra el periodo en x: la cara oeste de la columna 0 queda contra tierra, es decir, es la pared de hielo. Sin esto, `OBCS_CHECK` pide una OB oeste (`pkg/obcs/obcs_check.F:647-668`) |
| | y | `ny = 1`, dy = 10 m | Rebanada de 10 m como Sciascia 2013 (§2.1) |
| | z | `nz = 50`, dz = 10 m, H = 500 m, fondo plano | Sciascia 2013 (10 m); Xu 2012 usa 5 m |
| Tiempo | `deltaT` | **3 s** (provisional 2 s en Task 3) | Corrida corta: CFL máx. 0.20 (q0 = 0.03) y 0.32 (q0 = 0.1). Sciascia 2013 usa 1–20 s |
| | duración | 48 h → `nTimeSteps = 57600` | Brief |
| Dinámica | `nonHydrostatic` | `.TRUE.` | Xu 2012, Sciascia 2013 |
| | superficie libre | lineal, implícita, `exactConserv=.TRUE.` | |
| | `cg3dMaxIters` / `cg3dTargetResidual` | 500 / 1e-8 | Con 40 iteraciones (valor del tutorial) cg3d nunca convergía (residuo ~4e-5). Converge en 320–390 iteraciones |
| | `cg2dTargetResidual` | 1e-13 (≈170 iteraciones) | `tutorial_plume_on_slope` |
| Viscosidad | `viscAh` / `viscAr` | 0.25 / 1e-3 m² s⁻¹ | **Sciascia et al. 2013, Tabla 1** (2D NH, 10 m × 10 m en el frente) |
| Difusividad | `diffKhT = diffKhS` / `diffKrT = diffKrS` | 0.25 / 2e-5 m² s⁻¹ | Sciascia 2013, Tabla 1 (Pr horizontal = 1). Alternativa: Xu 2012 usa νv = 0.1 m² s⁻¹, biarmónica 300 m⁴ s⁻¹ y κ = ν |
| | `implicitViscosity/Diffusion` | `.TRUE.` | |
| Contornos | pared de hielo y fondo | no deslizamiento (`no_slip_sides`, `no_slip_bottom`) | Sciascia 2013 §2.1 ¶14 |
| Advección | `tempAdvScheme = saltAdvScheme` | 33 (DST3, con limitador) | `tutorial_plume_on_slope` |
| EOS | `eosType` | `'LINEAR'` | Idéntica a `theory/plume.py` (C2'') |
| | `tAlpha` / `sBeta` | 3.87e-5 K⁻¹ / 7.86e-4 psu⁻¹ | Jenkins 2011, Tabla 1 (Task 2) |
| | `rhoNil = rhoConst` | 1027 kg m⁻³ | Task 2 (`[verificar]`) |
| | `HeatCapacity_Cp` | 3974 J kg⁻¹ K⁻¹ | Task 2 / Sciascia 2013, Tabla 1 |
| | `gravity` | 9.81 | igual a la teoría |
| | `tRef`, `sRef` | 0 °C, 34 psu | Con EOS lineal solo desplazan ρ por una constante; g′ no cambia. sRef = 34 reduce la anomalía (redondeo) |
| icefront | archivos | `R_icefront.bin` (= H en la columna 0), `icefrontlength.bin` (= 1/dx₀ = 0.2 m⁻¹) | Task 3; `icefront_description.tex:12-16` |
| | ley de transferencia | **la de defecto**: γT = max(1e-2·√1.5e-3·\|wVEL\|, 1e-4) m s⁻¹, γS = 5.05e-3 γT | `icefront_thermodynamics.F:89-97` |
| | T del hielo | lineal: `ICEFRONTthetaSurface = −20 °C` en la superficie, 0 °C en la base del frente | `icefront_thermodynamics.F:116-123`; para C2'' la teoría debe usar `ice_depth = H` |
| | constantes (explícitas, iguales al defecto sin shelfice) | `rhoIcefront = 917`, `ICEFRONTlatentHeat = 3.34e5`, `ICEFRONTHeatCapacity_Cp = 2000`, `ICEFRONTkappa = 1.54e-6` (no se usa en la termodinámica) | `icefront_readparms.F:84-88` |
| | `applyIcefrontTendT/S` | `.TRUE.` | El agua de deshielo entra como flujo virtual de sal, sin volumen (`icefront_thermodynamics.F:151-154`) |
| addMass | `selectAddFluid` | 1 | |
| | `addMassFile` | `addmass.bin`: ρ0·q0·dy kg s⁻¹ en la celda (k = 49, i = 0), es decir, la celda mojada del fondo con `R_icefront` (fallo R8) | Task 3 |
| | `temp_addMass` / `salt_addMass` | 0 °C / 0 psu | brief |
| | requisitos (`model/src`) | `#define ALLOW_ADDFLUID` (`config_check.F:363-373`); sin `rigidLid` (`config_check.F:837-845`); `staggerTimeStep=.TRUE.` evita el aviso de desfase de 1 paso en T/S (`config_check.F:846-856`); `selectAddFluid` ∈ {−1, 0, 1, 2} (`config_check.F:828-835`). `convertFW2Salt` queda sin fijar y entonces vale −1 (`ini_parms.F:649-652`), así que `apply_forcing.F:505-507` usa `temp_addMass`/`salt_addMass` con superficie libre lineal. **No exige superficie libre no lineal ni `exactConserv`** (se usa `exactConserv` igual) | verificado en la fuente |
| OBCS | frontera este | `OB_Ieast = −2` → columna Nx−1 (`obcs_readparms.F:671-672`); Orlanski (`CMAX = 0.45`, `cvelTimeScale = 250 s`) + esponja U/V de 10 celdas (≈ 2 km), `Urelaxobcsbound = 60 s`, `Urelaxobcsinner = 1200 s`; sin esponja de T/S | Como la tesis (`BayIW_Rectan/run_expand/data.obcs`) |
| | `useOBCSbalance` | **apagado** | Solo equilibra EmPmR, no addMass (`obcs_balance_flow.F:300-315`), y atraparía la descarga. Comprobado: η medio sube 3.5 mm en 2 h (q0 = 0.03), contra 10.8 mm si el volumen quedara atrapado, y la subida se frena: el volumen sale por Orlanski |
| Salidas | diagnósticos MDS | instantáneas de `THETA`, `SALT`, `UVEL`, `WVEL`, `ETAN` cada 1 h; `ICFfwFlx`, `ICFhtFlx` promediados cada 1 h (`icefrontDiag`) | brief; nombres en `icefront_diagnostics_init.F:41,48` |
| | estado final | `T/S/U/V/W/Eta.<iter>.data` al inicio y al final (`dumpInitAndLast`) | Las instantáneas de diagnóstico se escriben un paso tarde (`diagnostics_write.F:62-72`): la de t = endTime **no existe** |
| | pickups | `pChkptFreq = 7200 s` | Para la prueba bit a bit C3 (Step 5) |
| E/S | | `readBinaryPrec = writeBinaryPrec = 64`, big-endian, `useSingleCpuIO`, sin MNC | restricciones globales |

## 2. Costo medido y extrapolación

Corrida corta oficial: `a_line2d_twolayer_q0.03`, 2 h de modelo, dt = 3 s (2400 pasos), np = 1, en `/Volumes/LaCie/portfolio_runs/exp1/a_line2d_twolayer_q0.03_short2h/`.

- `run_status.json`: `{"state": "OK", "last_iter": 2400, "n_timesteps": 2400, "wall_s": 65}`; sin NaN; `pgrep` vacío al terminar.
- Control de CFL del caso más fuerte, `a_line2d_twolayer_q0.1` (2 h, dt = 3 s): OK, 67 s.
- Compilación: ~65 s (serial).

| | por caso (48 h) | 6 casos |
|---|---|---|
| Extrapolación lineal (24 × 65–67 s) | ≈ 26–27 min | **≈ 2.6 h** (9 500 s) |

Supuesto: costo por paso constante. cg3d domina (320–390 iteraciones por paso) y podría crecer si el flujo se vuelve más turbulento. Margen: aun con +50 %, son ≈ 4 h, dentro de ≤ 8 h por corrida y ≤ 24 h para todo el Exp. 1 (que además incluye (b) y (c)).

Alternativa barata, a decidir en PD1: con `cg3dMaxIters = 40` (sin converger) la corrida corta tardó 34 s con dt = 2 s. La diferencia a 2 h frente a la versión convergida fue máx. 0.011 °C en T y 2 mm s⁻¹ en w. Se eligió la versión convergida.

## 3. Máximos de CFL (monitor, `advcfl_*`, 2 h)

| caso | dt | `advcfl_uvel_max` | `advcfl_wvel_max` | w máx. (m s⁻¹) | u máx. (m s⁻¹) |
|---|---|---|---|---|---|
| twolayer q0.03 | 3 s | 0.202 | 0.198 | 0.66 | 0.53 |
| twolayer q0.1 | 3 s | 0.247 | 0.318 | 1.06 | 0.69 |

Todos < 0.5. Con dt = 4 s, q0 = 0.1 llegaba a 0.42. Por eso se bajó a 3 s, igual para los 6 casos (no se ajusta por caso).

## 4. Figura rápida

`pd1_a_quicklook.png` (`pd1_quicklook.py`, `pfcommon.style`): T a las 2 h con contornos de w, en los primeros 1.5 km. La pluma sube pegada a la cara (w ≈ 0.2–0.6 m/s en las primeras celdas), se despega en la interfaz de dos capas (~200 m) y sale hacia el fjord en la interfaz, mezclando agua atlántica en la capa superior. El derretimiento medio en la cara durante la última hora es ICFfwFlx = −5.6e-3 kg m⁻² s⁻¹ (≈ 0.53 m d⁻¹ de hielo).

## 5. Puntos para que Esteban decida

1. dt = 3 s y cg3d convergido (≈ 2.6 h los 6 casos). ¿Se acepta, o se prefiere cg3d a 40 iteraciones (≈ 1 h)?
2. Viscosidad y difusividad de Sciascia 2013 (Laplaciana 0.25 / 1e-3). Xu 2012 (νv = 0.1, biarmónica) sería la alternativa de sensibilidad.
3. Pared sin deslizamiento (Sciascia) frente al arrastre cuadrático Cd u² de la teoría. Es una diferencia de modelo conocida para C2' (clase c).

## 6. Decisión PD1 y C3 (Step 5)

Esteban aprobó la configuración: cg3d convergido (500 iteraciones, 1e-8), dt = 3 s, viscosidad y difusividad de Sciascia 2013, pared sin deslizamiento (reportada como clase (c) frente al Cd·u² de la teoría).

**C3, reinicio bit a bit: cumplido.** `a_line2d_twolayer_q0.03`: 4 h seguidas (`c3_restart_4h`) frente a 2 h + pickup + 2 h (`c3_restart_2h_a`, `c3_restart_2h_b`). Los volcados finales `T` y `Eta` (iteración 4800) son idénticos con `np.array_equal`; también S, U y W.
