# PD1 — Puente con la tesis: `bay_thesis`, V1 y R1k (Task 3)

Fecha: 2026-11-11. Rama `exp4-arctic-bayiw`. Corridas en `/Volumes/LaCie/portfolio_runs/exp4/`.

## 1. Qué corrida de la tesis es la referencia

- La corrida guardada con bahía y estratificación "realistic" es `BayIW_Rectan/run_expand/` (21 de junio de 2026). Lo comprueban tres cosas:
  1. su `data` real (el que quedó en el directorio tras la etapa 2) es idéntico a `data.stage2`, y `data.stage1` es idéntico a la copia de `thesis_copy/`;
  2. el T inicial de `state.0000000000.t001.nc` coincide con `feb_temp_50zlev_560x352.bin` con diferencia 0.0 (la estratificación "realistic", febrero);
  3. el `dxF` de `grid.t001.nc` es el de `bahia_01_expand_dx.bin` actual (malla estirada, 5 a 13.27 km), no el de los `*.bak_uniforme`.
- Las carpetas `Output_IWBay/{realistic,linear,2layer}` están vacías; los productos derivados (`npz_out/`) no son campos completos.
- **Salida guardada:** MNC en float64, 8 teselas de 140 × 176, `state` (Eta, T, S, U, V, W) cada 900 s en 0–120 h (481 registros por variable y tesela), `diag_rho` (RHOAnoma) cada 900 s y los pickups MDS `pickup.0000007200` (fin de la etapa 1) y `pickup.ckptA` (iteración 14400, fin de la etapa 2).
- Por tanto **V1-exacto se puede evaluar tal como está preregistrado**, y además los pickups se pueden comparar byte a byte.

## 2. Configuración `bay_thesis` (R1)

R1 es la física de la tesis sin cambios. Solo cambia la E/S (ruling: la comparación de V1 es de valores):

- `code/`: el de la tesis (`thesis_copy/BayIW_Rectan/code/`), con dos cambios que no tocan la física:
  - `packages.conf` sin `mnc`;
  - `SIZE.h`: las mismas 8 teselas de 140 × 176, pero en 4 procesos (`nSx = 2`, `nPx = nPy = 2`) en vez de 8. La tesis corrió con 8 procesos; aquí el límite es 5. Con `GLOBAL_SUM_ORDER_TILES` (activo por omisión, `eesupp/inc/CPP_EEOPTIONS.h:127`), las sumas globales se hacen en el mismo orden de teselas con cualquier reparto entre procesos (`eesupp/src/global_sum_tile.F:161-191`). Por eso se espera el mismo resultado bit a bit. Con 5 procesos habría que cambiar las teselas, y con ellas el redondeo.
- Mismo compilador y opciones que la tesis: GNU Fortran 15.1.0 (Homebrew), `mpif77 -fconvert=big-endian -fallow-argument-mismatch -O3 -ftree-vectorize -funroll-loops` (el `Makefile` de la tesis y el nuestro coinciden en `FC`, `FFLAGS`, `FOPTIM` y `DEFINES`). Mismo árbol de MITgcm (checkpoint69j). El optfile difiere (`darwin_arm64_gfortran` frente a `darwin_arm64_gfortran_esteban`) solo en las rutas de NetCDF.
- Namelists (`config/bay_thesis/input_templates/`): `data` = `data.stage1` y `data.free` = `data.stage2`, con estos cambios de E/S, que la prueba `test_r1_namelists_differ_from_thesis_only_in_io` comprueba clave por clave:
  - `startTime`/`endTime` → `nIter0`/`nTimeSteps` (las mismas iteraciones: 0–7200 y 7200–14400);
  - `dumpFreq` 900 → 0: los volcados MNC de estado se sustituyen por diagnósticos MDS en los mismos instantes (`state3d`: UVEL, VVEL, THETA; `state2d`: ETAN; cada 900 s). `diag_rho` (RHOAnoma, 900 s) se conserva;
  - los nombres de los archivos de entrada (`bathy.bin`, `delX.bin`, …);
  - `data.pkg` sin `useMNC`.
- Entradas: batimetría, `delX`, `delY` y T/S "realistic" se copian byte a byte de la ruta de solo lectura de la tesis (prueba `test_r1_bathymetry_byte_identical_to_thesis`). El viento se regenera con `make_case.thesis_wind_patch` (la función de `wind_patch.py` con los argumentos del cuaderno `make_wind_forcing_local.ipynb`: `cx = 280`, `cy = 310`, `radius = 60`) y es idéntico byte a byte a `merid_local_bay_560x352.bin` (`test_r1_wind_*`).
- Dos etapas, como la tesis: la etapa 1 (60 h) con viento, `periodicExternalForcing` y pickup permanente en 216 000 s; la etapa 2 (60 h) sin viento ni condiciones iniciales y con `chkptFreq = 86400`. Reinicia de todos los `pickup*.0000007200.*`, incluidos los de Orlanski (E, S, W). `pfcommon.campaign.run_stages` admite ahora un plan de etapas explícito y un namelist por etapa.
- **Detalle de la tesis que se conserva:** el pulso de viento se muestrea en `t = linspace(0, 60, 180)` h (paso de 0.3352 h), pero MITgcm pone el registro i (desde 0) en (i + ½)·1200 s y los interpola linealmente (`eesupp/src/get_periodic_interval.F:118`). El registro máximo (i = 7, nominal 2.35 h, 0.0988 N m⁻²) cae así en t = 2.5 h del modelo, y el pulso del modelo dura un 0.56 % menos que el nominal. Es (d), error de la referencia, pequeño; se reproduce tal cual.

### Profundidad de la columna (PD0 §9, punto 4)

La suma del `delR` de la tesis es 599.9 m, pero la plataforma de su batimetría mide 1000 m. MITgcm trunca la columna en el fondo de la malla:
- `model/src/ini_masks_etc.F:112-113` acota la fracción de cada celda a [0, 1], y solo hay celdas k = 1…Nr;
- `ini_masks_etc.F:126-141` recalcula R_low = rF(1) − Σ drF·hFacC, que no puede bajar de 599.9 m;
- `model/src/ini_depths.F:119` lee R_low de `bathyFile`.

La profundidad efectiva fuera de la bahía es, entonces, de 600 m: el `Depth` máximo de `grid.t001.nc` de la tesis es 599.9 m. **No se cambia en R1** ni en R1k, R2 y R3, que usan la misma malla vertical. `case.json` lo registra (`column_depth_m`, `bathy_max_depth_m`, `effective_shelf_depth_m`).

## 3. R1b: control del build de KPP

R1k necesita un build distinto, `bay_thesis_kpp`:
- `kpp` en `packages.conf`;
- `OLx = OLy = 4`, porque con `KPP_SMOOTH_SHSQ`/`KPP_SMOOTH_DBLOC` (`pkg/kpp/KPP_OPTIONS.h:18,22`) KPP_CALC necesita uVel, vVel y el forzamiento en (−2:sNx+4) (`pkg/kpp/kpp_calc.F:100-104`; la parada de `kpp_check.F:161-167` está compilada fuera, así que el error no saltaría solo);
- `numDiags = 7*Nr` en `DIAGNOSTICS_SIZE.h` (tesis: 20*Nr). Con 20*Nr y 2 teselas por proceso, los datos estáticos del ejecutable con KPP pasan de 2 GB (2.12 GB de `__DATA`). Así chocan con la región compartida de dyld y el ejecutable falla al arrancar, de forma aleatoria (6 de 10 intentos, "syscall to map cache into shared region failed"). Con 7*Nr bajan a 1.83 GB y arranca 20 de 20 veces. Los diagnósticos de R1k usan 6*Nr + 2 niveles.

R1b es R1 (mismos namelists) corrido con ese build. **En 2 h (240 pasos), R1b es idéntico bit a bit a R1** en `state2d`, `state3d` y `diag_rho` (`pfcommon.campaign.compare_bitwise`, 8 registros cada uno). Compilar KPP, el halo de 4 y el tamaño de los diagnósticos no cambian la solución. Lo que distingue R1k de R1 son, por tanto, solo los namelists de la tabla siguiente.

## 4. R1k frente a R1

R1k = R1 + lo que KPP necesita (PD0 §6). La malla vertical de la tesis no cambia. La difusividad de fondo ya es la de L94 en la tesis (`diffKrT = diffKrS = 1e-5`), así que no cambia. La prueba `test_r1k_activates_kpp_and_differs_from_r1_only_in_listed_keys` compara los namelists clave por clave contra esta tabla:

<!-- R1k-R1 -->
| archivo | grupo | clave | R1 | R1k | motivo |
|---|---|---|---|---|---|
| `data.pkg` | PACKAGES | useKPP | — | `.TRUE.` | activa KPP |
| `data` | PARM01 | viscAr | `1.E-5` | `1.E-4` | viscosidad de fondo de KPP (PD0 §6, L94 ec. 29; `pkg/kpp/kpp_routines.F:1197`) |
| `data` | PARM01 | ivdc_kappa | `1.` | — | KPP no admite otra convección (`pkg/kpp/kpp_check.F:127-133`); la convección la hace KPP (difmcon) |
| `data` | PARM01 | implicitViscosity | — | `.TRUE.` | KPP la exige (`pkg/kpp/kpp_check.F:141-145`) |
| `data.free` | PARM01 | viscAr | `1.E-5` | `1.E-4` | ídem, etapa libre |
| `data.free` | PARM01 | ivdc_kappa | `1.` | — | ídem |
| `data.free` | PARM01 | implicitViscosity | — | `.TRUE.` | ídem |
| `data.free` | PARM03 | nTimeSteps | `7200` | `21600` | R1k corre 240 h (60 h forzadas + 180 h libres); H4d usa 0–120 h |
| `data.kpp` | KPP_PARM01 | * | — | present | grupo vacío: todos los parámetros de KPP por omisión (PD0 §6) |
| `data.diagnostics` | DIAGNOSTICS_LIST | fields(1,4) | — | `'KPPhbl'` | salida: capa límite de KPP, 900 s |
| `data.diagnostics` | DIAGNOSTICS_LIST | fileName(4) | — | `'kpp2d'` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | frequency(4) | — | `-900.` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | timePhase(4) | — | `0.` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | fields(1:2,5) | — | `'KPPviscA','KPPdiffT'` | salida: viscosidad y difusividad de KPP, cada 3 h |
| `data.diagnostics` | DIAGNOSTICS_LIST | fileName(5) | — | `'kpp3d'` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | frequency(5) | — | `-10800.` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | timePhase(5) | — | `0.` | salida |
<!-- /R1k-R1 -->

Además, el código: el build `bay_thesis_kpp` (§3), que R1b muestra neutro.

## 5. R2 y R3 (puente; generados y probados, se corren tras el PD1)

Ruling R9: se generan y se prueban aquí, y se corren en la campaña cuando el PD1 fije T_fin. En `cases.yaml` están con 280 h, de forma provisional. Son R1 con estos cambios (edit sets `kpp`, `pd0_mixing`, `f78n` y `eos_vmf`):
- **KPP**, como en R1k;
- **mezcla del PD0 §6:** sin `viscAh`, con `viscC2smag = 2.2`, `viscC4smag = 1`, `diffKhT = diffKhS = 0`, `bottomDragQuadratic = 2.5e-3` y `tempAdvScheme = saltAdvScheme = 77`;
- **f de 78°N en plano f:** `f0 = 1.4266e-4` y `beta = 0`. **La tesis no era un plano f** (§6), así que R1 → R2 también quita β;
- **EOS lineal ajustada** a Van Mijenfjorden (`configs.fjord_vmf.eos`): tα = 3.0456e-5, sβ = 8.4068e-4, `rhoNil = rhoConst = 1026.96` (`rhoConst` explícito; por omisión sería `rhoNil`, `model/src/ini_parms.F:483`), y `tRef = 50*0.4676`, `sRef = 50*33.6395`;
- **perfil de verano de VMF** interpolado a los 50 niveles de la tesis, repitiendo el último valor bajo 240 m, y estabilizado con la EOS ajustada. Coincide con `figures/profiles_vmf.json` (`summer_on_thesis_levels`, Task 2) con error < 1e-12 (`test_r2_f0_78n_and_summer_profile_on_thesis_levels`);
- **R3:** la batimetría sin bahía de la tesis (`nobahia_01_expand_bat.bin`; sus `dx`/`dy` son idénticos a los de la bahía).

R1 → R2 cambia a la vez f (y β), perfil, EOS y mezcla; está declarado en el PD0.

## 6. V1: R1 frente a la corrida guardada de la tesis

`analysis/v1_thesis.py --run /Volumes/LaCie/portfolio_runs/exp4/R1` → `figures/v1_thesis.json`. R1 pasó `require_ok` en las dos etapas (14 400 pasos; 3428 s + 3377 s de reloj con 4 procesos).

- **V1-exacto (congelado): PASA.** máx |Δη| = **0.0 m** en los **480 registros comunes** (iteraciones 0–14 370, cada 900 s, dominio completo). Tolerancia: 1×10⁻¹⁰ m. Los 480 registros son idénticos bit a bit.
- U, V y T en 41 de esos registros (uno de cada 12, más el último), en todo el dominio y los 50 niveles: máx |ΔU| = máx |ΔV| = 0.0 m s⁻¹ y máx |ΔT| = 0.0 °C.
- **Pickups, byte a byte:** `pickup.0000007200.data` (fin de la etapa 1) y `pickup.ckptA.data` (iteración 14 400, fin de la etapa 2) son idénticos a los de la tesis. El estado completo (U, V, T, S, η y las tendencias de Adams–Bashforth) es el mismo al final de cada etapa.
- **V1-espectral** (solo descriptivo, porque V1-exacto pasa): en la estación de la boca (Seiche2, j = 325, i = 280), el pico de la seiche en η está en el mismo bin (0.267 cph, Δf = 1/120 cph) y la DEP difiere un 0.0 %.
- **Lectura:** el mismo código, las mismas opciones del compilador y el mismo orden de las sumas globales dan el mismo resultado, aunque cambie el reparto de las teselas entre procesos (8 → 4). La E/S (MNC → MDS) no toca la solución. No hay nada que clasificar en (a)–(d).

## 7. La discrepancia de f₀ (Task 3, Step 6)

**Resuelta: la tesis corrió en un plano β, no en un plano f.**
- El `data` con que corrió la salida guardada tiene `f0 = 6.97E-5` y no fija `beta`. El eco del namelist en `run_expand/STDOUT.0000` (líneas 790–798) dice `f0 = 6.97E-05`, `beta = 1.0E-11` y `selectCoriMap = 1` ("Beta-Plane", línea 928).
- `beta = 1e-11 m⁻¹ s⁻¹` es el valor por omisión de MITgcm (`model/src/set_defaults.F:112`), y en malla cartesiana f = f0 + β·yC (`model/src/ini_cori.F:77`).
- Con `ygOrigin = −1621 km` y la malla estirada, y va de −1614 km a +215 km. f va de 5.36×10⁻⁵ s⁻¹ en la frontera sur (32.6 h) a 7.19×10⁻⁵ en la norte (24.3 h). El `fCori` de `grid.t001.nc` de la tesis vale 5.3556×10⁻⁵ en j = 0, igual al calculado.
- **En la bahía** (j = 324–347) f = 7.05–7.16×10⁻⁵ s⁻¹ (24.4–24.8 h). En la boca (j = 323), f = 7.04×10⁻⁵ (24.78 h), un 1 % más que f0. Los 25.04 h de f0 corresponden a y = 0, en j ≈ 308, 15 celdas al sur de la costa.
- **Los "~31 h" de los scripts de la tesis** (`CORIOLIS_F = 0.0322` cph en `visualizar_espectros_metodologia.py`, `figura_seiche_series_psd.py` y `visualizar_eta_frontera_sur.py`) corresponden al f local en j ≈ 39 (y ≈ −1350 km), cerca de la frontera sur, no a la bahía. Los scripts posteriores de la tesis ya lo marcan como error (`picos_seiche_inercial.py:50`, `figura_mapas_estaciones_psd.py:100`), pero lo corrigen a f0 suponiendo un plano f ("plano f, `data` &PARM01"). Es (d), error de la referencia.
- **Espectro de R1** (*post hoc*, descriptivo): el pico de la energía cinética superficial en 0.02–0.07 cph está en 0.0417 cph en las 7 estaciones cercanas a la bahía (j = 265–345). Queda a < 0.4 bins del f local y a 0.2 bins de f0, y a 1.1 bins de 0.0322 cph. Con 120 h (Δf = 0.0083 cph) el espectro no distingue f0 del f local de la bahía (difieren 0.06 bins). Lo que decide es el `data` y el `fCori` guardados. En las estaciones del sur (j = 60–200), lejos del parche de viento, no hay un pico inercial claro (picos en 0.025 y 0.058 cph).
- **Consecuencias:**
  - R1 y R1k conservan el plano β, que es la física de la tesis.
  - R2/R3 usan plano f a 78°N (`beta = 0`, plan, "Latitud y f"). El cambio R1 → R2 incluye, entonces, quitar β; se declara.
  - Los números de la teoría de la bahía de la tesis (`theory.json`, T_i = 25.0 h) usan f0; el f local de la bahía es un 1–3 % mayor. No cambia ninguna banda ni ningún umbral del PD0, pero se declara para el PD1.

## 8. R1k y costo medido

- **R1k corrió completa:** 240 h, 28 800 pasos en dos etapas (60 h forzadas y 180 h libres; build `bay_thesis_kpp`, 4 procesos), estado OK, `require_ok` en las dos etapas. Reloj: 4097 s + 12 242 s = 4.5 h. Tras cada corrida, `pgrep` no encontró ningún proceso de MPI ni del modelo.
- **Descriptivo**, no es H4d: η en la boca (Seiche2), 0–120 h:
  - desviación típica de 4.49×10⁻⁴ m en R1 frente a 4.50×10⁻⁴ m en R1k;
  - diferencia máxima entre las dos: 4.8×10⁻⁵ m.
- KPPhbl: la capa límite llega a 16 m durante el pulso (máximo en el dominio), baja a < 2 m en ~30 h y queda en 0.5–0.6 m en la fase libre.
- H4d (η = F_bc,out/W_viento, R1k frente a R1, 0–120 h) se evalúa en la Task 7 con `flux.section_flux` (Task 5). Las salidas que necesita están en las dos corridas: U, V, T, RHOAnoma y η cada 900 s, en todo el dominio.
- **Costo medido** (4 procesos, núcleos de rendimiento):
  - R1: 0.47 s por paso; 1.9 h para 120 h;
  - R1k: 0.57 s por paso; 4.5 h para 240 h;
  - corridas cortas de 2 h (R1, R1b, R1k): 124, 128 y 188 s. Incluyen el arranque y la E/S de 440 MB de entradas.

  Total de la Task 3: ~6.6 h de las 30 h del presupuesto.
  Para R2/R3 a 280 h, con KPP y la mezcla del PD0: ~5.3 h cada una, por extrapolación de R1k, que no lleva la mezcla del PD0. Están por debajo de las 8 h por corrida.
- **Disco:** R1 ocupa 150 GB y R1k 308 GB en LaCie (snapshots 3D cada 900 s, como la tesis). R2/R3 a 280 h ocuparán ~360 GB cada una.

## 7. R2 y R3 corridos (Task 6, parte B; 2026-11-12/13)

- **Smoke de 2 h (ruling R12):** `R2_smoke2h` y `R3_smoke2h` (240 pasos, build `bay_thesis_kpp`, 4 rangos) salen OK. `DIAGNOSTICS_SET_POINTERS` usa 302 niveles, con `numDiags` = 350 (7*Nr), así que alcanza. `STDERR.0000` vacío, máx |w| ≈ 1×10⁻⁵ m s⁻¹, CFL W_hf ≈ 9×10⁻⁴. Reloj: 213 y 219 s. **No hizo falta recompilar.**
- **280 h** (`run_campaign.py --config bay_thesis_kpp --cases 'R[23]' --np 4 --max-concurrent 1 --skip-ok`), una tras otra, las dos OK en las dos etapas (33 600 pasos):

  | Corrida | Etapa forzada (s) | Etapa libre (s) | Total | s/paso | Disco |
  |---|---|---|---|---|---|
  | R2 | 5892 | 20 396 | **7.30 h** | 0.78 | 357 GB |
  | R3 | 5416 | 19 913 | **7.04 h** | 0.75 | 357 GB |

- Las dos quedan por debajo del límite de 8 h por corrida, pero por encima de las ~5.3 h que se habían extrapolado de R1k (0.57 s/paso).
- Parte de la diferencia se debe a que la evaluación de V2–V4 leía LaCie a la vez durante R2. R3, casi sin lecturas, salió a 0.75 s/paso.
- Tras cada etapa, `pgrep` no encontró procesos de MPI ni del modelo. Quedan 972 GB libres en LaCie.
