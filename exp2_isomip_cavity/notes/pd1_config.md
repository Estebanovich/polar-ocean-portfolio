# PD1 — Configuración MITgcm del Exp. 2 (ISOMIP+), ley de fusión y trazador

Fecha: 2026-11-08. Task 2. Fuente MITgcm: `/Volumes/ESTEBAN_MIT/MITgcm`, commit `8e1c40762` (solo lectura). Las líneas `pkg/…:N` se refieren a ese árbol. AD16 = Asay-Davis et al. 2016.

## 1. Lo que hace `pkg/shelfice` sin modificar (Step 1)

### 1.1 Coeficientes de transferencia

**No existe** la forma de AD16 (ΓT constante y adimensional, con γT = ΓT·u*). El paquete ofrece dos:

- **γ constante en m s⁻¹:** `SHELFICEheatTransCoeff` y `SHELFICEsaltTransCoeff` se copian en `shiTransCoeffT/S` (`pkg/shelfice/shelfice_init_varia.F:80-81`). Por omisión valen 1e-4 m s⁻¹ (`shelfice_readparms.F:118`), y γS = 5.05e-3·γT (`shelfice_readparms.F:117, 203-205`).
- **Holland y Jenkins (1999)**, con `SHELFICEuseGammaFrict` (`shelfice_thermodynamics.F:148-163, 470-493` del original): γT = u*/(Γturb + Γmole,T), con Γturb = ln(u*² ζN η*²/(|f|·5ν))/κ + 1/(2ζN η*) − 1/κ, y Γmole función de los números de Prandtl y Schmidt. Esta no es la forma de AD16.

### 1.2 u_tidal

**No existe.** El original calcula u*² = Cd·max(1e-6, |u|²) (`shelfice_thermodynamics.F:470`), es decir, un piso de |u| = 1e-3 m s⁻¹ en lugar de u_tidal = 0.01 m s⁻¹ (AD16 Eq. 27). Cd es `shiCdrag` (`shelfice_init_varia.F:95`, valor por omisión 0.0015, `shelfice_readparms.F:137`).

### 1.3 Signo y unidades de `SHIfwFlx` (lo consume la Task 3)

- **Unidades:** kg m⁻² s⁻¹. Título: "Ice shelf fresh water flux (positive upward)" (`pkg/shelfice/shelfice_init_fixed.F:249-255`).
- **Valor:** `SHIfwFlx = shelfIceFreshWaterFlux = freshWaterFlux = rhoConst·γS·(1 − S/Sb)` (`shelfice_thermodynamics.F:645-647` y `:660`; se llena en `:701`).
- **Signo: negativo = fusión.** Con fusión, Sb < S y el flujo es negativo. Lo dice el comentario de `shelfice_thermodynamics.F:651-653`: "downward (negative) fresh water flux implies melting".
- **Conversión al protocolo** (m a⁻¹ de agua dulce, positivo = fusión, `pd0_preregistro.md`): **m = −SHIfwFlx / ρfw · 31 536 000**, con ρfw = 1000 kg m⁻³ (AD16 Tabla 4). Está implementada en `theory/melt.py::melt_m_per_yr`.
- **Consistencia con AD16 Eq. 26:** SHIfwFlx = −ρfw·m, porque rhoConst·γS·(Sb − S)/Sb = −ρfw·m. Aquí rhoConst hace el papel de ρsw.
- **`SHIhtFlx`** [W m⁻²] es el flujo de calor hacia arriba, hacia el hielo (`shelfice_init_fixed.F:257-262`).
- **`SHIuStar`** [m s⁻¹] es u* (en esta configuración incluye u_tidal). **`SHIgammT/S`** [m s⁻¹] son γT y γS = Γ·u*.

### 1.4 Punto de fusión

- **Original:** Tf = a0·S + c0 + b·p, con a0 = −0.0575 °C PSU⁻¹, c0 = 0.0901 °C y b = −7.61e-4 °C dbar⁻¹ (`shelfice_thermodynamics.F:184-188` del original).
  - p = `shelficeMass·g·SItodBar` (`:289`), en dbar: el peso del hielo, que con `shelficeLoadAnomaly` equivale a la integral hidrostática de Losch 2008 (`shelfice_init_fixed.F:86-88, 226-227`).
  - T se convierte antes a temperatura in situ con `SW_TEMP` (`:459`; `model/src/seawater.F:88`, UNESCO 1983).
  - `useISOMIPTD` usa una forma no lineal (`:192-196`), distinta también.
- **AD16 (Eq. 25, Tabla 4):** λ1 = −0.0573 °C PSU⁻¹, λ2 = 0.0832 °C y λ3 = −7.53e-8 °C Pa⁻¹ = −7.53e-4 °C dbar⁻¹. Da un punto de fusión **potencial**, sin conversión a in situ (AD16 §3.1.8).
- **Diferencia en Tf** a S = 34.5 y p = 500 dbar: original −2.27415 °C, AD16 −2.27015 °C, es decir 0.004 °C. A esto se suma la conversión in situ del original, del orden de centésimas de °C a 500 dbar [verificar].

### 1.5 Conclusión del Step 1

Hacen falta tres cambios pequeños, todos en **una sola rutina**. No requieren reestructurar nada, así que no hubo punto de decisión.

## 2. Cambios en `config/ocean0/code/`

Las copias empiezan diciendo de qué archivo son copia, y cada cambio está marcado `ISOMIP+ (Exp. 2)`.

| Archivo | Origen | Cambio |
|---|---|---|
| `shelfice_thermodynamics.F` | `pkg/shelfice/` | Bajo `#ifdef SHI_ISOMIP_PLUS_MELT`: (1) a0, c0, b = λ1, λ2, λ3 [°C/dbar] (constantes `PARAMETER` shiLam1/2/3dBar); (2) no se llama a `SW_TEMP`; (3) u*² = Cd·\|u\|² + u_tidal² (shiUTidal = 0.01 m s⁻¹), γT = `SHELFICEheatTransCoeff`·u* y γS = `SHELFICEsaltTransCoeff`·u*, en lugar de la forma HJ99. Requiere `SHELFICEuseGammaFrict=.TRUE.` (las rutas de \|u\| y de `SHIuStar` están dentro de ese bloque). |
| `SHELFICE_OPTIONS.h` | `pkg/shelfice/` | `#define SHI_ISOMIP_PLUS_MELT` |
| `ptracers_apply_forcing.F` | `pkg/ptracers/` | Llama a `SHELFICE_FORCING_MWTR` para el trazador 1. |
| `shelfice_forcing_mwtr.F` | nuevo | Fuente del trazador F_C = −SHIfwFlx/ρfw [m s⁻¹], repartida en la vertical con los mismos pesos que `SHELFICE_FORCING_S` (`pkg/shelfice/shelfice_forcing.F:180-206`). |
| `CPP_OPTIONS.h` | `model/inc/` | `#define NONLIN_FRSURF`, como `verification/isomip/code/CPP_OPTIONS.h:109`. |
| `DIAGNOSTICS_SIZE.h` | `verification/isomip/code/` | numDiags = 30·Nr |
| `PTRACERS_SIZE.h` | `pkg/ptracers/` | PTRACERS_num = 1 (sin cambio numérico, solo comentario) |
| `packages.conf` | — | gfd, shelfice, ptracers, rbcs, diagnostics |

- **El resto del cálculo es el original:** la rama de tres ecuaciones BRIOS (`shelfice_thermodynamics.F:560-683` del original).
  - Con `SHELFICEkappa = 0` y `SHELFICEadvDiffHeatFlux = .FALSE.`, eps3 = eps8 = 0, es decir, hielo aislante (κi = 0, AD16 Tabla 4).
  - Con `SHELFICEconserve = .TRUE.`, los flujos de calor y sal son los de AD16 Eq. 29–30 (forma conservativa de Jenkins et al. 2001).
  - El gemelo en Python es `theory/melt.py::three_equation`, con una cita de línea por operación.
- **Capa límite de Losch (2008):** `SHELFICEboundaryLayer=.TRUE.`.
  - T, S y u se promedian sobre una celda de espesor bajo el hielo (`shelfice_thermodynamics.F:341-389` del original).
  - Esto fuerza `SHELFICE_oldCalcUStar` (`shelfice_readparms.F:207-208`).
  - Los flujos se reparten en esa misma capa (`shelfice_forcing.F:74-101, 180-206`).
- **ρsw:** el código usa rUnit2mass = rhoConst = ρref = 1027.51 kg m⁻³ para ρsw (AD16: 1028). AD16 §3.1.7 pide usar ρref en ambos casos si el modelo los iguala. Diferencia de tipo (c), elección numérica, de 0.05 %.

### Binarios (Step 4)

Hay un `code/` común (`config/ocean0/code`) y un directorio por resolución con solo `SIZE.h` propio. Los demás archivos son enlaces, como en el Exp. 1 (`conv_dx5`).

| Config | Ejecución | nPx × sNx | sNy |
|---|---|---|---|
| `ocean0_2km` | MPI, np = 5 | 5 × 48 | 40 |
| `ocean0_4km` | serial | 1 × 120 | 20 |
| `ocean0_4km_mpi` | MPI, np = 5 | 5 × 24 | 20 |

- La brief pedía `SIZE.h.mpi`. `scripts/build_mitgcm.sh` usa un solo `-mods`, así que el 4 km MPI va en su propio directorio. `run_campaign.py` ahora rechaza np ≠ nPx·nPy (`pfcommon.campaign.check_np`).
- `ocean0_4km_mpi` no tiene casos en `cases.yaml`. Queda compilado para un caso de 4 km a np = 5, si se necesita.

## 3. Namelists (plantillas en `config/ocean0/input_templates/`)

### EOS, mezcla y límites

- **EOS lineal** (AD16 Eq. 23): rhoNil = rhoConst = 1027.51, tRef = −1, sRef = 34.2, α = 3.733e-5, β = 7.843e-4. Además, cp = 3974 y rhoConstFresh = 1000.
- **Mezcla** (AD16 Tabla 4):
  - νH = 6 y κH = 1 m² s⁻¹; νv = 1e-3 y κv = 5e-5, implícitas;
  - `ivdc_kappa = 0.1`.
  - **Nota:** AD16 pide también ν_unstab = 0.1, pero `ivdc` solo actúa sobre los trazadores (`model/src/calc_3d_diffusivity.F:96-110`). Diferencia de tipo (c).
- **Paredes laterales sin deslizamiento.** Fondo y base del hielo con deslizamiento libre y arrastre cuadrático, Cd = 2.5e-3 (`bottomDragQuadratic`, `SHELFICEDragQuadratic`; AD16 Tabla 4). El arrastre no incluye u_tidal (AD16 §3.1.8).
- **Advección:** scheme 77 para T, S y el trazador (context-config, Yung 2026).

### Superficie libre

- Superficie libre no lineal: `nonlinFreeSurf=4` y `select_rStar=2`, con hFacInf = 0.02 y hFacSup = 2 (`verification/isomip/input.htd/data`).
- `useRealFreshWaterFlux=.FALSE.`: flujo virtual de sal, AD16 Eq. 28–30.

### Celdas parciales: `hFacMin = 0.2`, `useMin4hFacEdges = .TRUE.`

- **Justificación:** de los envíos MITgcm a ISOMIP+, MITgcm-BAS usa celdas parciales de ≥ 20 % y la misma columna mínima de 40 m que nosotros. MITgcm-JPL usa 25 % (Yung et al. preprint 2025, `.md` L887–890; versión final Y26 L763–769). AD16 no fija el valor.
- **Regla de redondeo de MITgcm** (`model/src/ini_masks_etc.F:105-195`): una fracción menor que hFacMin/2 cierra la celda, y una en [hFacMin/2, hFacMin) se redondea a hFacMin. Se aplica primero contra el fondo y luego contra la base del hielo.
- **Lo que hace gendata:**
  1. reproduce esa regla (`geometry.mitgcm_hfac`);
  2. aplica la columna mínima de 2 celdas y la conectividad sobre las fracciones de MITgcm;
  3. escribe `bathy.bin` e `icetopo.bin` ya ajustados a ellas.

  Así, el modelo no cambia la geometría que lee, y `shelficeLoadAnomaly.bin` queda evaluado en la base del hielo que usa el modelo.
  - Pruebas: `test_snapped_geometry_is_a_fixed_point_of_mitgcm_rounding` y `test_load_anomaly_at_snapped_draft`.
- **Hallazgo diferido de la Task 1** (`_connected` trataba cualquier fracción > 0 como abierta): una columna conectada solo por una lámina menor que hFacMin/2 queda ahora seca, en vez de convertirse en un lago aislado en el modelo. Prueba: `test_sliver_connection_is_rounded_away`.
- **Caras:** `useMin4hFacEdges` da hFacW = min(hFacC) (`ini_masks_etc.F:243-262`), que es lo que supone la prueba de conectividad.

### Tiempo y pickups

- Se usan nIter0 y nTimeSteps (sin endTime ni startTime).
- pChkptFreq = 1 año y chkptFreq = 1 mes (rotativo). MITgcm escribe además un pickup permanente al final de cada corrida (`writePickupAtEnd`, valor por omisión .TRUE., `model/src/set_defaults.F:345`; `model/src/do_write_pickup.F:76`).
- `pickup_stems` = [pickup, pickup_ptracers]. No hay pickup de shelfice, porque SHELFICEMassStepping está apagado.
- **Para la Task 4:** las etapas de Ocean1 deben terminar en años enteros (lo exige run_stages con pChkptFreq) y, por los diagnósticos, en límites de mes.

## 4. Trazador de deshielo (TRAC01)

- **Variable:** C = fracción de volumen de agua de deshielo [adimensional]. Empieza en 0.
- **Fuente:** virtual, F_C = −SHIfwFlx/ρfw [m s⁻¹] por unidad de área de la base. Se reparte entre la celda superior mojada (kTopC) y la fracción de kTopC+1 que completa una celda de espesor, con los mismos pesos que la sal (`config/ocean0/code/shelfice_forcing_mwtr.F`).
  - Así, ∫C dV = volumen de agua de deshielo que ha entrado al océano.
  - No se usa la ruta `PTRACERS_EvPrRn` (`pkg/ptracers/ptracers_forcing_surf.F:113-132`): con flujo virtual no actúa.
- **Sumidero:** rbcs restaura C a 0 en la zona de esponja, con la misma máscara y τ que T y S. Sin `relaxPtracerFile`, el valor objetivo RBC_ptracers es 0 (`pkg/rbcs/rbcs_init_varia.F:72`).
- **Mezcla:** la misma advección (77) y la misma difusividad que T y S. `PTRACERS_Iter0 = 0`, para que las etapas ≥ 2 lean el pickup (`pkg/ptracers/ptracers_init_varia.F:96, 124`).

## 5. Diagnósticos (`data.diagnostics`)

| Archivo | Frecuencia | Campos |
|---|---|---|
| `shiDiag` | 2 628 000 s (media mensual) | SHIfwFlx, SHIhtFlx, SHIuStar, SHIgammT, SHIgammS, SHIForcT, SHIForcS, SHI_TauX, SHI_TauY, ETAN (2D) |
| `dynDiag` | 2 628 000 s | UVEL, VVEL, WVEL, THETA, SALT, TRAC01, UTRAC01 (3D, todos los niveles) |
| `snapDiag` | `snapFreq` < 0 (instantáneo; 1 año en `cases.yaml`), timePhase = 0 | SHIfwFlx, SHIuStar, SHIgammT, ETAN (2D) |
| `snap3D` | igual que `snapDiag`, mismo paso | THETA, SALT, TRAC01 (3D) |
| `shiFast` | 432 000 s (5 días) | SHIfwFlx, SHIuStar, ETAN |
| `fewLev` | 432 000 s | THETA, SALT, UVEL, VVEL, TRAC01 en los niveles 5, 15 y 25 (90, 290 y 490 m) |
| `dynStDiag` | 86 400 s | estadísticas globales de ETAN, THETA, SALT y TRAC01 (texto) |

- **Dos reglas de `pkg/diagnostics` aprendidas en el smoke** (primera corrida, ver §7):
  - una lista toma los niveles de su primer campo, así que nunca se mezclan campos 2D y 3D en una lista;
  - la fase por omisión de las instantáneas es −frecuencia/2 (`pkg/diagnostics/diagnostics_readparms.F:295`), por lo que se fija `timePhase = 0`.
- `dumpAtLast=.TRUE.`: una corrida que no termina en límite de mes, como el smoke de 30 días, escribe el promedio parcial.
- Cumple R4: medias mensuales 3D de UVEL, VVEL, THETA, SALT y TRAC01 (más UTRAC01 y ETAN) y los SHI* 2D durante toda la corrida, incluidos los 20 años de Ocean1.
- **Volumen en disco a 2 km**, float64: un campo 3D ocupa 2.76 MB, así que dynDiag ocupa ≈ 19 MB por mes, 0.46 GB por cada corrida de 2 años y ≈ 4.6 GB para Ocean1 (20 años).

## 6. Criterios del smoke (escritos ANTES de lanzar la corrida)

**Corrida:** `cal2_gT0.03`, 2 km, np = 5, 30 días (`--end-time 2592000 --run-suffix _smoke30d`), dt = 200 s (provisional de la Task 1), es decir 12 960 pasos.

**Única desviación del comando del controlador:** `--set snapFreq=-864000`, que da instantáneas a los 10 y 20 días. Con el valor por omisión (1 año), una corrida de 30 días no tendría ninguna instantánea para el criterio S5. Solo cambia la salida: el flujo 3 de `data.diagnostics`.

Las celdas se evalúan con el hFacC del modelo (`hFacC.data`). La plataforma son las columnas mojadas con base del hielo < 0.

- **S1 (= N3):** `require_ok` pasa: estado OK, last_iter = 12 960 = nIter0 + nTimeSteps, sin NaN ni Infinity en el log.
- **S2 (= N2):** |ΔInv − (∫fuente − ∫sumidero)| / ∫fuente < 1 % sobre la corrida.
  - ΔInv = Σ C_fin·hFacC·rStarFacC·drF·rA, con C_fin = `PTRACER01` del último paso (`dumpInitAndLast`) y rStarFacC = (η + Hcol)/Hcol (`model/src/calc_r_star.F:103-105`). El inventario inicial es 0.
  - ∫fuente = −Σ ⟨SHIfwFlx⟩·rA·T/ρfw, con la media de 30 días de `shiDiag` y T = 2 592 000 s.
  - ∫sumidero = Σ (máscara/τ)·⟨TRAC01⟩·hFacC·drF·rA·T, con la media de 30 días de `dynDiag`. Es exacto porque el término es lineal en C.
  - **Informativo:** la misma razón sin el factor rStarFacC.
- **S3:** T y S dentro de [mín, máx] del perfil inicial ± 0.1: T ∈ [−2.0, 1.1] °C y S ∈ [33.7, 34.8]. Se evalúa en todas las celdas mojadas del estado final (`T`, `S` de `dumpInitAndLast`), de la media de 30 días (`dynDiag`) y de las dos instantáneas.
- **S4:** ⟨m⟩ > 0, con ⟨m⟩ = media ponderada por área, sobre la plataforma, de m = −⟨SHIfwFlx⟩/ρfw·31 536 000 en la media de 30 días.
- **S5:** el melt de MITgcm coincide con `theory.melt.three_equation` con error < 1 % relativo en todas las celdas de plataforma con m > 1 m a⁻¹, en cada una de las dos instantáneas (`snapDiag`, 10 y 20 días).
  - Entradas: T y S de la capa límite, reconstruidos de THETA y SALT instantáneos con los pesos de Losch de `shelfice_thermodynamics.F:341-389` (original) y el hFacC del modelo; p = `shelficemassinit`·g en Pa (`shelfice_thermodynamics.F:289`); u* = SHIuStar instantáneo; ΓT = 0.03.
  - Se reportan el máximo y la mediana del error. Pasa si el máximo es < 1 %.
- **Informativo (sin umbral):** el hFacC del modelo contra `geometry.mitgcm_hfac` de la geometría escrita; el tiempo por paso [s/paso] (wall_s/12 960) y la extrapolación a 2 años.

Si un criterio falla, no se cambia: se diagnostica, se reporta y se para para que decida el controlador.

## 7. Resultados del smoke

### Corridas

Todas en `/Volumes/LaCie/portfolio_runs/exp2/`.

1. **Intentos 1 y 2:** errores nuestros de namelist, sin ningún paso de tiempo.
   - `useMin4hFacEdges` está en PARM04, no en PARM01 (`model/src/ini_parms.F:302`).
   - `delZ` no se puede mezclar con `viscAr`/`diffKr*`, que son nombres en r ("Cannot mix z, p and r"); se usa `delR`.
   - Corregidos, con pruebas en `tests/test_exp2_config.py`.
2. **Corrida 1** (`cal2_gT0.03_smoke30d_run1_mixedsnap`): OK, pero con `snapDiag` defectuoso. Mezclaba campos 2D y 3D, y guardó solo el nivel 1 de THETA/SALT/TRAC01. Además, las instantáneas salieron a los 5, 15 y 25 días por la fase por omisión. Con eso no se puede evaluar S5.
   - Se separó en `snapDiag` (2D) y `snap3D`, con `timePhase = 0`. La física es idéntica: el T final es igual bit a bit al de la corrida 2.
3. **Corrida 2** (`cal2_gT0.03_smoke30d`): la que se evalúa. Instantáneas a los 0, 10 y 20 días.
   - Evaluación: `exp2_isomip_cavity/analysis/pd1_smoke_check.py`. Resultado completo en `pd1_smoke_check.json`, en el directorio de la corrida.

### Veredictos (corrida 2)

| Criterio | Resultado | Valor |
|---|---|---|
| S1 (N3) | **PASA** | OK; 12 960/12 960 pasos; sin NaN |
| S2 (N2) | **PASA** | error relativo 1.9e-15. Inventario = 1.519e10 m³, fuente = 1.519e10 m³, sumidero rbcs = 8e-111 m³ (en 30 días el trazador no llega a la esponja). Sin el factor r*: 7.7e-5 (informativo) |
| S3 | **FALLA** | T ∈ [−2.52, 0.96] °C (límite [−2.0, 1.1]) y S ∈ [33.35, 34.69] (límite [33.7, 34.8]), en el estado final. La media de 30 días da T ≥ −2.40 y S ≥ 33.40; las instantáneas de 10 y 20 días, T ≥ −2.51 y S ≥ 33.36. El máximo nunca se pasa |
| S4 | **PASA** | ⟨m⟩ de plataforma = 17.2 m a⁻¹ (2689 columnas). Informativo: melt profundo (base < −300 m) = 53.9 m a⁻¹ en el primer mes, todavía en el ajuste inicial |
| S5 | **PASA** | 10 días: 2301 celdas con m > 1 m a⁻¹, error relativo máximo 0.29 %, mediana 3.7e-5. 20 días: 1941 celdas, máximo 0.17 %, mediana 5.6e-5. (t = 0: 2689 celdas, 9e-15) |

**Informativo:**
- El hFacC del modelo es igual a `geometry.mitgcm_hfac` de la geometría escrita (diferencia máxima 0.0).
- Tiempo: 247 s para 12 960 pasos, es decir 0.0191 s/paso a np = 5. Una corrida de 2 años a 2 km tardaría ≈ 1.67 h y Ocean1 (20 años) ≈ 16.7 h, en etapas de ≤ 8 h que terminan en años enteros.
- dt = 200 s se mantiene: CFL advectivo máximo de 0.084 (u), 0.071 (v) y 0.11 (w); \|u\|máx = 0.84 m s⁻¹ durante el ajuste; ≤ 127 iteraciones de cg2d.

### Diagnóstico de S3 (no se cambia el criterio; decide el controlador)

- **Dónde:** 5001 de 164 922 celdas mojadas (3 %), en el estado final. 1723 tienen T < −2.0 y 4932 tienen S < 33.7.
  - 3925 están bajo el hielo, sobre todo en las 1–2 celdas superiores (2612 en la celda superior); 1076 en mar abierto, en la pluma que sale.
  - **Todas tienen agua de deshielo:** C ≥ 0.0055, mediana 0.017, máximo en el dominio 0.034.
- **Magnitud coherente con mezcla de agua de deshielo** (línea de Gade):
  - con C ≈ 0.034, el flujo virtual da ΔS ≈ −S·(ρfw/ρ0)·C ≈ −1.1 PSU y ΔT ≈ −(L/cp + T − Tb)·(ρfw/ρ0)·C ≈ −2.9 °C;
  - partiendo de ≈ (0.5 °C, 34.5) en la base profunda, esto da ≈ (−2.4, 33.4), que es lo observado;
  - correlaciones con C: −0.70 (ΔS) y −0.46 (ΔT). Las pendientes ajustadas, −22 y −47, son menores que las de mezcla pura, porque ΔT y ΔS se toman contra el perfil inicial en la misma celda y las parcelas fueron advectadas.
- **Clasificación:** (b)/(c). El rango del perfil ambiental ± 0.1 no admite el agua modificada por el deshielo, que el protocolo espera bajo el hielo; el punto de fusión a 500 m ya es ≈ −2.2 °C. No es un error nuestro identificable.
- **Hallazgo aparte** (preocupación, no criterio): 1589 celdas están por debajo del punto de fusión local, Tf(S, p), con un mínimo de −0.31 °C, y 290 de ellas por más de 0.1 °C.
  - La más fría es una celda superior delgada (hFacC = 0.2, 4 m) a 496 m: T = −2.52 frente a Tf = −2.21.
  - Hipótesis, por verificar: con la capa límite de Losch, el melt se calcula con T promediado sobre 20 m (la celda delgada más 16 m de la siguiente) y el enfriamiento se reparte por igual. La celda delgada puede quedar sobreenfriada mientras el promedio sigue sobre el punto de fusión. Además, el agua que asciende por la base pierde presión. El modelo no tiene frazil.
