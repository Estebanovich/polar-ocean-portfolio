# PD1 — Configuración `fjord_vmf`, cronometraje y presupuesto (Task 4)

Fecha: 2026-11-12. Material para el PD1. Esteban aprueba aquí la configuración y el presupuesto, y **la campaña no se corre hasta entonces**. Esta nota no elige por él: donde hay opciones, da las medidas y, aparte y marcada como tal, una recomendación técnica.

- **Notas relacionadas:** `notes/pd1_seaice.md` (lectura de `pkg/seaice`, `pkg/kpp` y `pkg/ggl90`), `notes/pd1_batimetria.md` (Task 1b), `notes/pd1_puente.md` (Task 3) y `notes/pd0_preregistro.md` (congelado).
- **Números:** `figures/pd1_timing.json`, que escribe `analysis/pd1_timing.py` a partir de las corridas cortas.
- **Corridas:** `/Volumes/LaCie/portfolio_runs/exp4/<caso><sufijo>/`. **Builds:** `/Volumes/ESTEBAN_MIT/portfolio_builds/exp4/<config>/`.

---

## 0. Decisiones del PD1 (Esteban, 2026-11-12)

1. **Paso de tiempo:** `staggerTimeStep = .TRUE.` y **Δt = 120 s** en todas las corridas del fiordo (F1–F5, F1c, F1g, V3, V4); `cases.yaml` y la plantilla `data`.
   - Es un cambio numérico, clase (c). Motivo: el aviso de `gad_check.F:94-107` sobre las ondas internas con un esquema de trazador sin AB (§2).
   - **Comprobación con la configuración final** (escalonado, 120 s, `tracForcingOutAB = 1`, copias modificadas de `seaice` y `obcs` compiladas; `F1_pd1final_48h`, diagnóstico): **limpia**. 1440 pasos en 888 s (0.62 s/paso); máx |w| = 1.6×10⁻³ m s⁻¹, que baja a 6.0×10⁻⁴ a las 48 h, por debajo del umbral *post hoc* de 5×10⁻³ de §2; máx CFL 0.129 (W_hf, a las 4.3 h). `STDERR.0000` vacío: ya no aparece el aviso de `gad_check`.
2. **Hielo:** la opción 2 de `pd1_seaice.md` §6 (`DWATN` también sin reología) y la opción 1 de §7 (`obcs_check.F` relajado con `DISABLE_SEAICE_OBCS`), como copias de MITgcm en `config/fjord_vmf/code`.
   - F4 usa el mismo polígono que F2 (ruling R15).
   - Verificación en `pd1_seaice.md` §8: F2 y F4 arrancan con Orlanski; el esfuerzo bajo A = 1 coincide con ρ₀ C_w|Δu|Δu con el piso `SEAICEdWatMin` (F2 exacto; F4 ≤ 2.7 % del máximo de |τ| en el interior del hielo); el reinicio es idéntico bit a bit.
   - Comprobación obligatoria tras F2 y F4 (Task 6): el hielo a ≥ 10 celdas de toda frontera abierta durante las 280 h (`pd1_seaice.md` §8).
   - **Binario de la campaña (Task 6).** El ejecutable `fjord_vmf` con que corren F1–F5, F1c, F1g y V3 se compiló el 2026-11-12 a las 07:15, al arrancar F1, con el código del commit **afc3526** (las cabeceras de las copias modificadas cambiaron en ese commit; `config/fjord_vmf/code` no cambió después: `git diff afc3526 HEAD` vacío). sha256 `25bdc7f3e89e674f1653696f2792bc2057cf5703b541fb44b1c47066b10fed91`, el mismo que registra `binary_sha256` en el `case.json` de F1, F2, F3, F1c y V3. No se recompila durante la campaña.
     - `run_campaign.py --dry-run` (que ahora informa del plan; Task 6) confirma que **F1g, F5 y F4 no recompilan** (build al día con `--skip-ok`).
     - **V4 sí recompila** `fjord_vmf_400`: ese ejecutable es de las 01:48, anterior a las copias modificadas de 8ba509d/afc3526 (`OBCS_OPTIONS.h`, `SEAICE_OPTIONS.h`, `obcs_add_tides.F`, `obcs_check.F`, `seaice_dynsolver.F`). La recompilación va a otro directorio (no toca el binario de `fjord_vmf`) y deja V4 con el mismo código que F1, que es lo que pide V4; queda registrada como aviso para el controlador.
3. **T_fin = 280 h** en todas las corridas, sin recortes.
   - El presupuesto de cómputo sube a **~31–32 h en total**: Task 3 6.6 h, cronometraje y verificación de esta Task ~1.5 h, y campaña con Δt 120 s ≈ 23.0 h (9 corridas del fiordo 12.4 h + R2/R3 10.6 h; §3).
   - El límite de 8 h por corrida se cumple: la más larga del fiordo es F1g, ≈ 1.8 h.
4. **V3:** forzamiento suave con Δt 120 s escalonado: η₀ = **0.1 m** (U₀ = 0.34 cm s⁻¹) y rampa de medio coseno durante un ciclo M2 (12.42 h; `obcs_add_tides.F` modificado, `EXP4_TIDE_RAMP_S`).
   - **Estabilidad en 2 ciclos M2** (corridas de diagnóstico de 25 h, 750 pasos):

     | η₀ | Corrida | Resultado | Máx CFL horiz. | Máx CFL vertical (W_hf) | Máx corriente (Mariasundet / Akselsundet) |
     |---|---|---|---|---|---|
     | 0.2 m | `V3_pd1final_eta02_25h` | **diverge** en el paso ~640 (21.3 h, 1.7 ciclos) | 0.44 | 1.14 antes del NaN | 0.92 m s⁻¹ y creciendo |
     | 0.15 m | `V3_pd1final_eta015_25h` | completa los 2 ciclos | 0.41 | **0.93** (paso 480, 16 h) | 0.76 / 0.67 m s⁻¹ |
     | **0.1 m** | `V3_pd1final_eta01_25h` | completa los 2 ciclos | 0.27 | **0.548** (paso 480), ≤ 0.38 el resto | 0.56 / 0.45 m s⁻¹ |

   - **Elección.** Las dos amplitudes menores completan los 2 ciclos. Ninguna cumple estrictamente el CFL ≤ 0.5 del plan, pero 0.15 m llega a 0.93, cerca del 1.14 con el que divergió 0.2 m, y 0.1 m se queda en 0.55. Tomo **0.1 m** por ser la más cercana al criterio de CFL, con menos riesgo en 22 ciclos. Si se prefiere la lectura literal ("la mayor que completa 2 ciclos"), 0.15 m también está probada y es un cambio de una línea en `cases.yaml`.
   - **La estabilidad más allá de 2 ciclos no está probada.** La corrida de 280 h tiene ~22 ciclos M2: es un riesgo.
   - **Regla de respaldo fijada de antemano (ruling R18):**
     - si V3 diverge en la campaña, se repite una vez con el siguiente η₀ más bajo probado o la mitad del elegido (con 0.1 m: **0.05 m**), y se declara;
     - si vuelve a diverger, V3 se reporta como validación no realizada, declarada.
     - V3 no entra en ningún cociente de H4.
   - **La amplitud no es realista** (S&F14: 0.3–0.8 m fuera de Akseløya, p. 89). Con η₀ = 0.5 m V3 divergía por las corrientes de los estrechos (§1).
   - V3 mide la velocidad de fase de la onda de Kelvin interna frente al umbral congelado del PD0 (§4: 0.21–0.84 m s⁻¹, alrededor de 0.42 m s⁻¹, sin cambios), además de la presencia y el sentido ciclónico. La amplitud no es un criterio.
5. **Lo que queda como estaba** (regla del PD0):
   - el invierno sigue a 200 m, con la limitación declarada (R₁/Δx = 1.7–6.6; opción (c) de §5);
   - caja de 100 × 60 km, 200 m y 40 niveles;
   - Smagorinsky 2.2 + `viscC4smag` 1;
   - los pendientes de batimetría tal como están (§6, punto 10).

---

## 1. Configuración generada (los valores del PD0, sin cambios)

**Un solo build** (`config/fjord_vmf/code`):
- paquetes `gfd obcs diagnostics kpp ggl90 seaice`;
- sin `mnc` y sin `exf`; la termodinámica del hielo está apagada (`pd1_seaice.md` §4);
- `data.pkg` elige el cierre y el hielo.

Los builds de 400 m, de 100 m y de la subcaja de 100 m (`fjord_vmf_400`, `fjord_vmf_100`, `fjord_vmf_100s`) solo cambian `SIZE.h`. Una prueba exige que los demás archivos sean idénticos.

**Desviación del brief:** el brief nombra `SIZE.h` y `SIZE.h.mpi`. Seguí el patrón del repositorio: un `SIZE.h` por configuración (como `bay_thesis` y `exp2/ocean0_4km_mpi`), 5 rangos (`nPx = 5`, franjas a lo largo de x). No hay `SIZE.h.mpi`, porque el build no lo leería.

Lo que fija el PD0 sale de `cases.yaml` (`configs.fjord_vmf`) y de las plantillas `config/fjord_vmf/input_templates`:

| | Valor | Fuente |
|---|---|---|
| Caja | UTM 33N, x 455–555 km, y 8594–8654 km, 500 × 300 × 40 | Task 1b, ruling R5 |
| Malla vertical | 40 niveles, 1.5 m arriba, Σ = 240 m; hFacMin 0.2 y hFacMinDr 2 m | R1, R2 |
| f | 1.4266×10⁻⁴ s⁻¹, plano f (β = 0) | plan |
| EOS | lineal ajustada (tα 3.0456×10⁻⁵, sβ 8.4068×10⁻⁴, ρNil = ρConst 1026.96, tRef 0.4676, sRef 33.6395) | PD0 §1.3 |
| Cierre | KPP, viscAr 1×10⁻⁴, diffKr 1×10⁻⁵, `implicitViscosity`, sin `ivdc_kappa` | PD0 §6 |
| Horizontal | Smagorinsky 2.2 y `viscC4smag` 1; viscAh constante 0; diffKh 0 | PD0 §6 |
| Fondo y advección | `bottomDragQuadratic` 2.5×10⁻³; esquema 77 para T y S | PD0 §6 |
| F1c | KPP apagado, viscAr 1×10⁻⁵, `ivdc_kappa` 1 (el cierre de la tesis, `data.stage1`); lo demás igual que F1 | plan; prueba `test_f1c_thesis_closure` |
| F1g | GGL90 con sus valores por omisión en lugar de KPP | plan |
| Viento | pulso de la tesis (máx. 0.1 N m⁻², σ 1 h, centro 2.5 h), uniforme, **hacia 192.412°** desde el este (es decir, del 077.6°, el "078°" del PD0); F5, +90° | PD0 §7 |
| Hielo | h = 1 m y A = 1 en el polígono del fiordo interior (el de los perfiles): 13 286 celdas, 531.4 km² | spec §3.5 |
| C_w | 5.5×10⁻³ | PD0 §8 |
| Fronteras | O, N y S Orlanski (CMAX 0.45 y `cvelTimeScale` 250 s, de la tesis); esponja de 10 celdas en las velocidades; balance barotrópico | plan |
| Salida | §4 | — |

**Detalles que el PD0 no fijaba y que fijé aquí** (se pueden cambiar en el PD1):
- **Forzamiento.** Los registros del viento van en los instantes que usa MITgcm, (n − ½)·1200 s (`eesupp/src/get_periodic_interval.F:118-128`). Así el pulso del modelo es la forma nominal muestreada exactamente: el pico de 0.1 N m⁻² cae en el registro 8, a 9000 s.
  - La tesis muestreaba `linspace(0, 60 h, 180)` (Task 3).
  - Hay 181 registros (60 h + 1). El último está a cero y el ciclo periódico no vuelve a interpolar hacia el registro 1 dentro de la etapa forzada. Con 180 registros, a las 60 h volvería un viento de hasta 0.003 N m⁻². En la tesis pasa lo mismo, con 0.002 N m⁻² (la mitad del registro 1, 0.0044 N m⁻²).
- **Etapas.** Etapa forzada de 60 h (plan) y etapa libre. `data.free` es `data` sin el viento (prueba `test_free_stage_has_same_physics_without_wind`). Hay pickup permanente cada 60 h.
- **Tiempos de relajación de la esponja:** 1200 s en el borde interior (tesis) y max(60 s, 2Δt) en la frontera. En la tesis, 60 s con Δt = 30 s.
- **V3 (solo M2).** Velocidad normal barotrópica M2 hacia adentro, uniforme en las tres fronteras abiertas (`ALLOW_OBCS_TIDES`, `pkg/obcs/obcs_add_tides.F`), con forma u = U₀ sin ωt. U₀ sale del balance de volumen de una cuenca que sube η₀: U₀ = ω η₀ A_mar / A_frontera (A_mar = 2261 km², A_frontera = 9.37 km²). Primero era η₀ = 0.5 m (U₀ = 1.70 cm s⁻¹); **tras el PD1 es η₀ = 0.1 m (U₀ = 0.34 cm s⁻¹), con una rampa de un ciclo M2 (§0, punto 4).**
  - El η₀ = 0.5 m original estaba dentro de los 0.3–0.8 m de S&F14 (p. 89); el de 0.1 m no lo está, y se declara (§0).
  - **Punto PD1-V3:** es un diseño mínimo. La amplitud real en Akseløya se mide y se reporta, sin calibrarla.
  - **Falla en la corrida corta (`V3_smoke2h`, 2 h, Δt 120 s, escalonado): la marea no actuó.**
    - Máx η = 7.5×10⁻⁶ m a las 2 h. η medio ≈ 10⁻²² m y velocidad normal media en la frontera O ≈ 10⁻⁹ m s⁻¹ (`V3_smoke2h/STDOUT.0000:2337-3013`): flujo neto por las fronteras exactamente cero, cuando la teoría daba η ≈ 0.23 m.
    - **Diagnóstico** (fix round 1; corridas de 1 h, Δt 120 s, escalonado, `V3_diag_nobalance_1h` y `V3_diag_noorlanski_1h`): sin `useOBCSbalance`, η no responde (η medio 5×10⁻⁹ m). Sin Orlanski (fronteras con valor 0 más la marea), la velocidad normal en la frontera es **exactamente 0.0**: la marea no se sumaba a nada.
    - **Causa: error nuestro (a) en las entradas generadas, no Orlanski.** `pkg/obcs` lee por punto de frontera una **amplitud** (m s⁻¹, `OB*_*TidAmFile`) y una **fase en segundos** (`OB*_*TidPhFile`) (`pkg/obcs/OBCS_PARAMS.h:162`), y forma Cs = Am cos(ωPh) y Sn = Am sin(ωPh) (`pkg/obcs/obcs_init_fixed.F:536-575`). Yo escribía U₀ en el archivo de fase y no daba amplitud, así que Am = 0.
    - **Corregido** en `make_case.tide_v3` y `cases.yaml`: Am = ±U₀ y Ph = T/4, de modo que u = U₀ sin ωt. Prueba: `test_v3_tide_only`, que reconstruye u(t) como lo hace MITgcm.
  - **V3 corregido (`V3_diag_tidefix_*`, solo diagnóstico): la marea actúa, pero la corrida es inestable.**
    - η medio sigue el balance de volumen: 0.169 m a 1.67 h, frente a η₀(1 − cos ωt) = 0.167 m.
    - Pero las corrientes de marea en Mariasundet llegan a 1.3 m s⁻¹ a las 2 h (S&F14 citan > 2 m s⁻¹ en Akselsundet, p. 89).
    - Con **Δt 120 s** la corrida diverge (NaN entre 1.7 y 2 h).
    - Con **Δt 60 s**, el CFL vertical (W_hf) llega a 1.01 a las 2 h y a 4.1 justo antes de diverger, entre 2.7 y 3.3 h, cerca del máximo de entrada (T/4 = 3.1 h).
  - **Punto PD1-V3 — decidido: forzamiento más suave con Δt 120 s (§0).** Antes del PD1 la nota decía: V3 necesita un Δt propio mucho menor (≲ 20–30 s para CFL_W ≤ 0.5; costo ≈ 4.7–7 h por corrida de 280 h) o un forzamiento más suave (η₀ menor, rampa), además de la duración y la amplitud ya pendientes. El PD0 no fija un umbral de V3 que dependa de Δt, pero el plan pide el mismo Δt por malla en las corridas comparables, y V3 no se compara con las F en ningún cociente.
- **F4:** el hielo a la deriva usa el mismo polígono que F2, así que F4 frente a F2 aísla la dinámica del hielo. El spec dice "hielo libre en todo el fiordo" y el plan "fuera no hay hielo, salvo en F4". **Punto PD1-F4:** ¿mismo polígono o todo el fiordo, incluido el exterior? Además F4 tiene el problema de `pd1_seaice.md` §2.

---

## 2. Estabilidad y Δt: hay una inestabilidad numérica

**Resultado principal del cronometraje.** Con la configuración del PD0 (sincrónica, esquema 77), **F1 desarrolla ruido de escala de malla** en la plataforma sur y en Bellsund. Empieza al pie del talud sur de Bellsund (i = 62, j = 32; 217 m; pendiente ≈ 0.15), unas horas después del pico del pulso, y se extiende.

Qué es:
- w de varios cm s⁻¹ en damero, en toda la columna;
- U con estructura de modo 1;
- ruido de 0.2 °C en T (`F1_t6h`).

No es física resuelta. Lo que se observa:

| Δt (s) | Corrida | Inicio (máx \|w\| > 5×10⁻³ m s⁻¹)* | inicio × Δt (h·s) | máx CFL_W (W_hf) |
|---|---|---|---|---|
| 60 | `F1_t6h` | 2.67 h | 160 | 1.23 |
| 48 | `F1_diag_dt48` | 3.33 h | 160 | 0.66 |
| 40 | `F1_diag_dt40` | 4.33 h | 173 | 0.38 |
| 30 | `F1_dt30_48h` (abortada en 6.7 h) | 6.33 h | 190 | 0.11 (sube) |
| 60, cierre de la tesis (F1c) | `F1c_diag_4h` | igual que KPP | — | — |

\* El umbral es *post hoc*, fijado después de ver las corridas: la respuesta suave no pasa de 1.5×10⁻³ m s⁻¹.

**Lectura:**
1. No depende de KPP: F1c, con el cierre de la tesis, crece igual.
2. Reducir Δt no la quita, solo la retrasa. El inicio es ∝ 1/Δt (inicio × Δt casi constante), así que la tasa de crecimiento es ∝ Δt. Es la firma de una inestabilidad débil e incondicional del paso de tiempo, no de un límite CFL.
3. **El mecanismo lo señala el propio MITgcm.**
   - `pkg/generic_advdiff/gad_init_fixed.F:148-157` aplica Adams–Bashforth a la tendencia de T y S solo con los esquemas 2, 3 y 4 (centrado de 2.º orden, upwind de 3.er orden y centrado de 4.º orden; `GAD.h:23-33`). El 77 (y el 33) no usan AB.
   - Con paso sincrónico, sin `implicitIntGravWave` y con un trazador activo sin AB, `gad_check.F:94-107` avisa "potentially unstable time-stepping (Internal Wave)" y "need staggerTimeStep=.TRUE.".
   - El aviso está en `F1_t6h/STDERR.0000:1-2` (y en todas las corridas sincrónicas) y no aparece en las escalonadas (`F1_diag_st60_12h/STDERR.0000` vacío).
   - Es la inestabilidad de las ondas de gravedad internas cuando el trazador avanza hacia delante sin AB mientras el momento usa AB2. La documentación recomienda el paso escalonado: "essentially allows the gravity wave terms to leap-frog in time giving second order accuracy and more stability" (`doc/algorithm/algorithm.rst:705-715`).

**Prueba.** F1 con `staggerTimeStep = .TRUE.` y Δt = 60 s, 12 h: **sin ruido** (`F1_diag_st60_12h`).
- Máx |w| se queda en 1.0–1.7×10⁻³ m s⁻¹ durante las 12 h. En la sincrónica, a las 4 h ya valía 6.5×10⁻².
- Máx CFL 0.066, a 0.50 s por paso (el mismo costo).
- Coincide con lo que predice el aviso de `gad_check.F`.
- Prueba con Δt = 120 s escalonado: limpia, máx |w| ≤ 1.6×10⁻³ m s⁻¹ (`F1_diag_st120_12h`). Corrida de 48 h escalonada a **120 s** (el "CFL ≤ 0.5 en 2 días" del plan; `F1_st120_48h`): limpia. Máx |w| ≤ 1.6×10⁻³ m s⁻¹, que baja a 6×10⁻⁴ a las 48 h; máx CFL advectivo 0.129, alcanzado en el paso 130 (4.3 h, pulso).

**Punto PD1-Δt — decidido: opción (a), `staggerTimeStep` con Δt 120 s (§0).** Opciones medidas:
- (a) **[ELEGIDA]** `staggerTimeStep = .TRUE.` con **Δt = 120 s**: limpia en 48 h, CFL advectivo 0.13 y de ondas internas 0.42 (§3). La alternativa conservadora es Δt = 60 s (limpia en 12 h; cuesta ~60 % más por hora simulada). Es un cambio numérico, no físico, y se declara (c). Antes de la campaña conviene repetir la corrida de 48 h con la configuración final y el `tracForcingOutAB = 1` que `pkg/seaice` exige; las de esta Task se hicieron antes de añadirlo, salvo `F1_st120_48h`, `V4_smoke2h`, `F3_100s_t1h` y los diagnósticos de hielo.
- (b) Paso sincrónico con Δt más corto. **No sirve:** con 30 s también aparece a las 6 h.
- (c) **El esquema 33 (alternativa del PD0 §6) no lo arreglaría:** tampoco usa AB (`gad_init_fixed.F:148-157`) y daría el mismo aviso.
- (d) Un esquema de trazador con AB (2, 3 o 4: centrado de 2.º orden, upwind de 3.er orden o centrado de 4.º orden). Quita la condición del aviso, pero son esquemas sin limitador, que el PD0 §6 descartó en favor del 77. Exige enmienda del PD0 y no está probado.
- (e) `implicitIntGravWave = .TRUE.`. Exige compilar con `ALLOW_NONHYDROSTATIC` (`doc/getting_started/getting_started.rst:2141`) y una enmienda del PD0, y no está probado.

**CFL ≤ 0.5:** ver §3.

---

## 3. Cronometraje (5 rangos, núcleos P)

Todas las corridas usan np = 5, un caso a la vez, y escriben a LaCie. `wall_s` (de `run_status.json`) incluye el arranque y la lectura de entradas; "bucle" es el temporizador `MAIN LOOP` de MITgcm. Todas terminaron sin procesos sueltos (`pgrep` vacío).

| Corrida | Malla | Paso | Δt (s) | Pasos | wall (s) | s/paso (wall) | s/paso (bucle) | Máx CFL adv. | Estado |
|---|---|---|---|---|---|---|---|---|---|
| `F1_t6h` | 200 m | sincr. | 60 | 360 | 187 | 0.52 | 0.51 | 1.23 (ruido) | OK, con ruido |
| `F1_diag_dt48` / `dt40` / `dt30` | 200 m | sincr. | 48 / 40 / 30 | 450 / 540 / 480 | 220 / 253 / 219 | 0.49 / 0.47 / 0.46 | — | 0.66 / 0.38 / 0.03 | OK, con ruido salvo la de 30 s (4 h) |
| `F1_dt30_48h` | 200 m | sincr. | 30 | 800 de 5760 | 363 | 0.45 | — | 0.11 | abortada a las 6.7 h (ruido) |
| `F1c_diag_4h` | 200 m | sincr. | 60 | 240 | 116 | 0.48 | — | 0.97 (ruido) | OK, con ruido |
| `F1_diag_st60_12h` | 200 m | escal. | 60 | 720 | 363 | 0.50 | 0.50 | 0.066 | OK, limpia |
| `F1_diag_st120_12h` | 200 m | escal. | 120 | 360 | 222 | 0.62 | — | 0.129 | OK, limpia |
| `F1_st120_48h` | 200 m | escal. | 120 | 1440 | 881 | 0.61 | 0.61 | 0.129 | OK, limpia en 48 h |
| `F1g_smoke1h` | 200 m | escal. | 60 | 60 | 43 | 0.72 | 0.63 | 0.011 | OK |
| `V3_smoke2h` | 200 m | escal. | 120 | 60 | 46 | 0.77 | 0.69 | 0.001 | OK, pero sin marea por un error de entradas (§1; corregido) |
| `V3_diag_tidefix_2h` / `_dt60_4h` | 200 m | escal. | 120 / 60 | 60 / 240 | 47 / 133 | — | — | NaN | divergen (marea corregida; §1) |
| `V4_smoke2h` | 400 m | escal. | 120 | 60 | 13 | 0.22 | 0.17 | 0.048 | OK |
| `F2_smoke3h`, `F4_smoke3h` | 200 m | escal. | 60 | 0 | — | — | — | — | no arrancan: `OBCS_CHECK` (`pd1_seaice.md` §7) |
| `F3_100s_t1h` (subcaja 100 m) | 100 m | escal. | 60 | 60 | 87 | 1.45 | 1.32 | 0.039 | OK (ejecutable re-enlazado; §5) |

- **Memoria por rango (200 m):** RSS ≈ 0.90 GB (medido con `ps` durante `F1_t6h`). Datos estáticos del binario (`size -m`, segmento `__DATA`): 1.63 GB.
- **Cuánto pesa el Δt.** A 120 s el paso cuesta un 24 % más que a 60 s: 0.62 frente a 0.50 s, probablemente por más iteraciones de `cg2d`. Aun así, cada hora simulada cuesta un 38 % menos.

**Extrapolación a la campaña** (s/paso de corridas de ≥ 12 h; los casos que solo tienen una corrida corta usan el s/paso de F1 a su Δt, más un 25 % para F1g, como en la medida de 1 h):

| | Δt 60 s escal., 280 h | Δt 120 s escal., 280 h | Δt 120 s escal., 240 h |
|---|---|---|---|
| Por corrida a 200 m (F1–F5, F1c, V3) | 2.3 h | 1.45 h | 1.24 h |
| F1g | 2.9 h | 1.8 h | 1.6 h |
| V4 (400 m) | 0.8 h | 0.4 h | 0.35 h |
| **9 corridas del fiordo** (7 × la fila 1 + F1g + V4) | **20.0 h** | **12.4 h** | **10.6 h** |
| R2 + R3 (Task 3: 5.3 h cada una a 280 h) | 10.6 h | 10.6 h | 9.1 h |
| **Total de la campaña** | **30.6 h** | **23.0 h** | **19.7 h** |

- **Presupuesto (spec §4):** ≤ 8 h por corrida y ≤ 30 h en total. **Tras el PD1 (§0, punto 3) el total sube a ~31–32 h y se corre todo a 280 h.**
  - La Task 3 gastó 6.6 h (R1, R1k) y esta Task ~0.9 h de modelo (suma de `wall_s`, más los builds). Quedan ~22 h.
  - **Con Δt 120 s escalonado:** 280 h se pasa por ~1 h (23.0 h); 240 h cabe (19.7 h).
  - **Con Δt 60 s:** 280 h (30.6 h) no cabe sin recortar (F4, F5, V4, F1g).
- **Si el PD1 elige 280 h con Δt 120 s,** recortar F4 (que además no tiene sentido físico sin la opción 2 de `pd1_seaice.md` §6) deja 21.5 h, justo dentro.

**CFL ≤ 0.5 (criterio del plan, en 2 días).**
- Advectivo: con Δt 120 s escalonado, el máximo en 12 h es 0.13 (W_hf). En 48 h: 0.129 (también W_hf, a las 4.3 h).
- **Ondas internas:** c₁Δt/Δx con c₁ = 0.70 m s⁻¹ (modo 1 de verano en la columna más profunda, 238.8 m; `theory/modes.py` con el perfil y la EOS ajustada) da 0.21 con 60 s, **0.42 con 120 s**, 0.52 con 150 s y 0.84 con 240 s.
- Por tanto, **120 s es el mayor Δt que divide los 1200 s de los registros y deja los dos CFL ≤ 0.5**, y solo con paso escalonado.
- En invierno c₁ = 0.29 m s⁻¹ en la columna más profunda; a 100 m, el número de ondas internas con 120 s sería 0.34.

---

## 4. Salida en disco por corrida (ruling R14)

**Diagnósticos** (MDS, float64, instantáneas cada 1200 s salvo indicación; `config/fjord_vmf/input_templates/data.diagnostics` y los conjuntos de `cases.yaml`). Bytes por registro medidos en `F1_t6h`:

| Flujo | Campos | Para qué | MB/registro (200 m) |
|---|---|---|---|
| `state3d` | UVEL, VVEL, RHOAnoma (3-D) | flujo de Nash et al. 2005 en los estrechos; Ē_in; espectros en estaciones | 144 |
| `state2d` | ETAN | espectros y mapas de η | 1.2 |
| `wlev` | WVEL en k = 7, 18, 26 (9.8, 30.3 y 61.3 m) | mapas PSD de w | 3.6 |
| `surf` | oceTAUX, oceTAUY (después de `seaice`) | trabajo del viento W_NI | 2.4 |
| `kpp2d` (KPP) | KPPhbl | capa límite | 1.2 |
| `kpp3d` (KPP), cada 3 h | KPPviscA, KPPdiffT | diagnóstico del cierre (como R1k) | 96 |
| `ggl90` (F1g), cada 3 h | GGL90Kr, GGL90ArU | lo mismo para GGL90 | 96 |
| `seaice` (F2, F4) | SIarea, SIheff, SIuice, SIvice | hielo | 4.8 |

Más:
- estado inicial y final de cada etapa (`dumpInitAndLast`): ~0.3 GB por instante;
- pickups permanentes cada 60 h y rotativos: 0.29 GB cada uno;
- copias de la entrada: 0.51 GB × 5, porque el viento son 2 × 217 MB.

**Por corrida (GB):**

| Corrida | 280 h | 240 h |
|---|---|---|
| KPP a 200 m (F1, F3, F5, V3) | 143 | 123 |
| con hielo (F2, F4) | 147 | 127 |
| F1c (sin `kpp*`) | 133 | 115 |
| F1g | 142 | 122 |
| V4 (400 m) | 36 | 31 |
| F2 o F3 a 100 m, caja completa | 571 | 493 |
| F2 o F3 a 100 m, subcaja | 264 | 228 |
| R2 / R3 (Task 3: ~720 GB las dos a 280 h) | ~360 cada una | ~310 cada una |

**Total de la campaña** (9 corridas del fiordo a 200 m, V4 a 400 m, más R2/R3):
- **280 h: ≈ 1.17 TB del fiordo + 0.72 TB del puente = 1.9 TB;**
- 240 h: 1.01 + 0.62 = 1.6 TB.

LaCie tiene ~2.9 TB libres hoy, después de las corridas de esta Task.
- Con F2 y F3 a 100 m en la caja completa, 280 h: +0.86 TB → 2.75 TB. No cabe con margen.
- Con la subcaja: +0.24 TB → 2.1 TB.

Si hace falta, se puede reducir sin tocar el preregistro:
- escribir `kpp3d` cada 6 h (−4.5 GB por corrida);
- quitar los volcados de estado de `dumpInitAndLast`.

---

## 5. Malla de invierno a 100 m (PD0 decisión 1; ruling R13)

R13 pide cronometrar una corrida corta a 100 m con la configuración de invierno F2, np = 5. Lo que se pudo medir:

1. **Configuración F2: imposible por ahora.** F2 no arranca con fronteras de Orlanski en ninguna malla (`pd1_seaice.md` §7). El costo de MITgcm por paso no depende del hielo: `pkg/seaice` sin dinámica ni termodinámica solo reescribe `fu, fv`. Así que cronometré **F3 (invierno, agua abierta)** en la misma malla (`F3_100s`, `F3_100`).
2. **Caja completa a 100 m** (`fjord_vmf_100`: 1000 × 600 × 40; 5 teselas de 200 × 600):
   - **Batimetría:** el procedimiento de la Task 1b no corre a 100 m. `bathy.to_model_grid` falla porque 31 729 celdas no contienen ningún punto de IBCAO a 100 m. Para cronometrar uso la malla de 200 m partida 2 × 2. Una corrida de producción necesitaría otra interpolación (bilineal o la de la Task 1b a 50 m y luego promedio), que no existe.
   - **Memoria:** RSS a 200 m de 0.90 GB por rango. Escalada con el tamaño de tesela con halo ((200 + 8)(600 + 8) frente a (100 + 8)(300 + 8), ×3.80), da ≈ 3.4 GB por rango, **≈ 17 GB en total, frente a 18 GB de RAM**. Datos estáticos: **más de 4 GB por rango; el ejecutable no enlaza.** `ld: ARM64 ADRP out of range (4308848640 max is +/-4GB)`, en `/Volumes/ESTEBAN_MIT/portfolio_builds/exp4/fjord_vmf_100`. Los arreglos estáticos de MITgcm no caben en el modelo de código por omisión de arm64.
   - **No se puede compilar ni correr** en esta máquina sin reestructurar la memoria. Sería otro modelo de código, menos paquetes o más rangos, y lo último lo prohíbe la regla de ≤ 5.
   - Aunque cupiera: escalando el s/paso medido en la subcaja (1.32 s con 277 k columnas) a 600 k columnas da ≈ 2.9 s/paso. Una corrida de 280 h tardaría ≈ 13 h con Δt 60 s y ≈ 8 h con Δt 120 s; pasa o roza el límite de 8 h por corrida, y F2 + F3 se comerían el total de 30 h.
3. **Subcaja anidada a 100 m** (`fjord_vmf_100s`: x 482–555 km, y 8616–8654 km; 730 × 380 × 40; teselas de 146 × 380). Cubre el fiordo interior, los dos estrechos y el Bellsund interior. Sus bordes O y S cruzan Bellsund (fronteras abiertas); el N es tierra. La batimetría es un recorte de la de 100 m anterior.
   - **El binario por omisión no arranca:** `dyld ... map cache into shared region failed` en los 5 rangos (`F3_100s_t1h_defaultbase_execfail`). Sus datos estáticos (`__DATA`) son 2.93 GB, cargados en 0x1_0026_0000, y se solapan con la región compartida de dyld (0x1_8000_0000). Es el mismo fallo de la Task 3 con 2.12 GB.
   - Re-enlacé solo el ejecutable en el directorio del build, sin tocar la fuente ni el Makefile, con `make LIBS="-L/opt/homebrew/lib -lnetcdff -lnetcdf -Wl,-image_base,<base>"`.
   - Con base 0x10_0000_0000 (64 GB) y 0x8_0000_0000, el proceso muere al arrancar (`Killed: 9`; `F3_100s_t1h_base64GB_killed`). Con 0x3_0000_0000 y 0x4_0000_0000 arranca. Quedó 0x4_0000_0000 (16 GB): `__DATA` en 0x4_0026_0000. El ejecutable por omisión quedó guardado como `mitgcmuv.default_base`.
   - Con base 0x4_0000_0000 (16 GB) arranca, pero el desplazamiento aleatorio de ASLR hace que **~25 % de los arranques de un proceso fallen** (3 de 12 en una prueba en serie). `-no_pie` se ignora en arm64. Con 5 rangos, la probabilidad de que arranquen todos es ~0.75⁵ ≈ 0.24: un intento falló y el siguiente funcionó (`F3_100s_t1h_execfail_try1`, `F3_100s_t1h`).
   - **Resultado:** 1 h a Δt 60 s escalonado en 87 s de pared y 1.32 s por paso (bucle), CFL 0.039.
   - Extrapolado a 280 h: 6.2 h por corrida con Δt 60 s y ≈ 3.8 h con Δt 120 s (suponiendo el mismo +24 % por paso que a 200 m). A 100 m en invierno, el número de ondas internas con 120 s es 0.34 (§3).
   - Si el PD1 elige esta opción, el ajuste del enlazador tiene que pasar al build del repositorio (por ejemplo, un `LIBS`/`LDFLAGS` extra en `scripts/build_mitgcm.sh` por configuración) y declararse.
   - **Física:** la subcaja necesita valores de frontera abierta de una corrida madre a 200 m (anidamiento en un sentido, `OBCS_PRESCRIBE` con archivos de la madre), que no existe todavía. Con Orlanski sola en Bellsund cambia el dominio de F2 y F3 frente a F1 (H4c compara F1 con F3). Disco: 264 GB por corrida de 280 h (§4).
4. **Limitación declarada a 200 m** (lo que dice el PD0 si no cabe): R₁/Δx de invierno = 1.7–6.6 (PD0 §2.3). Costo: nada extra.

**Resumen para la decisión** (punto PD1-W; decidido: (c), §0):

| Opción | Costo extra frente a 200 m (F2 + F3, Δt 120 s) | Disco extra | Qué falta |
|---|---|---|---|
| (a) 100 m en la caja completa | **imposible en esta máquina**: no enlaza (> 4 GB estáticos) y pide ~17 GB de RAM | +0.86 TB | otra máquina o reestructurar el build |
| (b) subcaja anidada a 100 m | +4.7 h (280 h) / +4.1 h (240 h): F2 y F3 a ~3.8 h / 3.3 h cada una, más sus madres a 200 m, que ya están en la campaña | +0.24 TB | anidamiento (archivos de frontera de la madre, `useOBCSprescribe`), ajuste del enlazador en el build y reintentos de arranque; además, F2 no arranca por PD1-S2 |
| (c) **[ELEGIDA]** limitación declarada a 200 m | 0 | 0 | nada; se declara R₁/Δx = 1.7–6.6 en invierno (PD0 §2.3) y la advertencia de H4c |

- Con la opción (b), la campaña con Δt 120 s pasa de 23.0 a 27.7 h (280 h) o de 19.7 a 23.8 h (240 h). Con ~22 h disponibles, **no cabe sin recortar** (F4, F5, V4 y F1g suman 5.1 h a 280 h).
- **Recomendación técnica (la decisión es de Esteban):** (c), la limitación declarada. (b) añade un anidamiento que no existe, con su propia validación, y un build frágil.

---

## 6. Agenda del PD1 (para Esteban)

Cada punto lleva su sección. **Las decisiones de Esteban están en §0.**

1. **Caja** (R5): 100 × 60 km (§1).
   - Alternativa de la Task 1b: 100 × 45 km cerrando Van Keulenfjorden y Recherchefjorden en su boca (−25 % de celdas). Si se cambia, se vuelve a correr la Task 1b.
2. **Malla** horizontal: 200 m (500 × 300); vertical: 40 niveles (§1).
3. **Δt y esquema de tiempo** (§2): punto PD1-Δt. La configuración del PD0 tal cual (paso sincrónico) es inestable con cualquier Δt probado.
4. **viscAh:** Smagorinsky 2.2 + `viscC4smag` 1, sin viscAh constante (PD0 §6). La inestabilidad de §2 aparece igual con el cierre de la tesis y desaparece con el paso escalonado, así que no la atribuyo a la viscosidad (§2).
5. **T_fin:** 280 h o 240 h (§3 y §4: tiempo de pared y disco de las dos opciones).
6. **Resolución de invierno** (§5): 100 m en la caja completa, subcaja anidada a 100 m o limitación declarada a 200 m.
7. **Hielo, PD1-S1** (`pd1_seaice.md` §2, §5 y §6): con las banderas del spec, el océano no recibe arrastre bajo el hielo, ni fijo ni a la deriva. Está confirmado en corrida corta: τ = 0.0 bajo A = 1.
   - Para F4 hay una opción solo de namelist: `useHB87StressCoupling = .TRUE.`. El océano bajo el hielo recibe el viento sobre el hielo, sin el término de Coriolis del hielo (~0.03 N m⁻² a 0.2 m s⁻¹).
   - No sirve para F2: dejaría pasar todo el viento a través del hielo fijo.
   - La deriva libre de MITgcm no tiene inercia del hielo (spec: "losa con inercia"); clasificado (b).
   **Hielo, PD1-S2** (`pd1_seaice.md` §7): F2 y F4 **no arrancan** con fronteras de Orlanski (`pkg/obcs/obcs_check.F:400-405, 571-579`). Sin resolver esto no hay H4b.
8. **F4** (§1, punto PD1-F4): polígono del hielo a la deriva.
9. **V3** (§1, punto PD1-V3): forzamiento M2 mínimo (U₀ = 1.70 cm s⁻¹ en las tres fronteras, η₀ = 0.5 m). El error de las entradas (la marea era cero) está corregido. Con la marea activa, V3 diverge con Δt 120 y 60 s por las corrientes de los estrechos: hace falta un Δt propio ≲ 20–30 s o un forzamiento más suave.
10. **Pendientes de batimetría** (`notes/pd1_batimetria.md`):
    - (a) la caja de 100 × 60 km (punto 1);
    - (b) el control del ancho de Mariasundet **FALLA** (1.33 km de celdas mojadas frente a 0.7 km; clasificado (a), definición de la métrica, *post hoc*). A 400 m (V4) los dos estrechos siguen abiertos: 3 celdas en Akselsundet y 6 en Mariasundet sobre la sección;
    - (c) el corte de R7 en y = 8598 km deja un borde recto artificial en la cabeza de Van Keulenfjorden;
    - (d) Vallunden cerrada (R6);
    - (e) el 51 % de las celdas del fiordo interior no tiene dato medido en IBCAO: 50.0 % interpolado (TID 40–47) y 1.4 % con código de tierra. Va como limitación en el README (Task 8).
11. **Presupuesto** (§3 y §4): el recorte acordado si hace falta es F4, F5, V4 y F1g, en ese orden; R1k y F1c no se recortan.
