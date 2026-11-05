# Bitácora

## 2026-11-05 — Task 1: repositorio, paquete común y runner

- Se montó `pfcommon` (E/S binaria big-endian y MDS, estilo de figuras, estado de corridas), `scripts/build_mitgcm.sh` y `scripts/run_mitgcm.sh`, con pruebas (6 PASS).
- El runner copia (nunca enlaza) a LaCie, verifica el montaje (código 2 si falta), fija `ulimit -n 6000` y hilos a 1, escribe `run_status.json` y revisa `pgrep` al final.
- Hallazgo: en corridas seriales MITgcm no escribe `STDOUT.0000` (la bitácora del modelo sale por stdout), así que el runner usa `output.txt` como respaldo.
- Prueba de humo (`tutorial_plume_on_slope`, 20 pasos, np=1): compilación ~60 s en `/Volumes/ESTEBAN_MIT/portfolio_builds/smoke/`; corrida ~2 s en LaCie; `run_status.json` = OK, `last_iter` 20.
- Decisión: no se hizo `git init` (el repo ya existía); solo se amplió `.gitignore`.

## 2026-11-05 — Task 2: teoría de plumas y derretimiento de 3 ecuaciones

- Se implementaron `exp1_fjord_plumes/theory/melt.py` (3 ecuaciones; modos `mitgcm`, `jenkins`, `constant`; bandera `melt_off`, fallo R2) y `theory/plume.py` (pluma lineal de Jenkins 2011, semicono de Slater et al. 2016 y cono completo de MTT). Pruebas: 23 PASS en `exp1_fjord_plumes/tests`.
- **Hallazgo: MITgcm `pkg/icefront` NO usa transferencia constante.** `icefront_thermodynamics.F:91-97`: γT = max(1e-2·√(1.5e-3)·|wVEL|, 1e-4) m/s y γS = 5.05e-3·γT, todo fijo en el código. El piso gana para |w| < 0.258 m/s. Estos valores NO son los números de Stanton de Jenkins (2010, 2011) (1.1e-3, 3.1e-5) que reporta Xu et al. (2012, Tabla 1). El modo `mitgcm` reproduce el Fortran línea por línea (hay una prueba contra la transcripción).
- Signo de `ICFfwFlx` (kg m⁻² s⁻¹): negativo = derretimiento (`icefront_thermodynamics.F:50,141-146`). m_hielo = −ICFfwFlx/ρ_i; m_agua = −ICFfwFlx/rhoConst.
- Constantes y su origen:
  - Modo `mitgcm` (defaults de `MeltParams`): λ1 = −0.0575, λ2 = 0.0901, λ3 = 7.61e-4 (`icefront_thermodynamics.F:77-79`); L = 3.34e5 (`icefront_readparms.F:84`); c_i = **2000** (`icefront_readparms.F:85`; el brief decía 2009, que es el valor de Jenkins 2011); ρ_i = 917; T_hielo lineal desde −20 °C en superficie hasta 0 °C en la base del frente (`icefront_thermodynamics.F:121-123`, `icefront_readparms.F:88`). c_w = 3974, como en el brief y en `isomip/input/data:19`; ojo, el default de MITgcm es 3994 (`set_defaults.F:174`).
  - Modo `jenkins` (Jenkins 2011, Tabla 1): Cd = 2.5e-3, √Cd·ΓT = 1.1e-3, √Cd·ΓS = 3.1e-5. `MeltParams.jenkins2011()` da además λ1 = −5.73e-2, λ2 = 8.32e-2, L = 3.35e5 y c_i = 2009.
  - EOS lineal: tAlpha = 3.87e-5 y sBeta = 7.86e-4 están verificados en Jenkins 2011, Tabla 1. **No vienen de `verification/isomip`**, que usa `eosType='JMD95Z'` y `rhoConst=1030`. ρ0 = 1027 sigue `[verificar]` (no hay fuente local). Las corridas deben fijar `eosType='LINEAR'`, `tAlpha`, `sBeta` y `rhoNil = rhoConst = 1027`.
- **Cambio de constante respecto al brief (prueba de altura máxima).** MTT 1956 da x = 0.410 α^(−1/2) F0^(1/4) G^(−3/8) x1 (Ec. 14), con u = 0 en x1 ≈ 2.8 (p. 9) y α = 0.093 (perfiles gaussianos, p. 12). Con F0 = flujo total de flotabilidad g′0·Q0, z_max N^(3/4)/F0^(1/4) = 0.410·2.8·0.093^(−1/2) = **3.76**, no 5.0. El 5.0 = 3.76·π^(1/4) corresponde a F0 definido por unidad de π `[verificar]` (Turner 1973, no está en `literatura/`). Además, nuestras ecuaciones son de sombrero de copa: en forma de flujos son idénticas a las gaussianas de MTT si α_sombrero = √2·α_gauss. Por eso la prueba usa α = √2·0.093 y espera 3.76 ± 0.3; el modelo da 3.767. Comprobación propia del sistema adimensional de MTT: u = 0 en x1 = 2.805 y flotabilidad nula en x1 = 2.131 (el artículo da 2.8 y 2.125). Con el α = 0.093 literal del brief en ecuaciones de sombrero de copa saldría 4.48: no pasa la tolerancia del brief y no es lo que dice el artículo.
- Elección numérica: para la ley u ∝ z^(−1/3) el ajuste log-log se hace respecto al origen virtual (estimado del ajuste lineal de b), no de la fuente; con una fuente finita la pendiente medida desde la fuente no es −1/3 a 20 b0. La altura máxima también se mide desde el origen virtual 5b0/(6α) (corrección del 0.3 %).
- Elección numérica: el estado es (Q, M², QT, QS). Así u → 0 queda como un cruce regular de cero (evento terminal), sin la singularidad 1/M de las ecuaciones del cono.
- Diferencia de modelo (b, física faltante, a considerar en Task 5): la teoría usa T tal cual. MITgcm convierte θ a temperatura in situ (`SW_TEMP`, línea 114) y evalúa wVEL en la cara superior de la celda.

## 2026-11-05 — Task 3: generadores de entrada (gendata)

- Se implementaron `exp1_fjord_plumes/gendata/{profiles,grids,sources,make_case}.py` y `cases.yaml` (24 casos). 15 pruebas nuevas PASS (52 en total).
- **Valores de los perfiles cambiados respecto al brief.** Sciascia et al. (2013), Tabla 2 (estudio de invierno en Sermilik): agua polar −1.5 °C / 32.9 psu, agua atlántica 4 °C / 34.6 psu (el brief decía 33.0, 3.5 y 34.8). Se usan los de Sciascia. Straneo & Cenedese 2015 y Jackson 2014 no tabulan valores por capa. Jackson 2018 centra la picnoclina de Sermilik en 160 m (Δρ = 1.1 kg m⁻³); se mantienen 200 m y 75 m del brief como `[verificar]`.
- `icefrontlength` = 1/dx en la columna del frente (`icefront_description.tex:15-16`; es el interruptor de `icefront_thermodynamics.F:89`). `R_icefront` se lee en `icefront_init_fixed.F:73-75`.
- Decisión: la celda adyacente a la cara es la propia columna del frente (`i = front_i`), porque `icefront_thermodynamics.F` aplica el flujo a la celda mojada con `icefrontlength > 0`; el brief decía `front_i + 1`. Queda como parámetro `face_offset`.
- Ronda de correcciones Task 3: (c) escribe solo la fuente de fondo ρ0·Q0 (marcador `bottom_source_placeholder`; Task 7 pondrá sumideros de entrada y fuentes de salida con neto Q0). hFacC para las revisiones de descarga sale de la batimetría generada (opción `land_columns`). En (b) el perfil continuo de Sermilik se centra en 100 m (H − H2), con la misma forma tanh, ancho y valores de capa, para que la capa inferior mida 200 m como en el caso de dos capas; (a) y (c) mantienen 200 m. `make_case` valida los archivos nombrados antes de escribir nada.

## 2026-11-05 — Task 4 (Steps 1–4): configuración (a) `a_line2d`, corrida corta, PD1

- Configuración en `exp1_fjord_plumes/config/a_line2d/` (code: `SIZE.h` 193×1×50 serial, `packages.conf` gfd/obcs/icefront/diagnostics, `CPP_OPTIONS.h` con `ALLOW_ADDFLUID` y `ALLOW_NONHYDROSTATIC`, `ICEFRONT_OPTIONS.h`, `OBCS_OPTIONS.h` (solo este, Orlanski, esponja; archivo extra no listado en el brief porque `ALLOW_OBCS_SPONGE` viene apagado), `DIAGNOSTICS_SIZE.h`; plantillas `data*`, `eedata`). Tabla completa con fuentes en `exp1_fjord_plumes/notes/pd1_a_line2d.md`.
- **Viscosidad y difusividad: Sciascia et al. (2013), Tabla 1** (2D NH MITgcm, 10 m × 10 m en el frente): νh = κh = 0.25 m² s⁻¹, νv = 1e-3, κv = 2e-5. Xu et al. (2012) usan νv = 0.1, biarmónica 300 m⁴ s⁻¹ y κ = ν (Pr = 1); queda como alternativa de sensibilidad. Pared y fondo sin deslizamiento, como Sciascia 2013 ¶14.
- EOS lineal con los valores de la teoría (Task 2): tAlpha 3.87e-5, sBeta 7.86e-4, rhoNil = rhoConst = 1027, Cp 3974, g 9.81. icefront con la ley de transferencia por defecto y T del hielo lineal (−20 °C en la superficie, 0 °C en la base); constantes escritas explícitamente en `data.icefront`.
- Requisitos de `selectAddFluid=1`, verificados en la fuente: `ALLOW_ADDFLUID` (`config_check.F:363-373`), sin tapa rígida (`:837-845`), `staggerTimeStep` (aviso `:846-855`). `convertFW2Salt` queda en −1 automáticamente (`ini_parms.F:648-652`), así que la T/S del agua añadida es `temp_addMass/salt_addMass` también con superficie libre lineal (`apply_forcing.F:505-507`). **No** hace falta superficie libre no lineal ni `exactConserv`.
- **Hallazgo (cierre del dominio):** MITgcm es periódico en x. Con la cara en la columna 0 mojada y OB solo al este, el borde oeste de la columna 0 conectaba con la OB y `OBCS_CHECK` habría pedido una OB oeste (`obcs_check.F:645-669`). Solución: la última columna es tierra (`land_columns: [[0, -1]]` en `cases.yaml`) y `OB_Ieast = −2`. Así `front_i = 0`, dx₀ = 5 m y la fuente (fallo R8) no cambian.
- **Hallazgo (OBCS balance):** `useOBCSbalance` solo compensa EmPmR, no addMass (`obcs_balance_flow.F:300-315`), así que se dejó apagado. Comprobado: η medio = 3.5 mm a 2 h (q0 = 0.03), contra 10.8 mm si el volumen quedara atrapado; el agua sale por Orlanski.
- **Hallazgo (salidas):** las instantáneas de diagnóstico se escriben un paso tarde (`diagnostics_write.F:62-72`), así que no hay instantánea en t = endTime. El estado final está en `T/S/U/V/W/Eta.<nEnd>.data` (`dumpInitAndLast`) y en el pickup (`pChkptFreq = 7200`). El Step 5 (C3) debe comparar `T` y `Eta` de esos volcados o del pickup.
- Interfaz: `make_case(..., overrides=None)` (nuevo kwarg mínimo, p. ej. `{"endTime": 7200}` para la corrida corta). `cases.yaml` (a): `endTime` 172800 (48 h), `deltaT` 3, `land_columns`.
- `scripts/run_campaign.py` (glob de casos, compila una vez y luego corre en serie, se niega a arrancar si `pgrep` encuentra mpirun/mitgcmuv, verifica `SIZE.h` contra la malla, np ≤ 5). 12 pruebas (RED 12 fallos con un stub → GREEN 12; suite completa 66 PASS).
- Compilación serial en `/Volumes/ESTEBAN_MIT/portfolio_builds/exp1/a_line2d/`: ~65 s.
- **Intentos de la corrida corta** (`a_line2d_twolayer_q0.03`, 2 h, np = 1; los intentos previos se conservan como `..._short2h_tryN` en LaCie):
  1. dt = 2 s, `cg3dMaxIters` = 40: OK, 34 s, CFL máx. 0.134. **cg3d nunca converge** (siempre 40 iteraciones, residuo ~4e-5 frente a 1e-8).
  2. dt = 2 s, cg3d 300: OK, 88 s; converge en 280–300 iteraciones. Diferencia contra el intento 1 a 2 h: máx. 0.011 °C en T, 2 mm/s en w.
  3. dt = 4 s, cg3d 400: OK, 52 s, CFL 0.27, 322–390 iteraciones. Frente a dt = 2 s: ICFfwFlx medio igual a 0.04 %; T rms 0.008 °C (máx. puntual 0.37 °C en la interfaz, desfase de remolinos).
  4. Igual que el 3, con cg3d 500: OK, 52 s (configuración idéntica, salvo el techo de iteraciones).
  5. q0 = 0.1, dt = 4 s: OK, CFL w **0.42**, poco margen para 48 h.
  6. **Final: dt = 3 s para los 6 casos** (no se ajusta por caso). q0.03: OK, **65 s**, CFL u/w 0.20/0.20. q0.1: OK, 67 s, CFL 0.25/0.32. Sin NaN; `pgrep` vacío tras cada corrida.
- Costo extrapolado: 24 × 65–67 s ≈ 26–27 min por caso de 48 h, **≈ 2.6 h los 6 casos**. Dentro del presupuesto (≤ 8 h por corrida, ≤ 24 h el Exp. 1).
- Figura rápida `notes/pd1_a_quicklook.png`: la pluma sube pegada a la cara y se despega en la interfaz (~200 m). ICFfwFlx medio en la cara = −5.6e-3 kg m⁻² s⁻¹ (≈ 0.53 m d⁻¹ de hielo) en la segunda hora.
- **PD1 pendiente**: no se lanzó la campaña (fallo R4). Steps 5–7 quedan para después del visto bueno.

## 2026-11-05 — Task 4 Step 5: reinicio bit a bit (C3) en (a)

- PD1 aprobado por Esteban: cg3d convergido (500 iteraciones, 1e-8), dt = 3 s, viscosidad y difusividad de Sciascia 2013, pared sin deslizamiento (diferencia de clase (c) frente al Cd·u² de la teoría).
- Caso `a_line2d_twolayer_q0.03`, np = 1, en `/Volumes/LaCie/portfolio_runs/exp1/`: `c3_restart_4h` (4 h seguidas, 4800 pasos, OK, 120 s) contra `c3_restart_2h_a` (2 h, OK, 66 s) + `c3_restart_2h_b` (desde `pickup.0000002400` + `pickup_orlanskiE.0000002400.001.001`, `nIter0 = 2400`, 2 h más, OK, `last_iter` 4800, 58 s). Las entradas salen de `make_case(..., overrides={"endTime": ...})`; en el tramo b solo se cambió `nIter0` en `data`.
- Comparación de los volcados de estado finales (iteración 4800; la instantánea de diagnóstico en endTime no existe): `np.array_equal` es **True** para `T` y `Eta` (y también para S, U y W; máx. |diff| = 0). Sin avisos en el log del reinicio. **C3 cumplido.**
- Se borraron de LaCie los directorios de prueba `*_short2h_try*` y los de q0.1; queda solo `a_line2d_twolayer_q0.03_short2h` (lo cita PD1).

## 2026-11-05 — Task 5 Steps 1-2: módulo `plume_diag` (aún sin campaña)

- `exp1_fjord_plumes/analysis/plume_diag.py`: `melt_rate_m_per_s` (−ICFfwFlx/ρ_i, positivo = derretimiento; 1.0905e-7 m/s = 3.44 m/año), `detach_depth`, `plume_profile` (se niega con `require_ok`; usa snapshots de diagnóstico más los volcados T/S/U/W finales porque no hay snapshot en endTime). Ancho de pluma = de la cara al borde exterior de la racha contigua con w ≥ 10 % de w_max(z); `z_detach` = primer nivel sobre el máximo de w de la columna de la cara donde w < 5 % de su máximo (interpolado, niveles nativos de WVEL).
- `fig1_plume_validation.py` escrito (C2/C2'/C2'', `known_differences` en `fig1_metrics.json`); probado solo con la corrida corta de 2 h (ventana 1-2 h, salida en scratch): z_detach −186.6 m (teoría −190.6), derretimiento medio bajo el despegue 0.84 m/d (razón 1.07 frente a la teoría en modo mitgcm). Corrida de 2 h no representativa del estado casi estacionario; Steps 3-5 esperan la campaña.

## 2026-11-05 — Task 6 Step 1: configuración (b) `b_twoplume3d` (sin correr)

- `config/b_twoplume3d/{code,input_templates}` a partir de (a): misma física aprobada en PD1 (EOS lineal, icefront por defecto, addMass con selectAddFluid=1, cg3d 500/1e-8, Sciascia 2013, sin deslizamiento, MDS float64 big-endian, useSingleCpuIO, pickups). Cambios: malla 100 × 120 × 30 (10 m isotrópico), VVEL en los diagnósticos, snapshots cada 30 min, icefrontDiag promediado cada 15 min, pChkptFreq 5400 s.
- Fronteras: paredes laterales = filas de tierra j = 0 y Ny-1 (cierran también el período en y; ancho mojado 1180 m); última columna de tierra y OB_Ieast = -2 (Orlanski + esponja), como en (a). Nuevas claves `land_j_rows` / `land_i_cols` en `make_case` (5 pruebas nuevas; RED 2 fallos → GREEN; suite completa 75 PASS).
- `cases.yaml` (b): `endTime` 6 h → 3 h (plan, Step 2); dt = 1 s provisional.
- x0 = 0: las dos fuentes caen en la misma celda (una fuente de 2·Q0); `_single`: una fuente de Q0 en j = 60. Con x0 = 320 m las fuentes quedan a ≥ 425 m de las paredes.
- SIZE.h: sNx 100, sNy 24, nPy 5 (5 procesos en y). Compilación MPI OK en `/Volumes/ESTEBAN_MIT/portfolio_builds/exp1/b_twoplume3d/` (74 s). No se corrió MITgcm (campaña (a) en curso).
- Costo a priori (escalado de (a), 2.8e-6 s por celda·paso): 3.9e9 celda·pasos por caso ≈ 3.0 h serial, ≈ 0.6–0.9 h con 5 procesos → 14 casos ≈ 8.5–12 h (incertidumbre ×2 por cg3d en 3D). La corrida corta de Step 2 (30 min) debería tardar ≈ 6–9 min.

## 2026-11-05 — Task 7 Steps 1-3: PD3 (addMass como sumidero) y piezas en Python de (c)

- **PD3, spike** `config/spike_addmass` (10 × 10 × 10, hidrostática, np = 1, 300 pasos, 2 s por corrida; corridas en `/Volumes/LaCie/portfolio_runs/exp1/spike_addmass/`). Fuente: `apply_forcing.F:504-530/873-898` aplica `addMass·(temp_addMass − θ)` sin importar el signo, y la advección es en forma advectiva (`gad_calc_rhs.F:764-779`). **Con `temp_addMass` definido, `addMass < 0` retira agua a `temp_addMass`, no local**: el sumidero bajó de 2 a 1.94 °C y subió de 34 a 34.03 psu. Si `temp_addMass` queda sin definir y `salt_EvPrRn = UNSET_RL` (por defecto es 0, `set_defaults.F:265`), el sumidero queda exactamente en 2 °C / 34 psu, y la fuente con RBCS (máscara Q·τ/V) coincide con la de `temp_addMass` (7e-5 °C). El volumen se conserva (flujo del sumidero −2.000000 m³/s). Detalle y alternativas (i)/(ii)/(iii) en `notes/pd3_addmass_sink.md`. **Recomendación (iii); espera decisión de Esteban.**
- `gendata/plume_source_sink.py`: sumideros = −ΔQ en celdas enteras bajo la celda de z_neutral; toda la salida Q(z_lo) en esa celda (Cowton et al. 2015, §2.1.3), con la T y S de la pluma en su cara inferior. Suma = Q0 a 1e-13. No se reparte entre z_neutral y z_detach porque, en Sermilik con Q0 = 100, z_neutral = −181 m y z_detach ≈ −3 m (sobrepaso de 18 celdas). Q_out = 3725 m³/s, T/S de salida 3.42 °C / 33.53 psu. El derretimiento es ~1e-4 del arrastre.
- Hallazgo en la teoría (Task 2): con transferencia `mitgcm`/`constant` el semicono falla justo bajo z_detach, porque el derretimiento no se anula cuando u → 0 y b → ∞. `plume_source_sink` reintenta bajando el techo de integración; solo necesita z ≤ z_neutral.
- `analysis/exchange.py`: `exchange_flux` = ½∫|u|dA. Jackson et al. 2018 no la escriben así en §2; usan Q_in = −Q_out (§3c(3), Ec. 26). Con transporte neto nulo son iguales. `heat_toward_front` = ρ0 c_p ∫(−u)(T − T_ref)dA, con T_ref explícita.
- TDD: RED (módulos inexistentes) → GREEN; suite completa 82 PASS.

## 2026-11-05 — Task 2, seguimiento R15: término de la semicónica con derretimiento

- Falla: `point_plume(half_cone)` con derretimiento en modo `mitgcm`/`constant` se rompía justo bajo z_detach. El ancho de contacto 2b ∝ M^(−1/2) diverge cuando u → 0 y el flujo de derretimiento no se anula (piso γT = 1e-4 m/s), así que las fuentes 2b·m_w explotan.
- Solución (elección física y numérica): un evento terminal en u = u_stop = 1e-3 m/s, solo para la semicónica con derretimiento. Por debajo de ~1 mm/s la pluma ya no se distingue de las corrientes del fiordo y un radio mayor que la profundidad no tiene sentido. z_detach baja ~u_stop²/(2|g′|) ≈ mm; con u_stop = 1e-4 cambia 0.3 mm.
- El resto de los resultados queda bit a bit igual. El reintento de `plume_source_sink` ya no hace falta (se deja, es de Task 7); su q cambia ≤ 1.6e-6 relativo solo por la malla fina distinta.

## 2026-11-05 — Task 7 Step 4: configuración (c) `c_fjord3d` (compilada, sin correr)

- PD3 aprobado por Esteban, opción (iii). Con la corrección R15 de la teoría (1ae1414) se quitó el reintento de `plume_source_sink`. Nuevo `n_out_levels`: la salida se reparte entre la celda de z_neutral y las de encima.
- `make_case` (c) ya no usa el marcador R10, ahora es `plume_representation = "bpt_sink_source_rbcs"`:
  - addMass en la columna de la cara, repartido en 3 celdas a lo ancho (j = 14–16). Salida en k = 17–18 (z_neutral = −181 m, 3725 m³/s, 3.42 °C / 33.53 psu); arrastre de 3625 m³/s en k ≥ 19.
  - Neto = ρ0·Q0 exacto (C5).
  - RBCS con máscara q·τ/V = 0.93 (τ = 600 s) solo en las 6 celdas de salida.
  - `temp_addMass`/`salt_addMass` sin definir y `salt_EvPrRn = UNSET_RL`.
- Repartir la salida evita el recambio de 107 s de una sola celda. Con el reparto queda en 644 s, λ·dt = 0.09 con dt = 60 s.
- Física:
  - hidrostático; f = 1.37e-4;
  - viscosidad y difusividad de Cowton 2015, Tabla 1 (30 / 1e-5 m²/s, misma malla de 200 m × 10 m);
  - `ivdc_kappa = 10`;
  - free-slip con arrastre de fondo 2.5e-3;
  - EOS lineal como (a)/(b), MDS float64 big-endian.
  - Sin icefront: ≈ 10–15 m³/s de derretimiento, ≈ 0.3 % de la salida de la pluma, y ensuciaría los casos sin pluma.
- Forzamientos:
  - viento down-fjord de 0.2 N/m² durante 12 h, con el protocolo de 2 etapas de la tesis (`data` 0–24 h con `windx.bin`; `data.stage2` desde la iteración 1440 sin viento);
  - plataforma: OB con T/S prescritas (picnoclina ±30 m, 6 d, registros horarios) y u/v por Orlanski (`OBCS_PRESCRIBE_READ` va después de `ORLANSKI_EAST`), más la esponja de T/S de 2 km;
  - paredes de tierra en j = 0 y Ny−1, columna Nx−1 de tierra, `OB_Ieast = -2`, sin `useOBCSbalance`.
- Malla 200 × 30 × 50 = `SIZE.h` (nPx = 5, sNx = 40). Build MPI OK en `/Volumes/ESTEBAN_MIT/portfolio_builds/exp1/c_fjord3d/`.
- W/Rd = 0.67 (Rd = 8.4 km; con ny = 37 sería 0.83, no se cambió).
- Costo a priori (2.8e-6 s por celda·paso, de (a), no hidrostático): 10.1 h serial por caso, 2–3 h con 5 procesos, 8–12 h los 4 casos. Es una cota superior porque no hay cg3d. La corrida corta de 1 día debería tardar ≈ 4–6 min. Nota: `notes/pd1_c_fjord3d.md`.
- TDD: 9 pruebas nuevas, RED → GREEN; suite completa 92 PASS. No se corrió MITgcm.
