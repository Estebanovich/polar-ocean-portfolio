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

- PD3 aprobado por Esteban, opción (iii). Con la corrección R15 de la teoría (ae36d51) se quitó el reintento de `plume_source_sink`. Nuevo `n_out_levels`: la salida se reparte entre la celda de z_neutral y las de encima.
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

## 2026-11-05 — Task 7: correcciones previas a la corrida (R18, R19)

- R18: la salida de la pluma ahora se reparte solo a lo ancho: 6 columnas (j = 12–17), 1 nivel (k = 18). La mitad que antes iba en k = 17 era 0.068 kg/m³ más densa que el ambiente (inversión); en k = 18 es 0.049 más ligera. El recambio sigue en 644 s, u de salida ≈ 0.31 m/s y la máscara RBCS en 0.93.
- R19: `notes/pd1_c_fjord3d.md` documenta dos diferencias de clase (b): el forzamiento de la pluma es fijo (no responde a la picnoclina y no conserva estrictamente calor y sal), y ΣQ = Q0 incluye 0.39 m³/s de deshielo retirado del ambiente (también en `case.json`). `plume_source_sink` levanta `ValueError` si `n_out_levels` pasaría de la superficie. Se corrigieron los encabezados de `OBCS_OPTIONS.h` (c) y `SIZE.h` (spike).
- Fig. 3: la esponja cubre i = 188–197, así que la sección de Q_ex debe ir a media distancia del fiordo.
- Menores: el import de las pruebas pasa arriba y `_c_forcing` devuelve `externForcingCycle`. TDD: 3 pruebas nuevas, RED → GREEN; suite completa 95 PASS. Se recompiló (c) con MPI (cambió `code/`, solo comentarios). No se corrió MITgcm.

## 2026-11-05 — PD2 / fallo R20: preparación de convergencia y sensibilidad de (a) (nada corrido)

- Motivo: C2' falla en la evaluación provisional (Q/b = 0.52–0.56 de la teoría; w del núcleo 0.88–0.95; ancho de 20–24 m contra 7.5 m). Atribución propuesta: resolución más pared sin deslizamiento (clase c), posiblemente con entrainment extra (clase b). Se prueba la atribución; la métrica no se redefine. Plan y criterios fijados de antemano en `exp1_fjord_plumes/notes/pd2_convergence.md`.
- Casos nuevos (24 h, twolayer, q0 = 0.03): `conv_dx10` (dx0 = dz = 10 m, dt = 5 s, 170×50), `conv_dx5` (5/5 m, dt = 2.5 s, 193×100), `conv_dx2p5` (2.5/2.5 m, dt = 1.25 s, 216×200), `sens_freeslip` y `sens_xu2012` (malla de producción, dt = 3 s). dt = 0.5 s/m × min(dx0, dz), por CFL.
- Xu et al. 2012 (p. 3): νv = 0.1, biarmónica 300 m⁴/s con dx = 20 m, κ = ν. La biarmónica literal es inestable con dx = 5 m (límite L⁴/(32 dt) = 6.5 m⁴/s, `mom_calc_visc.F:48-49`). Se escala con el mismo tiempo de amortiguamiento de malla: 300·(5/20)⁴ = 1.171875 m⁴/s (viscA4 = diffK4T = diffK4S). Elección numérica documentada.
- `make_case`: la malla estirada se elige por la clave `dx0` y la fuente lineal por `discharge: line` (antes por el nombre de config). Nueva clave `template_config`. Las plantillas `config/a_line2d_study` (solo `data` tiene variables; el resto son enlaces simbólicos a `a_line2d`) reproducen la producción exactamente; lo verifica una prueba. Configs `conv_*` con su propio `SIZE.h` (el resto, enlaces) y `sens_line2d` (todo enlazado; nombre propio para que `--config a_line2d` siga eligiendo solo los 6 casos de producción). En `cases.yaml`, ancla `study_defaults` con los valores de producción.
- Seguridad de la campaña en curso: la entrada regenerada de `a_line2d_twolayer_q0.03` es idéntica byte a byte a la que escribió la campaña (los 16 archivos). `cases.yaml` se escribió de forma atómica.
- Pruebas: `tests/test_pd2_cases.py` (12). RED: 11 fallos + 1 pasa (la guardia de producción). GREEN: 12. Suite: 108 PASS.
- Builds: `/Volumes/ESTEBAN_MIT/portfolio_builds/exp1/{conv_dx10, conv_dx5, conv_dx2p5, sens_line2d}`, 64–66 s cada uno; `SIZE.h` verificado en cada build.
- Costo estimado (c = 2.81e-6 s/celda·paso): 7, 31 y 140 min, y 13 + 13 min; base 3.4 h. Con crecimiento de cg3d, 3.4–10.9 h. `conv_dx2p5` podría pasar de 8 h: cronometrar 1 h de modelo primero y, si hace falta, partir con pickup (C3 verificado).

## 2026-11-05 — PD2 ronda de correcciones 1 (fallos R20, R21): criterios mecánicos (nada corrido)

- `notes/pd2_convergence.md` reescrito con criterios mecánicos:
  - ventana de 12–24 h para todos los estudios y para la referencia de producción;
  - estacionariedad: |ΔR_Q(12–18 vs 18–24)|/R_Q < 0.02;
  - σ = desviación estándar de las 12 muestras horarias (t−1, t);
  - umbrales max(2σ_AB, 0.03) para R_Q y max(2σ_AB, 2.5 m) para b;
  - Richardson con r_eff variable (r21 = 1.881, r32 = 1.912; Celik et al. 2008 `[verificar]`), inconclusa si no es monótona o si p ∉ [0.5, 3];
  - tabla completa de clases de resultado (validez, ancho fijado por la malla, apoyada, parcial, solo tendencia, no apoyada, refutada por descenso o por plano, no monótona);
  - resultados aparte: sensibilidad a dz (producción frente a `conv_dx5`), pared (`sens_freeslip`) y cierre.
- **R21:** `cg3dMaxIters = 2000` en todos los estudios. Una corrida con algún paso en el tope es inválida. Como `cg3d_iters` solo se imprime en tiempos de monitor (`solve_for_pressure.F:435-447`, sin otro aviso por paso en `cg3d.F`), los estudios usan `monitorFreq = deltaT`: se pueden revisar cada paso y su CFL. Solo afecta a la salida; `output.txt` crece a ~0.1–0.4 GB por corrida (estimado).
- Refinamiento no uniforme: el campo lejano queda fijo; r_eff en 0–25 m es 1.91 y 1.88 (1.89/1.78 a 20 m; 1.63/1.45 a 100 m); L difiere < 0.4 %.
- `sens_xu2012` pasa a llamarse "cierre al estilo Xu (A4 reescalada a la malla)" (Xu et al. 2012, Ann. Glaciol. 53(60), p. 231). Con viscAh = 0 el arrastre sin deslizamiento sobre w casi desaparece (`calc_gw.F:456-465`, `mom_w_sidedrag.F:73-80`) y viscAr sube 100×. Por eso ΔR_Q no se atribuye al cierre; se solapa con `sens_freeslip`.
- Costo revisado (t_paso = celdas/9650 × (7.04 + 0.060·N_it) ms, ajustado a las corridas medidas): 3.5 h con N_it = 350; 6.2 h con 700; 15.9 h con 1999. `conv_dx2p5` pasa de 8 h si N_it > 1434. Sonda de 1 h; si da > 8 h, se para y se reporta. Partir en tramos solo con aprobación de Esteban.
- Menores: la prueba de free-slip ahora compara contra producción (solo cambian las dos líneas de pared, más el tope R21 y el monitor por paso); prueba R13 para C2' (pasa con el código actual; verificada por mutación: con C2p en modo jenkins falla); `case.json` registra `git_commit` y `git_dirty`.
- Pruebas: RED 7 fallos → GREEN; `test_pd2_cases.py` + `test_plume_diag.py` 19 PASS; suite 110 PASS. Producción: 15 de 16 archivos de entrada idénticos; `case.json` solo añade los campos git. La campaña en curso no se ve afectada: su proceso ya había importado `make_case`.

## 2026-11-05 — Task 5 Steps 3-5: Figura 1 y métricas C2/C2'/C2'' (campaña (a) completa)

- Seis casos (twolayer/sermilik × q0 0.01/0.03/0.1), 48 h, ventana 36–48 h, solo lectura; veredictos con la teoría en modo `mitgcm` y `ice_depth = H` (R13). Figura `figures/fig1_plume_validation.{png,pdf}` y `figures/fig1_metrics.json` (con `known_differences` por clase, `post_hoc` aparte y referencia a `notes/pd2_convergence.md`).
- **C2**: 5/6 PASA; sermilik q0.1 FALLA (z_detach −110 m contra −181 m de la teoría, 0.141 H). **C2'**: FALLA en los 6 casos con la definición pre-registrada (razón 0.51–0.65); diagnóstico: w de la celda de la cara y w_max están a 0.86–1.05 de la teoría, M/Q a 0.65–0.82; la pluma del modelo es 2–3 veces más ancha que la teórica (b = 7.5 m = 1.5 celdas) con 1.5–1.8 veces el flujo. Métricas post hoc solo etiquetadas como tales (R20). **C2''**: PASA 6/6 (razón 1.00–1.23).
- Pendiente: estudio de convergencia (`notes/pd2_convergence.md`), no evaluado aquí.

## 2026-11-05 — Task 5, ronda de correcciones 1 (R22–R24)

- Ancho de pluma con borde del 10 % interpolado entre centros de celda (R22; mismos veredictos, desaparece el zigzag); rango de promedio de C2'' declarado (R23) con razones post hoc sobre rango común y cara completa (sermilik q0.1: 1.50, en el límite); teoría reevaluada con el ambiente evolucionado a x = 1 km como diagnóstico post hoc (R24: brecha de despegue 0.141 H → 0.059 H en sermilik q0.1). Se añade w de la celda de la cara (máx. sobre la pluma). Figura rehecha y `figures/fig1_caption.md` (borrador de pie de figura).

## 2026-11-05 — Task 6 Step 2: corrida corta cronometrada de (b) (PD1)

- `b_3d_twolayer_x0_40`, 30 min de modelo, np = 5, en `/Volumes/LaCie/portfolio_runs/exp1/b_3d_twolayer_x0_40_short30m{,_dt2}/`. `OBCS_CHECK: set-up OK` y "Inside Mask and OB locations: OK" con las filas de tierra (pregunta del Step 1 resuelta).
  1. dt = 1 s: OK, **192 s**, CFL máx. u/v/w 0.044/0.042/**0.128** (≤ 0.15 → se probó dt mayor), cg3d 190–203 iteraciones (tope 500), sin NaN, `pgrep` vacío.
  2. dt = 2 s: OK, **106 s**, CFL w 0.256, cg3d 215–237, sin NaN, `pgrep` vacío. Frente a dt = 1 s: derretimiento total en la cara −1344.8 vs −1344.0 kg/s (0.06 %), W rms de la diferencia 8e-4 m/s.
- Figura `notes/pd1_b_quicklook.png`: plumas de núcleo de 1 celda, w ≈ 1.3 m/s, llegan a ~30–50 m (pasan la interfaz de 100 m), salida entre 50 y 150 m. Sin ruido de 2Δ (G ambiente ≤ 0.10; alternancia de signo 7e-5).
- Costo de 14 casos × 3 h: 4.5 h con dt = 1 s; **2.5 h con dt = 2 s (propuesta)**. Opciones reducidas: 10 casos 1.8 h, solo dos capas 1.2 h, dz = 15 m ≈ 1.1 h, caja chica ≈ 1.0 h. Nota completa: `notes/pd1_b_twoplume3d.md`. **PD1 pendiente; campaña no lanzada.** `cases.yaml` sigue con dt = 1 s hasta la decisión.

## 2026-11-05 — Task 6: preparación de la campaña (b) tras PD1(b) (fallo R29)

- PD1(b) aprobado por Esteban: dt = 2 s, los 14 casos y animaciones GIF (R29). `cases.yaml` (b) `deltaT: 2.0` (5400 pasos por caso de 3 h); se volvió a generar `b_3d_sermilik_x0_160` para comprobarlo.
- `data.diagnostics` (b): listas 8 `animW` (WVEL 3D cada 120 s; pkg/diagnostics solo filtra niveles, no cortes horizontales) y 9 `animT` (THETA en los niveles 8 y 16, z = −75 y −155 m). Las listas 1–7 no cambian. Extra: ≈ 0.28 GB por caso (≈ 0.35 GB en exFAT), ≈ 3.9 GB en total (≈ 4.8 GB en disco). Solo salida: `code/` no cambia, así que la física es idéntica.
- `analysis/b_checks.py::eta_retention`: fracción retenida de la descarga (η̄ contra Q t / A). Se marca si en la última hora queda retenido más del 50 % (salud de la OB; justificación en `notes/pd1_b_twoplume3d.md` §6). **Aviso:** en las corridas cortas, a los 30 min la retención es del 69 % con dt = 1 s y del 110 % con dt = 2 s (la OB metió agua entre los 20 y 30 min). El derretimiento no cambia; es un transitorio de Orlanski. Revisar el primer caso de la campaña.
- `analysis/animate_b.py`: GIF por caso (w en la cara, ΔT a −75 y −155 m, reloj; escalas fijas, cmocean). Probado con la corrida corta (2 cuadros); un cuadro revisado a ojo.
- Pruebas: RED 4 fallos (dt en el render, listas de animación, `eta_retention` ×2) → GREEN; suite completa 166 PASS. Campaña **no lanzada** (corre conv_dx2p5).

## 2026-11-05 — Task 6: ronda de correcciones 1 de la preparación de (b) (R30)

- La revisión encontró que R30 no estaba implementado: el comando `'b_3d_*'` corría los 14 casos sin gate. Arreglo:
  - `run_campaign.py --skip-ok`: salta los casos que ya están OK (exige el mismo deltaT/endTime; si no coinciden, se niega a arrancar) y la compilación si `mitgcmuv` es más nuevo que `code/` y coincide en MPI.
  - `scripts/b_gate.py <run_dir>`: (i) |retained_last_hour| ≤ 0.5; (ii) `icefrontDiag` y THETA/SALT/UVEL/VVEL/ETAN en la iteración 900 idénticos byte a byte a la corrida corta `_short30m_dt2`. Sale con 0 (PASS) o 1 (FAIL, con los motivos y la diferencia máxima).
- Procedimiento en `notes/pd1_b_twoplume3d.md` §7: primer caso → gate (+ `--check-scales`) → solo con PASS, `'b_3d_*' --skip-ok`. La explicación de Orlanski para el 110 % a los 30 min queda como hipótesis hasta leer el gate.
- Menores: `eta_retention` también marca una retención < −0.5. `animate_b`: aviso si el GIF sigue pasando de `--max-mb`, `write_within_budget` probado, `ffmpeg_path` se fija antes de crear el writer, `--check-scales` (fracción de píxeles saturados por panel; en la corrida corta < 0.7 %). La prueba de la corrida FAILED verifica el mensaje.
- TDD: RED 5 fallos (`run_campaign`) y 10 fallos (gate, flag negativo, presupuesto del GIF, saturación) → GREEN; suite completa 181 PASS. No se corrió MITgcm.

## 2026-11-05 — Task 7: preparación de la campaña de (c) (sin correr MITgcm)

- `run_campaign.py`:
  - casos de 2 etapas (si `make_case` escribe `data.stage2`): si la etapa 1 está OK y existen `pickup.<n1>.data/.meta` y `pickup_orlanskiE.<n1>.*`, las entradas y los pickups se **copian** a `input_stage2` y la etapa 2 corre en `<run>/stage2`;
  - cada etapa conserva su `run_status.json`; `run_status_final.json` está OK solo si ambas lo están y last_iter = nIter0 + nTimeSteps;
  - `--set KEY=VALUE` para sobrescribir parámetros; `--skip-ok` usa el estado final.
  - Probado con un mitgcmuv falso a través del `run_mitgcm.sh` real.
- `make_case` (c): una corrida de viento que no pasa de la etapa 1 queda en una sola etapa con viento.
- `analysis/c_runs.py` (etapas, sección, Q_ex y calor, C5', W/Rd), `fig3_exchange_partition.py` (i = 100, T_ref = −1.90 °C) y `animate_c.py` (T en el nivel de salida con flechas, sección a lo largo del fiordo con la sección i = 100 y la esponja, Q_ex con marcador; escalas fijas; `--check-scales`).
  - Nueva lista de diagnóstico `animV` (VVEL nivel 19, medias de 3 h): ≈ 12 MB por caso.
  - GIF ≈ 0.07 MB por cuadro a tamaño real (sintético).
- W/Rd con los insumos reales de (c) = 0.667 (Rd = 8.39 km).
- Comandos del Step 5 y Step 6 en `notes/pd1_c_fjord3d.md`. El build de (c) sigue al día (no cambió `code/`).
- TDD: 19 pruebas nuevas, RED → GREEN; suite completa 200 PASS.

## 2026-11-05 — Task 7: ronda 1 de correcciones a la preparación de (c) (R31)

- R31, por rigor, antes del Step 6 (los diagnósticos no se recuperan después):
  - el calor hacia el frente sale de la media de 3 h de `ADVx_TH` (degC·m³/s, `gad_advection.F:841`), con H = ρ0 c_p Σ[−ADVx_TH + T_ref·ū·dA]: exacto e incluye u'T';
  - Q_ex sale de instantáneas horarias de UVEL (lista `snap1h` con UVEL y THETA 3D);
  - los productos de medias de 3 h quedan como `cross_check` en `fig3_metrics.json`.
  - Pruebas sintéticas: el calor es igual al producto con u y T constantes y difiere si covarían; Q_ex horario no sesgado frente a la media de 3 h con u que cambia de signo; se rechaza `mean3h` sin `ADVx_TH`.
- Disco: ≈ 5.8 GB por caso, ≈ 23 GB en total (LaCie tiene 3.6 TB libres).
- Menores:
  - `--set` rechaza claves desconocidas y conserva enteros;
  - `wind_stage1_s` debe ser múltiplo de dt, 3 h y 1 h;
  - `animate_c` sale con 1 si se pasa del presupuesto;
  - `c_runs.step_rate` (temporizador `MAIN LOOP`, sin arranques) para extrapolar el Step 5;
  - un solo helper `is_two_stage`;
  - `--check-scales` en los cuatro casos en la nota.
- TDD; suite completa 207 PASS. No cambió `code/` (el build sigue al día). No se corrió MITgcm.

## 2026-11-05 — (c): frontera este Orlanski -> Stevens (Task 7 fix)
- Motivo: la corrida de 1 día (`c_fjord_all_short1d`) vació el dominio con Orlanski (η̄ −3.03 m a 12 h, −11.38 m a 24 h; obc_E_uVel_Int +28 900 y +58 300 m³/s frente a Q0 = 100 m³/s). Esteban aprobó Stevens (2026-11-05), igual que en (b).
- `make_case._stevens_east` (antes `_b_stevens_east`) lo comparten (b) y (c). En (c): OBEu = Q0/A_ob (A_ob = 2.8e6 m², u_ob = 3.57e-5 m/s en `c_fjord_all`), OBEv = 0, con los **mismos nrec registros** que OBEt/OBEs, todos iguales: (c) usa `periodicExternalForcing` y pkg/obcs lee todos los archivos OB con la lógica periódica (`obcs_fields_load.F:102-105, 159-173`).
- `data.obcs`: `useStevensEast`, sin OBCS_PARM02 (solo se lee con Orlanski, `obcs_readparms.F:514-518`), `OBCS_PARM04` con Trelax = Srelax = 86400 s. `OBCS_OPTIONS.h`: `ALLOW_OBCS_STEVENS`. Costo conocido: el modo barotrópico queda cerrado en la frontera.
- `run_campaign.py`: el pickup OB que exige la etapa 2 se deriva de `data.obcs` (`required_ob_pickups`): `pickup_stevensE` si `useStevensEast` (`obcs_write_pickup.F:225-231`, `obcs_read_pickup.F:241-246`), `pickup_orlanskiE` si `useOrlanskiEast`.
- TDD; suite 220 PASS. Build `c_fjord3d` (--mpi) OK. No se corrió MITgcm (campaña (b) en curso).

## 2026-11-05 — (c) Stevens, ronda 2: la esponja T/S pasa de OBCS a RBCS
- Defecto crítico (revisión de 2d83953): con Stevens, OBEt/OBEs guardan TENDENCIAS de trazadores (`obcs_calc_stevens.F:744-752`, sumadas en `obcs_apply_ts.F:150-160`), pero la esponja T/S de OBCS las lee como VALORES objetivo (`obcs_sponge.F:545-550`). Esponja T/S apagada; U/V se queda (OBEu/OBEv son velocidades, `obcs_calc_stevens.F:222`, `obcs_sponge.F:136-150`).
- Zona de relajación rehecha con RBCS (Esteban, opción B): `make_case.sponge_zone_mask`, máscara = rbcs_tau (L−isl)/(L τ(isl)), reproduce la esponja no lineal vieja (`obcs_sponge.F:540-565`) con error relativo < 1e-12 frente a una reimplementación independiente. Constantes (10, 600, 21600 s) en `cases.yaml`, compartidas con `data.obcs`.
- `rbcs_T/S.bin` con nrec = ciclo/período registros (mismos tiempos que OBEt: offset 0, mismo GET_PERIODIC_INTERVAL, `rbcs_fields_load.F:72-76`); celdas de la pluma sin cambio; error si la zona se traslapa con ellas. Guardia: `useStevens*` con `OBCSsponge_Theta/Salt` TRUE levanta ValueError.
- `required_ob_pickups` acepta .TRUE./.T./T/TRUE/true, ignora `#` y `!`, y falla si el valor no se reconoce.
- Con Stevens la frontera saca exactamente Q0: η̄ debe quedar ≈ 0 y obc_E_uVel_Int ≈ +100 m³/s.
- 238 PASS (la suite pasó de 16 s a ~90 s: rbcs_T/S pesan 330 MB por caso (c)). Build al día, dry-run OK, sin correr MITgcm.

## 2026-11-05 — (c) Stevens, ronda 3
- Pruebas (c) con entradas pequeñas (`externForcingPeriod`=6 h, `shelf_period_days`=1, nrec = 4): la suite pasa de ~90 s y 39 GB de temporales a ~16 s y 0 B tras una corrida (`pytest.ini`: `tmp_path_retention_policy = failed`, `count = 1`). 241 PASS.
- Runner: `check_free_space` (libre >= max(20 GB, 4 × tamaño de la entrada)) antes de cada caso. Guardia `_check_stevens_sponge` trata la bandera ausente como TRUE (`obcs_readparms.F:280-281`). Comentario del desfase de un paso (RBCS carga en myTime, `load_fields_driver.F:249`; OBCS en futureTime, `obcs_calc.F:413`) y de OBEt = 0 (`obcs_sponge.F:545`).

## 2026-11-05 — Exp. 1 (b): análisis, Figura 2, C4/C4′

- `plume_diag.py`: `face_total_melt`, `face_melt_field`, `face_w_profile`, `ambient_far_from_face` (+ pruebas 3D sintéticas). `fig2_two_plume_melt.py` y `fig2_metrics.json`. Puntos de laboratorio digitalizados a mano (`analysis/data/cenedese_gatto_2016a_fig2.csv`).
- **C4: FALLA** (twolayer 0.209 H, sermilik 0.188 H, tolerancia 0.15 H; robusto a la elección del ambiente y al perfil inicial). El modelo no despega: w en la cara llega a la superficie. **C4′: PASA** (0.521 y 0.525).
- **Hallazgo grave:** 4 corridas sermilik (x0_40, 80, 160, 320) divergieron (cg3d sin converger a t ≈ 2.67 h, NaN antes de 2.83 h) pero `run_status` dice OK. La ventana 2–3 h no se puede aplicar a esos casos; no entran en ningún veredicto (C4/C4′ usan solo `_single` y `x0_0`, sanos).
- GIFs de 4 casos (≤ 2.6 MB). Saturación de escalas ≤ 0.8 %.

## 2026-11-07/08 — Exp. 1 (b): subconjunto a 5 m y animación de (a)

- Corridas a 5 m (twolayer): fondo sin descarga (11.3 h de reloj; cg3d ~480 iteraciones/paso frente a ~300 con pluma), x0 = 20, 0 y 80 m (~5.3–6.3 h cada una).
- x0 = 20: `MPI_Finalize` falló tras la última escritura (OFI sobre utun10 = Tailscale, que se cayó cuando la laptop pasó a batería), así que no salió la línea de fin normal. **R41:** Esteban la aceptó. Se agregó `run_accepted.json`, que `require_ok` respeta solo con esa causa y la corrida completa; `run_status.json` no se tocó (9ea236a).
- **Criterios OB 1–3: PASAN las cuatro corridas.**
- **R_E a 5 m:** 0.977 / 1.204 / 0.947 para x0 = 0 / 20 / 80, frente a 0.762 / 0.953 / 0.927 a 10 m. Se cumple la condición fijada de antemano (diferencia > 0.1 en x0 = 0 y 20): la interacción depende de la resolución y la Figura 2 se presenta a 5 m.
- A 5 m, las plumas cercanas suben casi hasta la superficie (z_despegue −9 y −17 m, frente a −55 m de la pluma única). R_E > 1 en x0 = 20 queda sin explicar. Es estable entre las dos mitades de la ventana (post-hoc).
- Commits: análisis en `excess_dx5.py` (0d24952); Figura 2 con panel (b) y README (7127d06).
- Animación de (a): la nueva corrida es idéntica bit a bit a la de producción de 1 a 6 h (diferencia 0.0). GIF cada 5 min y sección «Animations» en el README (92f033f). 293 tests PASAN.

## 2026-11-08 — Exp. 2: preregistro, configuración, smoke y primera calibración (Tasks 0–4a)

- **PD0:** criterios K1, V1, V2a–d, V3, N1–N4, H1 y H2 fijados en `notes/pd0_preregistro.md` antes de cualquier corrida.
- **Tasks 1–3:** geometría desde el archivo GFZ para fondo y base del hielo (R7, regla de mayoría flotante R9), columna mínima de 2 celdas y conectividad; ley de fusión AD16 en `code/` (override de `shelfice_thermodynamics.F`), trazador de deshielo con fuente virtual, librería de métricas (`analysis/metrics.py`) con `require_ok`.
- **Smoke de 30 días:** S1, S2, S4, S5 pasan; **S3 FALLA** (T y S fuera del rango inicial; diagnóstico post hoc con la línea de Gade, R10). La revisión encontró que el dominio era periódico en y (R12): fila de tierra al norte y smoke repetido con paredes.
- **Task 4a:** cronometraje y dt por CFL (400 s a 2 km, 800 s a 4 km; R17); `checkIniTemp=.FALSE.` (R18). Primera cola de calibración (20:40–23:46 CST): 30 m a⁻¹ no quedó acotado con ΓT ∈ {0.015, 0.03, 0.06} → PD2.
- Todo esto corrió con la forma errónea de u\* y quedó **anulado** al día siguiente (R21); se conserva como registro.

## 2026-11-08 (tarde): correcciones de la revisión final (PD4)

- **Bloqueantes:**
  - Tabla de criterios completa (C1–C6, con resultado). C4 es ±0.15 H, no 0.10 H.
  - C4′ pasa solo por construcción.
  - C6 excedido: ~49 h sumadas en 74 corridas, con extensiones aprobadas.
- **(a):** `pd2_convergence_eval` usa la extensión pre-registrada `sens_xu2012_36h`, con ventana de 24–36 h.
  - R_Q pasa de 0.547 a 0.720 y b de 20.0 a 9.8 m.
  - El déficit de velocidad no es de resolución, pero depende del cierre.
  - La pluma es más ancha que el sombrero de copa; no es un chorro estrecho.
- **(b):** a 5 m la pluma se desprende por encima de la interfaz (−100 m), no por debajo. Se corrigió también en `pd1_b_twoplume3d.md`.
- **Etiquetas y límites:**
  - D1 (lectura del afloramiento) queda etiquetada como post hoc.
  - §4.1 acotada.
  - Nueva limitación: OB barotrópicamente cerrada en (c).
  - Conclusión 6 reescrita.
- **Citas verificadas contra los PDF:**
  - Cenedese & Gatto: la referencia es 2× una pluma con la descarga por fuente, igual que la nuestra.
  - Kimura: la comparación es contra una sola pluma con la misma descarga total.
  - Jackson 2018: lo que llaman simplificación es la condición en la boca; 30 m es «típico de Sermilik» y su control usa 14 m.
  - Carroll (viento) y Sanchez (Sermilik, la pluma duplica el intercambio) quedan separados.
- `.gitignore` sin entradas del asistente en todo el historial (reescrito; las entradas se movieron a `.git/info/exclude`). Las referencias de commits en las notas están actualizadas.

## 2026-11-09 (madrugada) — Exp. 2: calibración ampliada y serie diagnóstica R20 (anuladas por R21)

- **Enmienda de PD2 (R19):** se añadió ΓT = 0.0075. Segunda cola (00:41–03:13 CST): ΓT_cal = 0.0094686 (4 km) y 0.0099627 (2 km); K1 pasaba y V1 fallaba (< 0.011).
- **R20:** serie diagnóstica post hoc a 4 km (advección 33, hFacMin 0.25, `SHI_withBL_uStarTopDz`, ajuste convectivo), una opción por corrida, para explicar el ΓT_cal bajo. Ninguna lo explicó.
- Hacia las 08:05 CST se encontró el error de u\* (siguiente entrada): la calibración, la enmienda y la serie R20 quedaron anuladas.

## 2026-11-09 — Exp. 2: error (a) en la forma de u\* (R21) y relanzamiento de la calibración

- **Error (a), nuestro:** AD16 Eq. 27, verificada en la página del PDF (GMD 9, p. 2485 = p. 15 del PDF), es u\*² = Cd·(u_w² + u_tidal²), con piso √Cd·u_tidal = 5e-4 m s⁻¹. Implementamos Cd·u_w² + u_tidal² (piso 0.01 m s⁻¹). La fórmula del plan venía de una extracción de texto sin paréntesis.
- **Evidencia:** u\* medio profundo de MITgcm-JPL en Ocean0 = 0.0040 m s⁻¹; el nuestro, 0.0109 m s⁻¹.
- **Corrección con TDD:** prueba con anclas a mano (5.0e-4 y 7.0711e-4 m s⁻¹) que falla antes y pasa después; `melt.py`, override de Fortran (mínimo, marcado, con cita de la página), `metrics.ustar_check` (piso √Cd·u_tidal) y sus pruebas; prueba nueva que exige que Cd multiplique a los dos términos en el Fortran. Suite completa: 496 pasan.
- Binarios `ocean0_2km` (MPI), `ocean0_4km` y `ocean0_4km_mpi` recompilados.
- **Anuladas** (se conservan en disco): smoke, `cal2_*`/`cal4_*` (todos los ΓT), `ocean0_2km`, `ocean0_4km`, `tim4m_*`, `diag4_*`. La enmienda de PD2 (ΓT = 0.0075) y la serie R20 estaban condicionadas a la ley errónea: quedan como registro, no se reutilizan. K1, V1, N1–N3 y pd0 no cambian.
- **Relanzamiento** de la calibración de pd0 K1 desde cero, ΓT ∈ {0.015, 0.03, 0.06} en las dos mallas, con sufijo de corrida `_r21` (`scripts/queue_calibration.sh`). Si 30 m a⁻¹ no queda acotado en una malla, la cola se detiene ahí (PD2), sin extrapolar.
- **Resultados R21** (cola 08:11–12:28 CST): ΓT_cal = 0.029356 (4 km) y 0.030959 (2 km), acotados con los 3 puntos de pd0. K1 PASA en las dos mallas (30.02 y 30.00 m a⁻¹), V1 PASA (0.030959 ∈ [0.011, 0.2]), N1–N3 PASAN en las 8 corridas, CFL máximo 0.348. Detalle en `notes/pd1_calibracion.md`, «Resultados R21».
- **Candado:** `melt_law` en `case.json`; `--skip-ok` y `metrics.load_run` rechazan corridas con la ley anulada; sufijo por omisión `_r21`.
- **Comparación post hoc con MITgcm-JPL** (`notes/pd1_comparacion_jpl.md`): melt profundo 30.00 contra 31.35 m a⁻¹, u\* profundo 0.00365 contra 0.00397 m s⁻¹, correlación de los mapas de melt r = 0.91.
- **Producción (PD1a aprobada):** ΓT_cal en `cases.yaml` (4 km 0.029356 para el ensamble y `ocean0_4km`; 2 km 0.030959 para Ocean0 y Ocean1); `binary_sha256` en `case.json` al correr; `scripts/queue_production.sh` (ensamble de 28 casos, N1–N3, luego Ocean1 en 4 etapas de 5 años) lanzada.

## 2026-11-09 (noche) — Exp. 2: ensamble a 4 km, evaluación de H1 y H2 (Tasks 5 y 6)

- Cola de producción (desde las 16:07 CST, con `caffeinate` desde las 20:10 tras ver pausas por reposo del Mac): ensamble de 28 corridas y referencia; fin a las 22:03 CST. Aclaración R22 (⟨m⟩ de H1 = melt de toda la plataforma) subida a las 16:41, antes de leer ninguna salida.
- **N1–N3** (`ensemble_n_check.json`, 22:04:40): N3 y N2 pasan en las 29 corridas; N1 falla en 5 (casos fríos o de termoclina profunda), que se quedan en los ajustes y se señalan.
- **H1 FALLA** en los tres z_c: R² de la cuadrática ≥ 0.9995, pero n = 2.43, 2.78 y 2.70 > 2.0. **H2 PASA** en los tres z_c (ρ = −1). Figuras 2 y 3, y GIF del trazador (commit 8cc4547, 22:15 CST).
- Ocean1 a 2 km arrancó a las 22:04 (fin el 2026-11-10, 06:45).

## 2026-11-10 — Exp. 2: Ocean1 a 2 km, Figura 1 y comparación post hoc con JPL (Task 4 fase B)

- `ocean1_2km_r21` terminó (20 años, 4 etapas, 06:45 CST). Cargado con `metrics.load_run`: 240 registros contiguos, costuras de etapa sin huecos.
- **pd0, año 20:** V2a PASA (10.34 m a⁻¹), **V2b FALLA** (overturning 0.185 Sv < 0.20), V2c PASA (0.258 Sv), V2d PASA (0.492 Sv), V3 barotrópica PASA, **V3 overturning FALLA**, N1 (0.11 %), N2 (1.3e-6) y N3 PASAN. CFL máx. 0.23. V2b y V3 overturning van a PD2; clasificación pendiente. Detalle en `notes/pd1_ocean1.md`.
- Figura 1 (`analysis/fig1_ocean0_validation.py`, `figures/fig1_ocean.{png,pdf,json}`). El pie no usa el acuerdo de melt 2 contra 4 km como evidencia independiente.
- Comparación post hoc con MITgcm-JPL Ocean1: melt profundo 30.37 contra 31.55 m a⁻¹, u\* profundo 0.00371 contra 0.00396 m s⁻¹, overturning 0.185 contra 0.203 Sv, r de los mapas = 0.91.
- **PD2 (Esteban, 2026-11-10):**
  - V2b sigue como FALLA pre-registrada (0.185 Sv en la cavidad). Se añade una fila **post hoc** con la definición de Y26 (máximo en todo el dominio): 0.222 Sv, dentro de 0.20–0.25; JPL 0.234. Error (a) nuestro al operacionalizar el criterio: x < 640 km es de la Fig. 13.
  - V3 overturning: FALLA marginal, sin clase (4.9 % bajo la banda extendida, 1.8 % bajo la fijada en el anclaje de 10).
  - Figura 1 regenerada, con los tramos extrapolados de la banda rayados.

## 2026-11-10 — Exp. 2: análisis post hoc de H1 (R23, R24) y predictores cuadráticos

- **Contexto:** H1 FALLA en los tres z_c (n = 2.43, 2.78 y 2.70 > 2.0; la cuadrática ajusta con R² ≥ 0.9995) y H2 PASA (ρ = −1). Todo lo de abajo es **post hoc**; el veredicto preregistrado de H1 no cambia. Detalle en `notes/pd2_ensamble.md`.
- **R23**, forzamiento en la base (`analysis/h1_posthoc_tf.py`, `figures/fig2_posthoc_tf.*`):
  - Con el forzamiento a priori TF_eff, n dentro de cada z_c sale aún más alto: 2.36, 3.44 y 6.49.
  - En el conjunto, n = 1.55 con R² = 0.85.
  - Con T_w − T_b diagnosticado, n ≈ 1.45, pero esa medida es en parte circular.
- **Predictores cuadráticos** de Favier 2019 (Ec. 4 y 5) y Burgard 2022 (Ec. 14 y 18) (`figures/fig2_posthoc_quad.*`):
  - El no local es idéntico a TF_eff².
  - El local da q = 1.18–1.83 dentro de cada z_c, y 1.02 en el conjunto. Esta última pendiente mezcla la variación entre z_c y no se toma como evidencia.
- **R24**, prueba de exposición con perfil tipo Holland 2008 (z_c = 150 m, d = 100 m):
  - La predicción y la regla se subieron antes de correr (874205e, 11:00:55).
  - Se corrieron 5 casos `hol_T*` a 4 km, de 11:01 a 11:34, con `caffeinate`. N1–N3 pasan.
  - n = 2.16 (LOO 2.11–2.19): **no apoya**, según la regla.
  - Solo ~50 % de la base queda por debajo de 200 m.
  - Dentro de `hol`, n contra TF_eff = 1.98 y q local = 0.99 (`figures/fig2_posthoc_r24.*`).

## 2026-11-10 — Exp. 2, Task 7: exportación para el Exp. 3, README, reproducibilidad

- **Exportación para el Exp. 3** (`analysis/export_exp3.py`, TDD en `tests/test_export.py`: forma, NaN en celdas secas, exclusión de las corridas que fallan N3 o faltan, partición, N1, perfiles, escritura reproducible byte a byte). Archivo fuera del repo: `/Volumes/LaCie/portfolio_runs/exp2/exp3_dataset.npz`, 29 corridas (26 de entrenamiento, 3 de prueba con Tbot = 2.0), 0 excluidas, 5 marcadas por N1; sha256 `3d2eceb1…`. N1 y melt coinciden exactamente con `ensemble_table.json`. Esquema en `notes/exp3_dataset.md`.
- **README del Exp. 2** en estructura IMRaD, en inglés, renderizado desde `README.template.md` con `scripts/render_readme.py` (mismo mecanismo `{{name}}` del Exp. 1, con su propio `values()`). Todos los números salen de `figures/*.json`; la única excepción es el u\* de 0.0109 m s⁻¹ de la corrida anulada, citado de `pd1_config.md` y señalado así en el texto. Tabla de criterios con resultado: V2b FALLA (con la fila post hoc), V3 overturning FALLA marginal, H1 FALLA, H2 PASA, 5 fallos de N1. `render_readme.py --check` pasa; pruebas en `tests/test_exp2_render_readme.py`; suite completa: 569 pasan.
- **Registros anulados:** los JSON anteriores a R21 (`calibration_*`, `*_3pt`, `verify_*`, `pd1_calibration_eval.json`, `diag_bl_*`) pasaron con `git mv` a `figures/void_pre_r21/`, con un README. El renderizador se niega a leerlos. Rutas actualizadas en `pd1_calibracion.md` y `pd1_diagnostico_capa_limite.md`.
- **`reproduce.sh`:** `--quick` regenera calibración, evaluación, comparación con JPL, Figuras 1–3, análisis post hoc, la exportación y el README; `--gifs` añade la animación; el modo completo lanza las colas en su orden (calibración → producción → R24) con `caffeinate`. Con `SOURCE_DATE_EPOCH=0` (valor por omisión), dos `--quick` seguidos dan archivos idénticos (shasum de 27 archivos y del .npz). Frente a lo ya subido, solo cambian los PDF, y solo en su fecha de creación; JSON, PNG, GIF y README salen idénticos.
- README principal y `RUNNING.md` con el Exp. 2 (rutas `EXP2_RUNS`/`EXP2_BUILDS`, colas, `reproduce.sh`).

## 2026-11-10 (tarde) — Exp. 2: ronda final de correcciones antes del merge y prueba de reinicio (R25)

- **R25, reinicio bit a bit de Ocean1** (expectativa subida antes de correr, ad882cd): 2 meses a 2 km, np = 5, continua contra 2 etapas de 1 mes con pickup. Un primer intento quedó FAILED por la regla del runner (monitor diario que no cae en el último paso; el modelo terminó bien); repetido con monitor mensual. **Idéntico bit a bit** (máx. |dif| = 0.0 en pickup, trazador y estados finales), con el ejecutable de producción.
- **V3 barotrópica:** se reporta como PASA marginal (2.3 % dentro de la banda extendida, 4.5 % con la banda en el ancla), simétrica con la FALLA marginal de la overturning; `fig1_ocean.json` da ahora el margen interior y la fila post hoc de V2b usa `inside_range`.
- **Rutas:** los scripts de análisis toman `EXP2_RUNS` de `pfcommon.paths`; `jpl_compare.py --runs-root`; prueba nueva. `RUNNING.md` aclara que las colas fijan la ruta de LaCie.
- **README:** R9 declarado (s escala solo el hielo flotante; el borde mojado se mueve por la columna mínima), menos repeticiones de Holland 2008, §3.4 más corto, referencias de Rosevear 2022, Yung 2025 y Gade 1979 verificadas en Consensus (ficha bibliográfica; también en `readinglist.md`), frases matizadas (N1, sobreenfriamiento, Gade en la corrida anulada).
- Notas: §7b de `pd1_config.md` marcado como anulado por R21; lectura de R24 en `pd2_ensamble.md` suavizada. `.gitignore` con `*.nc`, `*.log`, `*.npz`. Entradas de la bitácora del 2026-11-08 y de la noche del 2026-11-09 completadas.

## 2026-11-10 (noche) — Exp. 2: Figura 0 (contexto, esquema de procesos, geometría)

- **(a)–(b) satélite:** MODIS Terra (GIBS, 1 dic 2025) del embayment del mar de Amundsen con localizador Blue Marble y el dominio ISOMIP+ a la misma escala (contorno leído del archivo de geometría); Landsat 9 C2 L2 (USGS, dominio público, vía Planetary Computer) de la plataforma de Pine Island, reproyectado a una estereográfica polar girada para que el hielo fluya a la derecha, como en el esquema. `figures/src/fetch_sat.py` documenta URL, escena, fecha y licencia.
- **(c) esquema VectorCraft** con texturas (grano, texturizer, roughen, zigzag), procesos 1–6 numerados y recuadro de la capa límite con las Ec. 24–27 de AD16; z_c y Tbot marcados como parámetros de diseño. Los dos regímenes de salida del agua de deshielo (Holland et al. 2008, §4) se dibujan como posibilidades.
- **(d) render Blender** de la geometría real (borrador y lecho, VE ×100, mitad cercana recortada por la línea central, simétrica a < 8 m), anotado con VectorCraft en puntos proyectados desde la cámara; proceso 7 (celda ciclónica intensificada en y > 50 km, Yung et al. 2026) y recuadro de los perfiles rampa. Error propio encontrado y corregido: la matriz de la cámara no se actualizaba antes de proyectar las etiquetas.
- Ninguna figura muestra resultados del modelo. README con la Figura 0 y nota de reproducción; `reproduce.sh` la redibuja solo con `VECTORCRAFT_CLI` (y `BLENDER` para d). PDF recomprimido con Ghostscript (solo las imágenes) para quedar bajo 1 MB; el SVG ya no se escribe.

## 2026-11-10 (noche, 2) — Exp. 2: introducción reescrita y correcciones de la Figura 0

- **§1 reescrita al estilo de artículo** (sin subsecciones): lo que está en juego (Rignot et al. 2013; Holland et al. 2008) → el CDW y la termoclina (De Rydt 2014, Rintoul 2016, Nakayama 2023, Jenkins 2018, Park 2024) → recorrido de la Figura 0c (procesos 1–6) → lo conocido y el hueco en sensibilidad (Holland 2008, Favier 2019, Burgard 2022, De Rydt 2014, Hattermann 2014) y en el destino del agua de deshielo (Sheehan 2023) → ISOMIP+ y Figura 0d (proceso 7) → "aquí" con H1 y H2 neutrales. Se quitó la cita a la página ARIES; §1.2 → §1 en dos referencias cruzadas. Cinco referencias nuevas (Consensus, verificadas en el resumen o en `literatura/`); todos los resultados de las búsquedas van a `readinglist.md`.
- **Figura 0a:** los puntos de Pine Island y Thwaites coinciden con las entradas US-ACAN del SCAR Composite Gazetteer (placenames.aq); se rotula la segunda zona de agua abierta como Amundsen Sea Polynya (73.0–74.2° S, 110–118° O, al norte del Dotson; Jenkins et al. 2018).
- **Figura 0d:** el rótulo de la pared y = 80 km apunta ahora a la cara de la pared en el océano abierto; el hielo apoyado se dibuja recortado a 200 m sobre el nivel del mar (encuadre de cámara sin cambio), anotado en la figura y en el pie.

## 2026-11-10 (noche, 3) — Merge del Exp. 2, correos y revisión de literatura del Exp. 3

- **Merge del Exp. 2 en `main` y push** (ab87c33), con la aprobación de Esteban (PD3). Se verificó que el link `tree/main/exp2_isomip_cavity` responde.
- **Correos a Nakayama y Sheehan** actualizados con el Exp. 2, fuera del repo. Sheehan et al. (2023) se leyó completo: la frase citada está en sus *Key Points*. Esteban enviará los correos después del Exp. 3.
- **Revisión de literatura del Exp. 3** (`exp3_melt_emulator/notes/rev_literatura.md`): 11 búsquedas en Consensus, ~153 artículos únicos, 30 citados.
- **Hallazgo principal:** DeepMelt-GL (Ockenden et al. 2026) ya publicó el ensamble de redes como indicador de confianza y la sensibilidad fuera de la distribución de entrenamiento. El Exp. 3 lo reproduce y lo cita, no lo reclama.
- **Aporte propuesto:**
  - la cuadrática local frente a la no local en un caso donde la cuadrática falla de forma conocida (H1 del Exp. 2);
  - una comparación con líneas base simples: regresión lineal, proceso gaussiano y las dos cuadráticas, siguiendo a Lütjens et al. 2024 y a Myren & Lawrence 2021.
- `readinglist.md` con 18 entradas nuevas en la sección del Exp. 3.

## 2026-11-11 — Exp. 3: validación cruzada y prueba final (evaluada una sola vez)

- **Validación cruzada:**
  - Comando: `run_protocol.py --cv --wall-limit-hours 20 --amendment "Enmienda 1 del PD0 …"` (entorno anemoi, 5 procesos, 1 hilo cada uno).
  - El primer intento (19:20 del 10-nov) se detuvo a las 21:04 por la regla congelada de 8 h, con una proyección de 14.1 h y 9 unidades escritas.
  - Esteban aprobó la Enmienda 1: el presupuesto de reloj pasa a 20 h y el protocolo científico no cambia.
  - Se relanzó a las 22:06 y terminó a las 09:14 del 11-nov: 476/476 unidades, 892 entradas, 40 100 s en esta sesión.
  - Resultados en `cv_results.json` (4b020d0). Las 9 unidades anteriores a la enmienda no llevan la clave `amendment`; el registro de la enmienda está en el nivel superior del JSON.
- **Revisión antes de la prueba final:** había un hallazgo importante: una caída durante `--final` habría gastado la prueba sin dejar registro. Se corrigió en 4b6e3a7:
  - cualquier falla queda escrita en el marcador;
  - `--resume-after-crash` solo se acepta después de una falla registrada;
  - todos los intentos aparecen en `final_test.json`.
  - La re-revisión no dejó hallazgos.
- **Prueba final:**
  - Comando: `run_protocol.py --final`, intento 1, de 09:42 a 10:33, 3061 s (M7 tomó 3052 s).
  - Se entrenó con las 25 corridas y se predijeron las 3 corridas con Tbot = 2.0.
  - Resultado en `final_test.json` (9d79dae). Mapas en LaCie (`final_maps.npz`, sha 20ada7e54cbc).
- La evaluación formal de H3a–H3c y las figuras van en la Tarea 9.

## 2026-11-11 — Exp. 3: diagnóstico post hoc de `profile_test` (M7 en `ocean0_4km`)

- **Pregunta:** si el +413 % de M7 en `ocean0_4km` es un error nuestro o un comportamiento fuera de distribución. Nota: `exp3_melt_emulator/notes/diag_profile_test.md`; números: `figures/profile_test_diag.json`; script: `analysis/profile_test_diag.py`.
- **Comprobaciones de entrada:** el perfil, el TF, el `draft`, la máscara y las unidades coinciden exactamente. El reajuste de M7 reproduce bit a bit los integrados de la unidad. **No se encontró ningún error.**
- **Mecanismo:**
  - el 58 % de las celdas de Ocean0 tiene un TF mayor que cualquier TF del entrenamiento a la misma distancia de z_c, porque d = 720 frente a 300 m;
  - en esas celdas, g_θ ≈ 4, cuando haría falta ≈ 0.14;
  - el 84 % del exceso integrado está en ellas, y por debajo de z_c − 50 m el error integrado por banda es de −6 % a +1 %.
- **Clasificación:** (c) elección del método, con un componente de (b). No cambia ningún resultado preregistrado.

## 2026-11-11 — Exp. 3: README, `reproduce.sh` y documentación general (Tarea 10)

- **README del Exp. 3** en formato IMRaD, generado por `scripts/render_readme.py` desde `README.template.md` y los JSON de `figures/` (mismo mecanismo que el Exp. 2). Una prueba nueva (`test_every_number_from_json`) falla si la plantilla tiene un dígito fuera de los marcadores, salvo etiquetas (M0–M7, H3a–c, Ec., §, Figura), años de citas y dos constantes de diseño del Exp. 2 (4 km, d = 300 m). `--check` sale con 0.
- **Resultado principal, como estaba preregistrado:** ganan los métodos simples. Por `rmse_map`, M4 (regresión por celda) gana en `loro`, `warm_out` y la prueba final, y M2 en `zc_out`, donde M4 es el peor. Por \|error integrado\|, M6 gana en `loro`, `warm_out` y la final. H3a FALLA, H3b PASA, H3c PASA para M6 y FALLA para M7 (2 casos). DeepMelt-GL se cita como el resultado que se replica; no se dice "el primero".
- **Matices escritos en el README:** el diagnóstico post hoc de `profile_test` (M7 +413 %, M2 −4.8 % con signo, M6 +53 %; M7, como M4/M5, hereda d = 300 m por `draft_minus_zc`); las advertencias de la Figura 3c (REF y +1 dentro de la muestra, DA sobre (+1, +2) como extensión de LB25 Ec. 3, el punto de M4 a z_c = 550 mal condicionado, umbrales de despliegue post hoc); los ajustes de potencia de (b) como post hoc respecto al Exp. 2; la sección de divulgación del PD0, la Enmienda 1 y la prueba final evaluada una vez. La sensibilidad sin las corridas N1 cambia el orden por mapa dentro del rango (M4 pierde la ventaja): se reporta como tal.
- **Hallazgos al redactar:** (1) Favier et al. 2019 Ec. 4 sí lleva (ρ_sw c_p/ρ_i L_i)² según el PDF (p. 2259), contra lo que dice la Decisión 5 del PD0; el README solo atribuye la forma de M2 a Burgard 2022 Ec. 14 y no afirma nada sobre la potencia de c en Favier. Pendiente de que Esteban lo revise; el resultado no cambia porque K absorbe la constante. (2) El PD0 pedía intervalos *bootstrap* sobre las corridas y no se calcularon: se declara como limitación. (3) Lütjens et al. es de 2025 (JAMES 17), no 2024: corregido en el README y en `readinglist.md`.
- **`reproduce.sh`:** `--quick` regenera `hypotheses.json`, las Figuras 1–4 y el README desde los JSON (con LaCie montado los JSON y PNG salen idénticos byte a byte; los PDF solo cambian la fecha de creación); el modo completo recalcula el piso de ruido en un archivo aparte y lo compara, corre `--cv` con la Enmienda 1 y el diagnóstico post hoc, y nunca `--final`; comprueba al final que `final_test.json` y el marcador no cambiaron.
- **Otros:** `pytest_anemoi.sh` añade la ruta por defecto también cuando solo se pasan opciones (`-q`); `RUNNING.md` documenta el entorno anemoi, el shim de pytest, los modos de `run_protocol.py` y `reproduce.sh`; el README general tiene la entrada del Exp. 3; `readinglist.md` con Berdahl et al. 2021 (*The Cryosphere* 15:2683, DOI verificado en la copia local).
- Pruebas: 145 pasan y 1 se omite bajo anemoi; 648 pasan y 3 se omiten bajo pyocean.

## 2026-11-11 (tarde) — Exp. 3: correcciones del README, intervalos *bootstrap* y errata de Favier (ronda 1 de la Tarea 10)

- **Intervalos *bootstrap* (R17).** Estaban preregistrados (PD0, "Métricas") y no se habían calculado. `analysis/bootstrap.py` remuestrea las corridas de prueba con reemplazo, sin reajustar nada, con los mismos índices para todos los métodos: B = 10 000 y semilla 0 (elección del analista, hecha al analizar, porque el PD0 no los fijaba). Resultado en `figures/bootstrap.json`; intervalo percentil del 95 % de la media por método de `rmse_map`, \|error integrado\| y `rmse_gl`, por tipo de pliegue y variante. En `loro`, `warm_out` y la final los intervalos de M4–M7 se solapan; solo en `zc_out` el de M2 queda entero por debajo del de M4. Los de la final (3 corridas) son indicativos. H3a–H3c no cambian: se definen sobre medias.
- **Matiz N1 en el resumen, la conclusión 1, §4.1 y el README general:** sin las 5 corridas N1, M4 pierde la ventaja (en `loro` gana M6, en `warm_out` M5).
- **`reproduce.sh`, modo completo:** `--cv` y el diagnóstico de `profile_test` reescriben `cv_results.json` y `profile_test_diag.json` con metadatos de sesión nuevos (el README habría dicho "0.0 h"). Ahora `scripts/rerun_compare_restore.sh` respalda el archivo, re-ejecuta, compara con `analysis/compare_records.py` (ignorando tiempos) y restaura el registro versionado. Con pruebas sintéticas.
- **Errata R16:** Favier 2019 Ec. 4 sí lleva (ρ_sw c_po/ρ_i L_i)² (PDF local, p. 2259). El README cita M2 como F19 Ec. 4 = B22 Ec. 14; el PD0 tiene una sección "Erratas" con la línea fechada (la Decisión 5 no se reescribe). El docstring de `QuadLocal` ya era consistente.
- **Menores:** §2.6 completa (ajustes dentro de la muestra, n y K de M7, H1 del Exp. 2); orden de `zc_out` rotulado; "sin override" precisado (`--resume-after-crash` tras una falla registrada); tope k ≤ m − 1 de M5 y Berdahl 2021 §2.3 citado por el *nugget*; `pytest_anemoi.sh` respeta `PYOCEAN`/`ANEMOI`; `-h` de `reproduce.sh`; más `assert` sobre las frases del README; Lütjens 2025 en `rev_literatura.md`; marcas ✅ de `readinglist.md` para los PDF ya en `literatura/`.
