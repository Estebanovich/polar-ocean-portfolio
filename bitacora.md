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
