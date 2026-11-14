# PD1 — Lectura de `pkg/seaice`, `pkg/kpp` y `pkg/ggl90`: qué esfuerzo recibe el océano (Task 4)

Fecha: 2026-11-12. Fuente: `/Volumes/ESTEBAN_MIT/MITgcm` (solo lectura; commit 8e1c407, 2026-01-12). Todas las citas son `archivo:línea` de ese árbol. Configuración que se lee: `config/fjord_vmf/code` (C-grid, `SEAICE_CGRID` definido y `SEAICE_BGRID_DYNAMICS` sin definir en `SEAICE_OPTIONS.h:90-95`; `SEAICE_ALLOW_FREEDRIFT` definido para F4), sin `pkg/exf` y sin termodinámica.

**Decisión del PD1 (Esteban, 2026-11-12):** se aplican la opción 2 de §6 (`DWATN` sin reología) y la opción 1 de §7 (`obcs_check.F` relajado con `DISABLE_SEAICE_OBCS`). F4 usa el mismo polígono que F2 (ruling R15). La implementación y su verificación están en §8. Las secciones 1–7 describen MITgcm sin modificar y siguen valiendo como diagnóstico de por qué hizo falta el cambio.

**Resumen para el PD1.** Los dos modos de hielo del spec (§3.5) **no hacen lo que el spec dice** en el esfuerzo que recibe el océano:

| Modo (spec §3.5) | Lo que el spec espera | Lo que hace MITgcm (lectura del código) | Comprobación empírica |
|---|---|---|---|
| Hielo fijo, `SEAICEuseDYNAMICS = .FALSE.` (F2) | "el océano siente solo el arrastre hielo–agua, sin viento": τ = ρ₀ C_w \|0 − u\|(0 − u) | **τ = 0** bajo hielo de A = 1: ni viento ni arrastre. Es una tapa sin fricción. | confirmado: τ = 0.0 en 13 072 puntos U durante 3 h (§5) |
| Hielo a la deriva, `SEAICEuseFREEDRIFT = .TRUE.` (F4) | losa con inercia y arrastre: el hielo se mueve con el viento y arrastra al océano | el hielo se mueve (deriva libre), pero **τ = 0** bajo A = 1: el océano no recibe ni el viento ni el arrastre del hielo | confirmado: hielo a 0.2 m s⁻¹ y τ = 0.0 bajo A = 1 (§5) |

La causa es la misma en los dos casos: el coeficiente de arrastre `DWATN` vale 0 al empezar y solo se calcula dentro de los resolvedores de la reología (LSR, EVP, Krylov, JFNK), que ninguno de los dos modos llama. No lo arreglo en silencio (instrucción del controlador). Las opciones van en §6, y es el punto **PD1-S1**.

**Segundo problema, que bloquea (descubierto al correr; §7):** F2 y F4, tal como están configurados, **no arrancan**.
- `pkg/obcs` rechaza las fronteras de Orlanski si `useSEAICE` está activo (`pkg/obcs/obcs_check.F:400-405`).
- También exige `pkg/exf` para cualquier corrida con `seaice` y `obcs` (`:571-579`).
- Las dos comprobaciones son incondicionales.

Es el punto **PD1-S2**.

KPP y GGL90 sí ven el esfuerzo que deja `pkg/seaice` (§3): el problema no está en el cierre, sino en el esfuerzo que llega.

---

## 1. Qué recibe el océano: el camino de `fu, fv`

1. **Viento.** `external_fields_load` interpola `fu, fv` de `zonalWindFile`/`meridWindFile` en cada paso (`model/src/external_fields_load.F:69-72` para el intervalo, `:114-125` para la lectura, `:268` para la interpolación de `fu`). En la etapa libre no hay archivo de viento y `fu, fv` no se recargan.
2. **`pkg/seaice` corrige `fu, fv` en el sitio** antes de que los use el océano. El orden de llamadas en `model/src/do_oceanic_phys.F` es: `OBCS_CALC` (`:322`), `SEAICE_MODEL` (`:453`), `EXTERNAL_FORCING_SURF` (`:579`), `KPP_CALC` (`:956`), `GGL90_CALC` (`:1010`) y el diagnóstico `oceTAUX/oceTAUY` (`DIAGS_OCEANIC_SURF_FLUX`, `:1119`).
3. `SEAICE_MODEL` llama siempre a `SEAICE_DYNSOLVER` (`pkg/seaice/seaice_model.F:174-188`; el comentario de `:174-176` dice que el esfuerzo superficial se actualiza "even if ... SEAICEuseDYNAMICS = F").
4. `SEAICE_DYNSOLVER`:
   - calcula el esfuerzo del viento sobre el hielo con `SEAICE_GET_DYNFORCING` en cada paso (`pkg/seaice/seaice_dynsolver.F:129-132`);
   - **solo si `SEAICEuseDYNAMICS`** (`:135-137`) resuelve el momento del hielo: deriva libre (`:303-321`) o LSR (`:343-346`);
   - siempre, si `SEAICEupdateOceanStress` (por omisión `.TRUE.`, `seaice_readparms.F:265`), llama a `SEAICE_OCEAN_STRESS` (`:374-385`).
5. **`SEAICE_OCEAN_STRESS`** (`useHB87StressCoupling = .FALSE.` por omisión, `seaice_readparms.F:264`), en los puntos U y V (`pkg/seaice/seaice_ocean_stress.F:105-140`):

   fu ← (1 − A_u) fu + A_u · ½(DWATN_i + DWATN_{i−1}) cos θ_w (u_hielo − u_1) (+ término de giro, θ_w = 0 por omisión, `seaice_readparms.F:546`),

   con A_u = ½(A_i + A_{i−1}) · `SEAICEstressFactor` (= 1, `seaice_readparms.F:510`) y u_1 la velocidad de la primera celda (`kSrf = 1`, `:54-59`). DWATN tiene unidades de kg m⁻² s⁻¹ (ρ₀ C_w |Δu|).
6. `EXTERNAL_FORCING_SURF` pasa ese `fu` al océano: `surfaceForcingU = fu · mass2rUnit` (`model/src/external_forcing_surf.F:222-224`).

## 2. El coeficiente `DWATN` y por qué vale cero en los dos modos

- `DWATN = 0` al inicio (`pkg/seaice/seaice_init_varia.F:78`). No se lee de la pickup (no está en `seaice_read_pickup.F`).
- Se calcula solo en `SEAICE_OCEANDRAG_COEFFS` (`pkg/seaice/seaice_oceandrag_coeffs.F:84-107`): DWATN = max(ρ₀ C_w |u_hielo − u_1|, `SEAICEdWatMin`), por la máscara de hielo, con `SEAICEdWatMin = 0.25` kg m⁻² s⁻¹ (`seaice_readparms.F:403`).
- Esa rutina solo se llama desde los resolvedores de la reología: `seaice_lsr.F:378`, `seaice_evp.F:738`, `seaice_krylov.F:201` y `seaice_calc_residual.F:79` (JFNK).

**Hielo fijo (`SEAICEuseDYNAMICS = .FALSE.`).** Se salta todo el bloque de `seaice_dynsolver.F:135-371` y no se llama ningún resolvedor. Entonces DWATN = 0 y u_hielo = 0 (`seaice_init_varia.F:73-74`). Con A = 1, fu = 0. **Bajo el hielo el océano no recibe nada**: ni el viento (que sí se quita, como pide el spec) ni el arrastre hielo–agua (que el spec pide y que no está). En el borde del hielo, los puntos U/V con A_u = ½ reciben la mitad del viento.

**Deriva libre (`SEAICEuseFREEDRIFT = .TRUE.`, con `SEAICEuseDYNAMICS = .TRUE.`, que es obligatorio porque la deriva libre está dentro del `IF` de `:135`).**
- `SEAICE_FREEDRIFT` (`pkg/seaice/seaice_freedrift.F:65-140`) resuelve el equilibrio estacionario τ_aire + τ_agua − m f k×u_hielo = 0. Usa `SEAICE_waterDrag·rhoConst` solo dentro de la rutina (`:98-120`) y **no escribe DWATN**.
- Con `SEAICEuseFREEDRIFT` se apagan LSR y EVP (`seaice_readparms.F:830, 856-857`), así que tampoco se llama a `SEAICE_OCEANDRAG_COEFFS`.
- Resultado: DWATN = 0 durante toda la corrida. El hielo se mueve, pero el océano bajo A = 1 recibe fu = 0. **El momento que el viento da al hielo se pierde**; no llega al océano.

## 3. KPP y GGL90: qué esfuerzo usan para u* (Step 1a)

- **KPP.** `KPP_CALC` llama a `KPP_FORCING_SURF` con `surfaceForcingU/V` (`pkg/kpp/kpp_calc.F:419-430`; comentario de `:377-379`). Ahí u* = (|τ|/ρ₀)^½, con τ promediado a los centros (`pkg/kpp/kpp_forcing_surf.F:186-207`). Si |τ| ≈ 0, u* = (½ ε Δz₁)^½ con ε = `phepsi` (`:198-204`).
  - Como `surfaceForcingU = fu·mass2rUnit` sale después de `SEAICE_MODEL` (§1, puntos 2 y 6), **KPP sí usa el esfuerzo hielo–agua cuando `pkg/seaice` lo calcula**. Con DWATN = 0 (los dos modos, §2), ese esfuerzo es 0 y KPP ve u* ≈ ε bajo el hielo: sin mezcla por esfuerzo superficial.
  - Lo que el spec (§3.6) pide, "bajo hielo, KPP recibe el esfuerzo hielo–agua como esfuerzo superficial", se cumple en el acoplamiento, pero el esfuerzo que llega es 0.
- **GGL90.** u*² se estima con `surfaceForcingU/V` promediados a los centros (`pkg/ggl90/ggl90_calc.F:787-792`, rama `calcMeanVertShear = .FALSE.`, la que se usa por omisión). La rama `.TRUE.` (`:764-782`) también usa `surfaceForcingU/V`. Es el mismo camino que KPP y la misma conclusión.
- **Convivencia con `seaice` y `obcs`.**
  - Ni `pkg/kpp` ni `pkg/ggl90` mencionan `seaice` ni `obcs` (grep vacío): leen solo `surfaceForcing*` y el estado.
  - OBCS se aplica antes (`OBCS_CALC`, `do_oceanic_phys.F:322`) y fija las velocidades en las fronteras.
  - `pkg/seaice` aplica sus propias condiciones de frontera (`OBCS_ADJUST_UVICE`, `seaice_model.F:190-195`; `OBCS_APPLY_SEAICE`, `:309-314`). Como no hay hielo en las fronteras abiertas (`ice.mask`, prueba `test_real_mask_inner_fjord_only`), no hacen nada.
  - KPP necesita OLx = OLy = 4 con `KPP_SMOOTH_SHSQ/DBLOC` (`pkg/kpp/kpp_calc.F:100-104`; `SIZE.h`).

## 4. Termodinámica apagada y otras banderas

- **Termodinámica:** `usePW79thermodynamics = .FALSE.` (por omisión `.TRUE.`, `seaice_readparms.F:268`).
  - Con ella apagada, `SEAICE_MODEL` no llama a `SEAICE_GROWTH`: solo calcula la carga `sIceLoad = h ρ_hielo` (`seaice_model.F:267-300`).
  - Tampoco hace falta `pkg/exf`: la comprobación que lo exige está dentro de `IF (usePW79thermodynamics)` (`seaice_check.F:686-691`; el comentario dice "Dynamics can be used without EXF").
  - El build `fjord_vmf` compila sin `exf`.
- **Carga del hielo:** `sIceLoad` solo entra en la presión de superficie con `useRealFreshWaterFlux = .TRUE.` (`model/src/external_forcing_surf.F:364-378`, con `ATMOSPHERIC_LOADING`), que aquí es `.FALSE.` (por omisión). El hielo no hunde la superficie. Si la hundiera, un escalón de 0.9 m en el borde del hielo forzaría un ajuste geostrófico espurio.
- **Viento sobre el hielo:** sin `exf`, τ_hielo = (`SEAICE_drag`/`OCEAN_drag`) · fu · máscara (`pkg/seaice/seaice_get_dynforcing.F:221-233`). Los dos valen 1×10⁻³ por omisión (`seaice_readparms.F:399-400`) y en `data.seaice` se ponen iguales de forma explícita. Así el hielo recibe el esfuerzo prescrito (F4).
  - Con `SEAICEscaleSurfStress = .TRUE.` (por omisión, `:266`), la fuerza sobre el hielo se escala por A (`seaice_dynsolver.F:266-282`).
- **Arrastre hielo–agua:** `SEAICE_waterDrag = 5.5×10⁻³` (PD0 §8; por omisión, `seaice_readparms.F:401`), escrito en `data.seaice`. Hoy solo lo usa la deriva libre para mover el hielo (§2).
- **Advección del hielo:** queda activa (`SEAICEadvHeff` y las demás, por omisión `.TRUE.`, `seaice_readparms.F:302`; se llama en `seaice_model.F:221-233`). Con hielo fijo u_hielo = 0 y no mueve nada; en F4 transporta h y A.
- **`tracForcingOutAB = 1`:** `pkg/seaice` lo exige (`seaice_check.F:82-86`). El valor por omisión es 0, porque `forcing_In_AB = .TRUE.` (`model/src/ini_parms.F:969, 1099-1102`).
  - Lo descubrió la primera corrida diagnóstica con hielo, que paró en `SEAICE_CHECK`.
  - Va en el `data` de todos los casos del fiordo, así que F2 y F3 solo difieren en el hielo (prueba `test_f2_f3_differ_only_by_the_ice`).
  - Sin flujos de calor ni de sal en superficie, el forzamiento de T y S es cero y el cambio no tiene efecto en los casos sin hielo.
  - Las corridas de diagnóstico de Δt de §2 de `pd1_config.md` se hicieron antes de añadirlo.
- **Pickup:** `pkg/seaice` escribe `pickup_seaice.<iter>`. Los casos con hielo lo exigen para la etapa libre (`pickup_stems` en `cases.yaml`; `run_campaign.py`).
- **Etapa libre:** sin archivo de viento, `fu` no se recarga (§1, punto 1). En agua abierta sigue en 0 (valor inicial tras la pickup). Bajo hielo, `SEAICE_OCEAN_STRESS` reescribe `fu` en cada paso. En los puntos de borde (A_u = ½), fu ← ½ fu_anterior + ½ τ_hielo, que converge a τ_hielo. Es inocuo mientras DWATN = 0.

## 5. Comprobación empírica (corridas cortas)

F2 y F4 no arrancan con las fronteras de Orlanski (§7). Para comprobar el esfuerzo, corrí los dos casos con `useOBCS = .FALSE.`: caja cerrada, **solo diagnóstico**, que no forma parte de la campaña. Datos:
- corridas `F2_noobcs_diag3h` y `F4_noobcs_diag3h` (3 h, Δt 60 s escalonado, `tracForcingOutAB = 1`, §4);
- `figures/pd1_timing.json`, campo `ice_stress`, de `analysis/pd1_timing.py`.

Se mide `oceTAUX`, el esfuerzo que recibe el océano después de `pkg/seaice` (`model/src/diags_oceanic_surf_flux.F:79-85`, que se llena en `do_oceanic_phys.F:1119`, después de `SEAICE_MODEL` en `:453`), en los puntos U con A = 1 a los dos lados. Se usa el A del mismo registro (`seaice`) si el hielo se mueve.

| Caso | Puntos U con A = 1 | máx \|τ_océano\| bajo hielo (N m⁻²) | máx \|τ_océano\| en agua abierta (N m⁻²) | Hielo |
|---|---|---|---|---|
| F2 (hielo fijo) | 13 072 | **0.0 en los 9 registros (0–3 h)** | 0.003 → 0.095 (pico del pulso) | u_hielo = v_hielo = 0; A total sin cambio (13 286) |
| F4 (deriva libre) | 13 072 al inicio → 2804 a las 2.7 h | **0.0 en los 9 registros** | igual | \|u_hielo\| hasta 0.19 m s⁻¹ y \|v_hielo\| hasta 0.21 m s⁻¹; A total baja de 13 286 a 12 395 (divergencia) |

**Confirmado:**
- Con A = 1, el océano recibe esfuerzo cero en los dos modos, aunque bajo el hielo fijo hay corrientes superficiales de hasta 2.0 cm s⁻¹ (percentil 99: 0.12 cm s⁻¹). Con arrastre, ρ₀ C_w \|u\|² sería ≤ 2×10⁻³ N m⁻², y con el piso lineal `SEAICEdWatMin`, hasta ~5×10⁻³ N m⁻².
- En F4 el hielo recibe el viento y se mueve a 0.2 m s⁻¹, pero no transmite nada al océano donde A = 1. Donde el hielo diverge (A < 1), el océano recibe (1 − A) veces el viento.

## 6. Opciones para el PD1 (punto PD1-S1; elegida la 2, §8; ver también §7)

1. **Aceptar el comportamiento de MITgcm y declararlo.**
   - F2: viento nulo y superficie sin fricción bajo el hielo, en lugar de una superficie con arrastre. El efecto sobre H4b va en contra de la hipótesis: sin arrastre, el hielo amortigua menos la energía interna del fiordo interior.
   - F4: el momento del viento se pierde en el hielo. F4 no tendría sentido físico y se recortaría (de todos modos es el primero del recorte).
   - Clasificación: (b) física ausente del modo elegido.
2. **[ELEGIDA, §8]** **Calcular `DWATN` también sin reología** con una copia modificada de `seaice_dynsolver.F` en `config/fjord_vmf/code/`. La fuente de MITgcm no se toca; `-mods` la sustituye.
   - Añadir, fuera del bloque `IF (SEAICEuseDYNAMICS ...)`, la llamada `SEAICE_OCEANDRAG_COEFFS(uIce, vIce, HEFFM, DWATN, ...)` antes de `SEAICE_OCEAN_STRESS`.
   - Con hielo fijo: τ = ρ₀ C_w |u_1| u_1 hacia el reposo, con el piso lineal `SEAICEdWatMin` (0.25 kg m⁻² s⁻¹, equivalente a |Δu| = 4.4 cm s⁻¹ con C_w = 5.5×10⁻³).
   - Con deriva libre: el arrastre consistente con el hielo que se mueve.
   - Es una línea de código, con una prueba de reinicio y una corrida corta. Riesgo: código propio en el build, que se declara.
3. **Hielo "casi fijo" con la reología:** `SEAICEuseDYNAMICS = .TRUE.` con LSR y una resistencia del hielo muy alta (`SEAICE_strength`). LSR calcula DWATN.
   - Solo vale para F2, y el hielo no queda exactamente inmóvil.
   - Cuesta un resolvedor LSR por paso (no medido).
4. **Solo para F4, solo con el namelist: `useHB87StressCoupling = .TRUE.`** (por omisión `.FALSE.`, `seaice_readparms.F:264`).
   - Con esa bandera, `SEAICE_OCEAN_STRESS` usa fu ← (1 − A) fu + A·τ_viento-sobre-hielo + divergencia del esfuerzo interno (`seaice_ocean_stress.F:64-97`).
   - En deriva libre la divergencia es cero (`seaice_dynsolver.F:315-316`), así que bajo A = 1 el océano recibe el viento sobre el hielo (τ_hielo = fu, §4).
   - El equilibrio de deriva libre da, como esfuerzo hielo–agua, τ_aire − m f k×u_hielo (más la inclinación de la superficie; `seaice_freedrift.F:82-84`). HB87 pasa τ_aire sin el término de Coriolis del hielo: con m = 910 kg m⁻² y \|u_hielo\| = 0.2 m s⁻¹ (§5), m f \|u_hielo\| ≈ 0.03 N m⁻², frente a un pico de 0.1 N m⁻².
   - No necesita código, pero **no sirve para F2**: sin dinámica, el término del viento pasaría entero a través del hielo fijo, lo contrario de lo que se quiere.
5. **Sin `pkg/seaice`:** tapa con arrastre cuadrático (por ejemplo, `pkg/shelfice` con un calado mínimo).
   - Es un cambio grande de configuración, choca con la celda superficial de 1.5 m y no lo recomiendo.

**Inercia del hielo en F4.** La deriva libre de MITgcm resuelve un equilibrio estacionario, sin du_hielo/dt: el lado derecho de `seaice_freedrift.F:82-84` solo tiene viento y Coriolis, y la solución es algebraica, `:96-135`. El spec (§3.5) pedía "losa con inercia y arrastre".
- Es **(b) física ausente** del modo elegido. Se declara.
- Escala: el tiempo de ajuste del hielo es m/(ρ₀ C_w \|Δu\|) ≈ 910/(1027 × 5.5×10⁻³ × 0.2) ≈ 800 s, menos que la escala del pulso (σ = 1 h). La aproximación es razonable para la respuesta forzada, pero no representa las oscilaciones inerciales del hielo.

**Recomendación técnica (la decisión es de Esteban):**
- **F2:** la opción 2. Es la que reproduce lo que el spec y el PD0 (§8, "τ = ρ₀ C_w |u_hielo − u|(u_hielo − u) con u_hielo = 0") daban por hecho, con el C_w ya aprobado.
- **F4:** la opción 2 también lo resuelve, con el arrastre coherente con el hielo que se mueve. Si se prefiere no tocar código, la opción 4 (HB87) es aceptable para F4, declarando el término de Coriolis que falta y la falta de inercia.
- En cualquier caso hace falta resolver antes PD1-S2 (§7).

## 7. `pkg/seaice` con `pkg/obcs` de Orlanski: F2 y F4 no arrancan (punto PD1-S2; elegida la 1, §8)

- **Qué pasó.** La corrida corta de F2 (`F2_smoke3h`, Δt 60 s, escalonado) y la de F4 terminan en la inicialización con `STOP ABNORMAL END: S/R OBCS_CHECK`. `STDERR.0000` dice:
  - "useOrlanski* OBC not yet implemented for SEAICE" (`pkg/obcs/obcs_check.F:386-406`; bloque `IF (useSEAICE)` en `:400-405`);
  - "for SEAICE OBCS, use pkg/exf to read input files" (`:571-579`, `IF (useSEAICE .AND. .NOT. useEXF)`).
- **No tiene salida por namelist.** Las dos comprobaciones están fuera de cualquier opción CPP. `DISABLE_SEAICE_OBCS` (`SEAICE_OPTIONS.h:266`) solo quita las llamadas de `seaice_model.F:190-195` y `:309-314`, no las comprobaciones. El hielo empieza a 29 km (146 celdas) de las celdas de frontera abierta, pero `obcs_check` no lo mira.
- **Lo que no hice:** cambiar las fronteras de Orlanski en los casos con hielo. F2 y F3 dejarían de ser comparables (H4b).

**Opciones (no elegidas):**
1. **[ELEGIDA, §8]** **Copia modificada de `obcs_check.F` en `config/fjord_vmf/code`.** Se relajan las dos comprobaciones cuando `DISABLE_SEAICE_OBCS` está definido, y se define en `SEAICE_OPTIONS.h`. La fuente de MITgcm no se toca.
   - Justificación: con `DISABLE_SEAICE_OBCS` no se aplica ninguna condición de frontera al hielo, y el hielo está lejos de las fronteras.
   - Riesgo: código propio, que se declara.
   - Se combina bien con la opción 2 de §6 (otra copia, de `seaice_dynsolver.F`).
2. **Compilar `pkg/exf` y llevar el viento por `exf`** (`ustressfile`/`vstressfile`) en **todas** las corridas del fiordo, para que F2 y F3 sigan siendo comparables.
   - Resuelve el segundo error, pero no el de Orlanski: también habría que cambiar el tipo de frontera.
   - Es el cambio más grande y el que más se aleja de la tesis (Orlanski, viento por `external_fields_load`).
3. **Fronteras prescritas o solo esponja** (sin Orlanski) en todas las corridas del fiordo. También necesita `exf` (segundo error) y cambia la radiación por las fronteras abiertas, que es parte del diseño de la tesis.
4. **Sin `pkg/seaice` para F2.** Como con `SEAICEuseDYNAMICS = .FALSE.` el océano recibe esfuerzo cero bajo el hielo (§2 y §5), el F2 que MITgcm haría es exactamente "viento × (1 − A) en los puntos U/V" sin ninguna otra física.
   - Se puede generar con el viento enmascarado y sin código nuevo.
   - Pero da la tapa sin fricción, no el arrastre del spec.
   - F4 no tiene equivalente sin `pkg/seaice`.

**Recomendación técnica (la decisión es de Esteban):** opción 1 de esta sección más la opción 2 de §6. Son dos copias pequeñas de archivos de MITgcm en `code/`, con una prueba de reinicio y una corrida corta. Darían el F2 del spec (arrastre bajo el hielo fijo) y un F4 coherente, con las fronteras de Orlanski en todas las corridas.

## 8. Implementación del PD1 y verificación (2026-11-12)

**Copias modificadas en `config/fjord_vmf/code`** (también en las de 400 m, 100 m y la subcaja; una prueba exige que sean idénticas). La fuente de MITgcm no se toca.
- Cada copia empieza con una cabecera: archivo original, versión (checkpoint69j, git 8e1c4076) y sha256 del original. Los cambios están entre líneas `C--- Exp. 4 BEGIN` / `C--- Exp. 4 END`.
- `tests/test_mitgcm_mods.py` comprueba que, quitando la cabecera y esos bloques, cada copia es el original byte a byte, y que el sha256 de la cabecera es el del original.

| Copia | sha256 del original | Cambio |
|---|---|---|
| `obcs_check.F` (`pkg/obcs`) | b19d4e2e…6bcc1 | incluye `SEAICE_OPTIONS.h`; con `DISABLE_SEAICE_OBCS` se saltan las dos paradas (líneas originales 400-405 y 571-579) |
| `seaice_dynsolver.F` (`pkg/seaice`) | 3730d29b…5a660 | antes de `SEAICE_OCEAN_STRESS`, si no corrió ningún resolvedor de reología (`.NOT.SEAICEuseDYNAMICS .OR. SEAICEuseFREEDRIFT`), llama a `SEAICE_OCEANDRAG_COEFFS(uIce, vIce, HEFFM, DWATN, …)`, como `seaice_lsr.F:378` |
| `obcs_add_tides.F` (`pkg/obcs`, V3) | bb3ca468…245ea8 | con `EXP4_TIDE_RAMP_S` (en `OBCS_OPTIONS.h`), rampa de medio coseno de la velocidad de marea durante un ciclo M2 |

Además, `SEAICE_OPTIONS.h` define `DISABLE_SEAICE_OBCS`.
- Eso quita las llamadas de frontera de `seaice_model.F:190-195` y `:309-314` (`OBCS_ADJUST_UVICE`, `OBCS_APPLY_SEAICE`).
- **No quita** la de `seaice_dynsolver.F:324-328`. En deriva libre (F4) sigue llamando a `OBCS_APPLY_UVICE`, que pone la velocidad del hielo en la frontera igual a `OB*uice`/`OB*vice`, que valen 0 (`pkg/obcs/obcs_init_variables.F:221, 241, 283`).
- Es inocuo mientras el hielo no llegue a la frontera: empieza a 29 km (146 celdas) de las celdas de frontera abierta (`analysis/pd1_timing.ice_distance_to_open_boundaries`).
- **Comprobación obligatoria después de cada corrida con hielo (Task 6):** como con `DISABLE_SEAICE_OBCS` el contacto sería silencioso, el análisis tiene que verificar que en las 280 h el hielo (SIarea > 0) se queda a ≥ 10 celdas (el ancho de la esponja) de toda celda de frontera abierta (`open_boundary_distance.passes` en `figures/pd1_timing.json`). Si no se cumple, la corrida se declara afectada.
- En las corridas cortas de 3 h: F2 y F4 a 146 celdas, la pasan.
- **Implementación (Task 6, parte A):** `analysis/ice_check.py` (CLI `analysis/validate.py ice-check RUN`). Recorre los registros `seaice` de **todas las etapas** (la versión de la Task 4, `pd1_timing.ice_distance_to_open_boundaries`, solo leía el directorio principal, es decir, la etapa forzada de 0–60 h) y sale con código 1 si falla. El resultado se guarda en `figures/ice_ob_check.json` con el sha256 del estado final de la corrida, y `load.run` **rechaza** cualquier corrida con hielo sin una comprobación vigente que pase. Así ningún análisis (V2–V4, H4) puede usar F2 o F4 sin ella. Sobre las salidas reales se corre en la parte B.

**Verificación** (`figures/pd1_timing.json`, `ice_stress.drag_check`; `figures/pd1_restart.json`):
- **F2 y F4 arrancan con Orlanski** (`F2_pd1final_3h`, `F4_pd1final_3h`: 3 h, Δt 120 s escalonado, OK).
- **Esfuerzo bajo A = 1 frente a la fórmula.** `analysis/pd1_timing.drag_check` compara oceTAUX con τ = ½(DWATN_i + DWATN_i−1)(u_hielo − u₁), con DWATN = max(ρ₀ C_w |Δu|, `SEAICEdWatMin` = 0.25 kg m⁻² s⁻¹), en cada registro:
  - **F2:** 13 072 puntos U; error máximo **0.0** en los 8 registros. Todos los puntos están en el piso lineal (|Δu| ≤ 1.7 cm s⁻¹ < 0.25/(ρ₀C_w) = 4.4 cm s⁻¹). |τ| crece hasta 4.4×10⁻³ N m⁻² a las 2.7 h. Antes del cambio era exactamente 0 (§5).
  - **F4:** hielo a 0.22 m s⁻¹; |τ| bajo A = 1 de hasta 0.097 N m⁻², en régimen cuadrático (fracción en el piso < 2 %).
    - La velocidad del hielo del paso comprobado no se guarda: `SIuice/SIvice` se escriben al empezar el paso (`model/src/do_statevars_diags.F:96`), antes de la deriva libre de ese paso.
    - Por eso la comprobación **recalcula la deriva libre del paso** con el álgebra de `seaice_freedrift.F:65-170` (`analysis/pd1_timing.free_drift_velocity`), a partir del estado diagnosticado: u₁, η (inclinación), h y A del registro, y el viento interpolado a t = it·Δt.
    - **En el interior del hielo** (≥ 2 celdas del borde del hielo y de la costa; 2630–3274 puntos U): error máximo ≤ 1.6×10⁻³ N m⁻², **≤ 1.6 % del máximo de |τ| del registro durante el pulso y 2.7 % cuando se apaga** (registro de las 2.7 h), con mediana ≤ 9×10⁻⁵ N m⁻² (≤ 0.2 %).
    - En el borde del hielo y junto a la costa, la reconstrucción de las máscaras (`SIMaskU`, A_u, h en el borde) no es exacta, y los errores llegan a 0.5–0.7 del máximo en unos pocos puntos (0.51–0.67 en los 8 registros, `max_rel_err_to_max_tau` de `figures/pd1_timing.json`; antes decía "~0.5", corregido en la Task 6). No es una discrepancia del modelo.
    - Con `SIuice` del inicio del paso, en cambio, el error es casi uniforme (~3.4×10⁻³ N m⁻², 5 % del máximo): el desfase de un paso.
- **Reinicio** (`scripts/restart_test.py`): 20 pasos seguidos frente a 10 + pickup + 10. **Idéntico bit a bit** en F2 y F4: pickups del modelo, de `seaice` y de Orlanski del paso 20, y diagnósticos `state3d`, `surf` y `seaice` del paso 10.
- Pendiente para declarar: la deriva libre sigue sin inercia del hielo, (b) (§6).
