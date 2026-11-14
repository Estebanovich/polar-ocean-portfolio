# PD3 — Hipótesis H4a–H4d y Figuras 2–4 (Task 7)

Fecha: 2026-11-13. Material para que Esteban decida en el PD3.

**Archivos:**
- veredictos: `figures/hypotheses.json`, de `analysis/hypotheses.py`, que solo lee los JSON;
- datos: `figures/partition.json` (bandas, estaciones, energía de la caja, trabajo del viento) y `figures/budget.json` (balance de Kang y Fringer de cada corrida), de `analysis/partition.py`;
- figuras: `figures/fig2_spectra_partition`, `fig3_strait_flux_budget` y `fig4_ice_maps` (`.png`, `.pdf` y `.json` con sus números), de `analysis/figures.py`.

**Ventanas reales** (`partition.json` → `runs.*.record`; `budget.json` → `runs.*.window_h`): las corridas del fiordo cubren **0–279.67 h**, 840 registros cada 1200 s (el estado tras el último paso no es un registro de diagnóstico); R1 y R1k en 0–120 h, 481 registros cada 900 s; R2 y R3 en **0–280 h**, 1121 registros cada 900 s. "0–T_fin" significa eso en todo el texto.

**Corridas:** las 13 de la campaña (F1–F5, F1c, F1g, V3, V4, R1, R1k, R2, R3), todas con `require_ok`. F2 y F4 tienen la comprobación del hielo. Los umbrales son los del PD0 §5, sin cambios (`validate.THRESHOLDS`, comprobados contra el texto del PD0).

**Definición del flujo:** decide el trabajo de presión de Kang y Fringer, F_bc = Σ u_n′ p′ dz (Enmienda 2 para H4d; **Enmienda 3**, 2026-11-13, para H4a, H4b(ii) y H4c(ii)). En H4b(ii) se compara con el trabajo de presión barotrópico F_bt. La Enmienda 3 se escribió y se versionó (commit 1c1f45c) antes de calcular ningún número de H4. F + A se reporta como descriptivo.

Signo: "+ hacia fuera" en `*_out_J` y `F_bc_pressure_out`; "+ hacia dentro" en `budget.json` (`totals_in_J`). 1 GJ = 10⁹ J.

## Resumen

| Hipótesis | Criterio (PD0 §5) | Resultado | Desenlace |
|---|---|---|---|
| **H4a** (F1) | banda baja > 25 % → no concluyente; si no, φ_NI > 0.5 apoyada, φ_K > 0.5 alternativa | banda baja = **0.93** del total en valor absoluto | **no concluyente con T_fin** |
| **H4b(i)** (F2/F3) | Ē_in(F2)/Ē_in(F3) ≤ 0.5 | 0.578 / 2.541 GJ = **0.228** | **se cumple** |
| **H4b(ii)** (F2) | F_bc,in ≥ \|F_bt,in\| | 0.0067 GJ frente a 0.0206 GJ | **no se cumple** |
| **H4c** (F1/F3) | apoyada si (i) y (ii) ≥ 1.5; refutada si las dos ≤ 1/1.5 | (i) W_NI: **0.83**; (ii) 119.6, con las dos sumas resueltas < 0 (bandera) | **no distinguible** |
| **H4d** R1k/R1 (0–120 h) | η(KPP)/η(tesis) ≥ 1.5 | η = 0.32 % frente a 0.082 %: **3.93** | **apoyada** |
| **H4d** F1/F1c (0–T_fin) | ídem | η(F1) = −6.07 % (bandera: numerador < 0) frente a 0.19 %: **−32.7** | **refutada** (aritmética literal, con bandera) |

La regla del factor 1.5 se aplica tal como está congelada. Su respaldo era la cota de resolución del 25 % de V4a, y **V4 falló** (V4a y V4b, `figures/validation.json`): esa cota no está demostrada.

## H4a — partición del flujo radiado (F1)

`hypotheses.json` → `H4a`; `partition.json` → `runs.F1.straits.total.F_bc_out_J`.

- Flujo neto por los dos estrechos, por banda (GJ, + hacia fuera): baja −1.144; Kelvin +0.018; cercana a la inercial −0.043; superinercial +0.010; seiche −0.015.
- Fracción de la banda baja: |F_low| / Σ₅|F_b| = **0.930 > 0.25**. El desenlace es "no concluyente con T_fin", por la regla del PD0, que se evalúa primero.
- Las φ_b se reportan igual: Kelvin −0.60, cercana a la inercial +1.45, superinercial −0.34, seiche +0.49. La suma de las cuatro bandas resueltas es **negativa** (−0.030 GJ: entra al fiordo), y por eso las fracciones salen fuera de [0, 1]. La bandera `resolved_net_inward` está activa; la regla no cambia.
- Descriptivo, F + A: banda baja 0.93 del total, resueltas −0.018 GJ.
- Es lo que anticipaba el PD2: el neto lo domina la banda baja (> 140 h). Esa banda comparte la deriva de la referencia fija (`hypotheses.json` → `general.low_band_caveat`), así que no se puede leer como intercambio físico sin la comprobación *post hoc*.

## H4b — hielo fijo (F2 frente a F3)

`hypotheses.json` → `H4b`; `partition.json` → `runs.F2/F3.energy.totals`; `budget.json` → `runs.F2.totals_in_J`.

**(i) se cumple.** Ē_in sin la banda baja (KE_bc + APE exacta, campos filtrados con paso alto en 2Δf):
- F2: 0.578 GJ (KE_bc 0.0006 y APE 0.578);
- F3: 2.541 GJ (KE_bc 1.239 y APE 1.302);
- cociente **0.228 ≤ 0.5**.

**(ii) no se cumple.** Trabajo de presión neto hacia dentro por los dos estrechos, 0–279.67 h:
- baroclínico F_bc,in = 0.0067 GJ;
- barotrópico F_bt,in = 0.0206 GJ;
- 0.0067 < 0.0206. Descriptivo, F + A: 0.0078 frente a 0.0189 GJ, mismo resultado.

**Residuo del balance** (pedido por el PD0): D_bc relativo a las fuentes brutas es −52 en F2 y +0.78 en F3 (`budget_residual_D_bc_rel`). En F2 el residuo es enorme porque casi no hay fuentes (P ≈ 0, F ≈ 0.008 GJ, C ≈ −0.014 GJ) y aun así E_bc crece 1.16 GJ (ver "Deriva de la APE").

**Deriva de la APE exacta (hallazgo; afecta a la lectura de H4b(i)).**
- En F2, sin viento sobre el fiordo, la APE exacta de la caja crece de forma monótona de 0 a 1.15 GJ en 280 h, mientras KE_bc termina en 0.014 GJ (`budget.json`, `runs.F2.energy_J.ke_bc.end`).
- El perfil medio de la caja se suaviza durante la corrida, igual en F2 que en F3 (`figures/posthoc_moving_reference.json` → `profile_drift`, *post hoc*). En F2, la densidad del primer nivel pasa de 1027.756 a 1027.787 kg m⁻³, y el escalón neutro de los niveles 4–7 se convierte en un gradiente. El cambio máximo es 0.031 kg m⁻³ en F2, 0.039 en F3, 0.22 en F1 y R2, y 0.003–0.006 en R1 y R1k. `pkg/seaice` corre sin termodinámica, así que no hay salmuera.
- Contra el ρ_b inicial fijo que manda la Enmienda 1, esa mezcla aparece como APE que crece. La de F1 llega a 12.7 GJ al final, y la de R2 a 149 GJ.
- Clase **(c)**: es la elección del método (referencia fija en el perfil inicial) frente a un perfil de fondo que se difunde. El modelo no está en cuestión.
- **Sensibilidad *post hoc*** (escrita después de ver la deriva; no es el criterio): el cociente de KE_bc sin la banda baja es 0.0006 / 1.239 = **0.0005** (`H4b.i.post_hoc_ke_bc_only`). Esa es la única cifra libre de la deriva. La APE exacta no es lineal: la del campo filtrado puede superar a la total (F2: 0.578 frente a 0.538 GJ, `partition.json` → `runs.F2.energy.totals`), así que la deriva no se puede restar del Ē_in congelado. El desenlace de (i) se mantiene con las dos cifras (0.228 y 0.0005).

## H4c — estación (F1 frente a F3)

`hypotheses.json` → `H4c`.

- **(i)** W_NI = ∫∫ τ·u_NI dA dt en la caja: F1 6.05 GJ, F3 7.33 GJ, cociente **0.83**. No llega a ninguna de las dos colas.
- **(ii)** Trabajo de presión saliente, suma de las cuatro bandas resueltas: F1 −0.0297 GJ, F3 −0.00025 GJ. Cociente 119.6, pero **los dos son negativos** (entran): banderas `f_res_f1 <= 0` y `f_res_f3 <= 0`. Descriptivo, F + A: −0.0187 y +0.0018 GJ.
- **No distinguible.** Con (i) = 0.83 el desenlace es "no distinguible" se lea como se lea (ii): no hay reinterpretación que lo cambie.
- **Advertencia del PD0 §5:** el modo 1 de invierno está mal resuelto a 200 m (R₁/Δx = 1.7–6.6).

## H4d — cierre (KPP frente al cierre de la tesis)

`hypotheses.json` → `H4d`; `budget.json` → `runs.R1`, `runs.R1k_0_120h`, `runs.F1`, `runs.F1c`.

**R1k frente a R1 (0–120 h, boca de la bahía): apoyada.**
- R1k: η = 0.333/103.5 GJ = **0.32 %**. R1: η = 0.284/346.3 GJ = **0.082 %**. Cociente **3.93**.
- R1 reproduce el 0.08 % de la tesis (bitácora del 11 de noviembre) y los 0.28 GJ del contraste de la Task 5.
- **Lectura descriptiva:** el cambio viene sobre todo del denominador. Con KPP el viento hace 3.3 veces menos trabajo (103.5 frente a 346.3 GJ), porque el esfuerzo se reparte en la capa de mezcla y la velocidad superficial es menor. El numerador solo crece ×1.17.
- η′ = F/P_bc: 0.325 % frente a 0.082 %. F + A saliente: +0.268 GJ en R1k frente a −5.40 GJ en R1, por la advección que entra en R1 (Task 5).

**F1 frente a F1c (0–279.67 h, estrechos): refutada, por aritmética literal, con bandera** (decisión de Esteban, 2026-11-13: se mantiene así, con la bandera).
- η(F1) = −1.171/19.30 GJ = **−6.07 %**: el trabajo de presión baroclínico neto **entra** al fiordo, y casi todo es de la banda baja (−1.144 GJ).
- η(F1c) = +0.0679/36.55 GJ = **0.19 %**.
- Cociente **−32.7 ≤ 1/1.5**: "refutada" por la regla aplicada al número tal como sale. La bandera `kpp: numerator F_bc_pressure_out <= 0` está en el JSON. Los rulings piden reportar y marcar, sin reinterpretar. **La refutación viene del signo del término de intercambio neto de baja frecuencia, no de que KPP radie menos.** La banda baja del trabajo de presión por los estrechos (−1.144 GJ en F1, −0.0097 GJ en F1c) comparte la deriva de la referencia fija (ver "Comprobación *post hoc*"), así que no se puede atribuir a un intercambio físico sin esa comprobación.
- *Post hoc* (decisión de Esteban; nunca un veredicto; `hypotheses.json` → `post_hoc_h4d_resolved_bands`): η con solo las cuatro bandas resueltas (referencia fija). F1: −0.0297 GJ / 19.30 GJ = **−0.154 %**. F1c: +0.0775 / 36.55 = **+0.212 %**. Cociente **−0.72**. Con la referencia fija, también la suma resuelta de F1 entra. Con la referencia móvil cambia de signo (ver la comprobación).
- Descriptivo: η′(F1) = −6.37 % y η′(F1c) = 0.19 %. F + A saliente: −1.19 GJ (F1) frente a +0.35 GJ (F1c).

## Descriptivos *post hoc* (no son veredictos; `hypotheses.json` → `post_hoc`)

### Flujo por estrecho y por banda (R24)

Trabajo de presión, GJ, + hacia fuera; cuatro bandas resueltas (K, NI, SI, S) y banda baja aparte (`post_hoc.band_flux_per_strait`):

| Corrida | Estrecho | K | NI | SI | S | baja |
|---|---|---|---|---|---|---|
| F1 | Akselsundet | +0.052 | −0.034 | +0.011 | −0.009 | −0.997 |
| | Mariasundet | −0.034 | −0.009 | −0.001 | −0.005 | −0.147 |
| F1c | Akselsundet | +0.159 | −0.027 | +0.005 | −0.014 | +0.102 |
| | Mariasundet | −0.033 | −0.011 | +0.001 | −0.003 | −0.112 |
| F1g | Akselsundet | −0.178 | −0.013 | +0.006 | −0.013 | −1.138 |
| | Mariasundet | −0.138 | −0.014 | 0.000 | −0.007 | −0.985 |
| F5 | Akselsundet | −0.048 | −0.019 | +0.014 | −0.002 | −0.709 |
| | Mariasundet | +0.010 | −0.005 | −0.001 | −0.002 | −0.712 |
| F3 | ambos | 0.000 | −0.0002 | 0.000 | 0.000 | +0.022 |
| F2 | ambos | −0.0005 | −0.0002 | 0.000 | 0.000 | −0.006 |
| F4 | ambos | −0.0004 | −0.0002 | −0.0001 | 0.000 | +0.035 |
| V4 (400 m) | ambos | +0.151 | −0.013 | +0.013 | −0.001 | +0.172 |

- La banda baja domina el neto en todas las corridas de verano con KPP o GGL90 (con la salvedad de la deriva de la referencia). F1c es la excepción: con el cierre de la tesis, el neto lo lleva la banda de Kelvin hacia fuera.
- Con hielo (F2), el flujo por los estrechos es dos órdenes de magnitud menor que en verano.

### Comparación al estilo de Arneborg y Liljebladh (2001)

En Gullmar, A&L filtraron en **0.01–0.07 cph** (p. 2561, biblioteca de la tesis, l. 601–605) y encontraron que la energía baroclínica neta acumulada **entra** casi siempre por la boca (l. 631–636). Aquí, con el trabajo de presión en esa misma banda y + hacia dentro (`post_hoc.arneborg_liljebladh_net_in_J`):
- entra: F1 +0.038 GJ, F2 +0.0004, F3 +0.0003, F4 +0.0004, F5 +0.015, F1g +0.201, V3 (solo marea) +0.460;
- sale: F1c −0.021, V4 (400 m) −0.076, R1 −0.112, R1k −0.531 (0–240 h; −0.144 en 0–120 h), R2 −0.204.

El signo coincide con Gullmar en las corridas del fiordo a 200 m con KPP o GGL90, en las de invierno (con valores minúsculos) y en V3. No coincide con el cierre de la tesis (F1c), a 400 m (V4) ni en la bahía. Con la resolución sin demostrar (V4 cambia el signo), no es una conclusión robusta. Es una comparación de signo en un fiordo distinto, sin umbral.

### Limitación de los estrechos (Johnson y Garrett 2006)

`post_hoc.strait_limitation`:
- Mariasundet: 0.50 km (el más ancho de sus dos pasos, NPI), frente a R₁ = 3.62 / 2.47 km en verano (H máx. / H medio). W/R₁ = 0.14–0.20 en verano y 0.38–1.46 en invierno.
- Akselsundet: 1.07 km. W/R₁ = 0.30–0.44 en verano.
- Johnson y Garrett (2006), citado por Jackson et al. (2018, `literatura/Jackson_2018 … .md`, l. 704–709): salvo que el estrecho sea "considerably narrower than the deformation radius", la mayor parte de la onda de Kelvin incidente se transmite.
- En verano, Mariasundet es entre 5 y 7 veces más estrecho que R₁, así que cabe esperar una transmisión reducida de Kelvin por ese estrecho. Además, la malla lo abre a 1.33 km de agua (Task 1b, clase (a) de la métrica).
- J&G 2006 no está en las bibliotecas locales [verificar] (`readinglist.md`).

### Puente (descriptivo, spec §7; Fig. 2)

Fracción con signo del trabajo de presión saliente, F_b / Σ₅|F_b| (`post_hoc.bridge`):

| | banda baja | cercana a la inercial |
|---|---|---|
| R1 (tesis, 120 h) | +0.49 | +0.36 (banda **no resoluble** a 120 h: 1.7 bins, PD0 §3) |
| R1k (240 h) | +0.24 | +0.61 |
| R2 (bahía a 78°N) | +0.38 | +0.53 |
| F1 (fiordo) | −0.93 | −0.035 |

- En la bahía, a f de la tesis o a 78°N, la banda cercana a la inercial lleva la mayor parte del flujo resuelto, y el neto sale.
- En el fiordo real, el neto lo domina la banda baja y entra.
- La energía cinética superficial en las estaciones va en la misma dirección. En el interior de R2 es 0.995 cercana a la inercial; en la costa sur de F1 es 0.56 de Kelvin y 0.09 cercana a la inercial.
- R3 (sin bahía) no tiene caja ni flujo; su estación de plataforma es 0.999 cercana a la inercial.
- Como predecía el PD0 §2.3, punto 3, el cambio respecto a la tesis lo pone la geometría, no f: R2, con f de 78°N y el perfil de verano, sigue siendo cercano a la inercial.

### Fracción de energía atrapada (PD0 §3, secundario)

KE_bc integrada en la vertical, en celdas a menos de 2R₁ de tierra (`post_hoc.trapped_fraction`):
- **F1:** con R₁(H máx.) = 3.6 km, 2R₁ = 7.2 km abarca casi todo el fiordo (W = 8.4 km), y la fracción es 1.00 en todas las bandas. No es informativa. Con R₁(H medio) = 2.5 km: Kelvin 0.96, cercana a la inercial 0.84.
- **F2** (R₁ de invierno, 1.32 km): Kelvin 0.59, cercana a la inercial 0.69.
- **F3:** Kelvin 0.62, cercana a la inercial 0.34.

## Métodos (fijados antes de calcular cualquier número de H4)

- **Flujos:** las piezas de `budget.kf_budget`, por cara y por nivel (`partition.face_series`); reproducen a `kf_budget` con un error relativo máximo de 8×10⁻¹⁶ en todas las corridas (`check_face_series_vs_kf_budget_max_rel`). Partición por bandas: el reparto coespectral exacto de la Task 6 (`validate.band_flux`), sin ventana. Las bandas suman −Δt Σ F.
- **Bandas** (PD0 §3, `theory.json`):
  - el fiordo, V3 y V4 con 280 h (el registro cubre 0–279.67 h);
  - R2 y R3 con f de 78°N, 280 h y el límite de seiche de la bahía (0.149 cph);
  - R1k con las de la bahía de la tesis a 240 h;
  - R1 y R1k en 0–120 h con la misma fórmula a 120 h, donde la banda cercana a la inercial tiene 1.7 bins y **no es resoluble** (marcado en `resolvable`).
- **H4b(i):** paso alto en 2Δf con FFT sin ventana, punto a punto, con los mismos bordes de banda. Una ventana Hann sesgaría la media temporal de un cuadrado. "La misma maquinaria Hann/FFT" se leyó como los mismos bordes y la FFT.
- **W_NI:** puntos U/V de la caja con los pesos de `kf_budget`. Sobre todas las bandas, W es igual a ∫(P_bt + P_bc) dt de `kf_budget` (`wind_work_J.check_all_vs_kf_P_J`).
- **Esfuerzo:** en agua abierta, el de los archivos de forzamiento (`load.forcing_stress`); en F2 y F4, `oceTAUX/oceTAUY`, que es lo que recibe el océano bajo el hielo.
- **Estaciones:** `cases.yaml` no define estaciones, así que es una elección de método. En el fiordo: la de la plataforma de V2a y las de la costa sur de V2b (reglas de la Task 6). En la bahía: las de la tesis (`v1_thesis.STATIONS`). Espectro de energía cinética superficial con una ventana Hann.
- **R1k frente a R1:** la serie de R1k se corta en t = 120 h. R3 no tiene caja: la bahía es tierra.
- **Lectura:** cada corrida se leyó de LaCie una sola vez, a un almacén local de la caja con 2 celdas de halo (`build_store`). Los núcleos dan resultados idénticos bit a bit a la lectura directa (corrida corta de 48 h y sintético).

## Limitaciones

- **V2a falló** (c): la causa del desplazamiento al azul (mediana 1.09 f) no está identificada. Las bandas clasifican frecuencias, no mecanismos.
- **V2b falló** (a): la métrica no está definida con viento uniforme. V3 mide la onda de Kelvin del modelo, 0.48 m s⁻¹.
- **V4 falló:** el neto por los estrechos no converge (+0.32 GJ a 400 m frente a −1.17 GJ a 200 m). Los cocientes de H4c y H4d no tienen respaldo de resolución.
- **Invierno a 200 m:** R₁/Δx = 1.7–6.6; afecta a H4b y H4c.
- **Deriva de la APE exacta** con ρ_b inicial fijo (c). Domina dE_bc y el residuo D en todas las corridas largas (F1 12.7 GJ, R2 149 GJ de APE final). Los residuos D de `budget.json` no miden disipación.
- **Transmisión por estrechos estrechos** (J&G 2006), y Mariasundet abierto de más en la malla.
- **T_fin = 280 h:** la banda de periodo > 140 h no se resuelve y domina el neto por los estrechos (H4a no concluyente; H4d de F1 negativo). Esa banda comparte la deriva de la referencia fija; ver la comprobación *post hoc*.
- Un fiordo y un evento de viento, sin marea en la física de estudio, y hielo sin termodinámica.

## Qué necesita Esteban decidir

1. ~~H4d F1/F1c~~ **Decidido** (Esteban, 2026-11-13): se mantiene "refutada" con la bandera, más la explicación y el η *post hoc* de las bandas resueltas.
2. **Deriva de la APE:** ¿se reporta solo como limitación (c), o se pide un diagnóstico aparte, por ejemplo con un perfil de referencia que se actualice? Sería *post hoc*.
3. **H4b(i):** se cumple con o sin la deriva. ¿Basta con la sensibilidad de KE_bc ya reportada?

## Comprobación *post hoc*: trabajo de presión con referencia móvil (ronda de correcciones 1)

**Por qué.** p′ se construye con p̃ = ρ₀gη + g∫ρ′dz, y ρ′ se toma respecto al mismo ρ_b inicial fijo cuya deriva infla la APE. Por eso la deriva del fondo también entra en F_bc, con un término ∫p′_deriva(z) Q′(z) dz allí donde hay un intercambio baroclínico medio Q′. Ese término es de baja frecuencia: cae en la banda que hace a H4a "no concluyente" y fija el signo de η(F1).

**Método, fijado aquí antes de calcular** (escrito y versionado antes de la primera corrida de la comprobación):
1. **Referencia móvil.** En cada registro y en cada nivel k, ρ̄_caja(z_k, t) es la media de ρ en las celdas mojadas de la caja, ponderada por volumen (Δz·rA). La caja es el fiordo interior, o la bahía para R1 y R1k. En los niveles sin celdas mojadas de la caja se usa el ρ_b fijo; ninguna cara de los estrechos llega a esos niveles.
2. ρ′_m = ρ − ρ̄_caja(z, t) en las celdas mojadas.
3. p̃_m = ρ₀gη + g Σ ρ′_m dz, con la misma discretización que `budget.pressure_tilde`. p′_m = p̃_m menos su media vertical en la cara; u′ como siempre.
4. F_bc,m = Σ u_n′ p′_m dz por cara, sumado por estrecho. Se reportan:
   - la integral en el registro (trapecio, + hacia fuera);
   - el reparto por bandas, con el mismo `validate.band_flux` y las mismas bandas;
   - η_m = F_bc,m,out / (P_bt + P_bc), con el mismo denominador de `budget.json`.
5. **Corridas:** F1 y F1c en 0–279.67 h; R1 y R1k en 0–120 h.
6. **Perfil de fondo:** ρ̄_caja(z) en el primer y en el último registro de F1, F1c, F2, F3, R1, R1k y R2.
7. **Salida:** `figures/posthoc_moving_reference.json`, de `analysis/posthoc_reference.py`. Es descriptivo; ningún veredicto congelado cambia.

**Resultado** (`figures/posthoc_moving_reference.json`; GJ, + hacia fuera; η sobre el mismo P_bt + P_bc):

| Corrida | Banda baja, fija → móvil | Suma resuelta, fija → móvil | η, fija → móvil | η resuelto, fija → móvil |
|---|---|---|---|---|
| F1 | −1.144 → **−0.750** | −0.030 → **+0.047** | −6.07 % → **−3.63 %** | −0.154 % → **+0.244 %** |
| F1c | −0.010 → +0.106 | +0.078 → +0.056 | 0.19 % → 0.44 % | 0.212 % → 0.152 % |
| R1 (0–120 h) | +0.183 → +0.186 | +0.103 → +0.137 | 0.082 % → 0.093 % | 0.030 % → 0.040 % |
| R1k (0–120 h) | +0.200 → +0.197 | +0.135 → +0.159 | 0.322 % → 0.343 % | 0.131 % → 0.154 % |

**Lo que muestra, sin adornos:**
- **En F1, la deriva explica cerca de un tercio de la banda baja hacia dentro, no su signo.** Con la referencia móvil, la banda baja sigue entrando (−0.750 GJ) y sigue dominando: |F_low| / Σ₅|F_b| = 0.83, todavía > 0.25. η(F1) sigue siendo negativo (−3.63 %). Con ninguna de las dos referencias H4a sería concluyente, ni el η de F1 sería positivo.
- **La suma de las bandas resueltas de F1 sí depende de la referencia:** pasa de −0.030 a +0.047 GJ, sobre todo por la banda de Kelvin (+0.018 → +0.092 GJ). En un registro de 280 h, la deriva no es solo de frecuencia cero: una rampa se filtra a las primeras bandas resueltas. Las fracciones φ_b de F1 no son robustas a la elección de la referencia.
- **En la bahía (R1, R1k) la referencia apenas importa:** el perfil cambia ≤ 0.006 kg m⁻³ y η cambia un 6–14 %. El cociente de H4d R1k/R1 sería 3.68 con la referencia móvil, frente al 3.93 congelado; es descriptivo.
- **Alcance de la comprobación:** la referencia móvil solo quita la parte horizontalmente uniforme de la deriva (la media de la caja). Una deriva con estructura horizontal, o un intercambio medio real, quedan dentro. Lo que queda en la banda baja de F1 (−0.75 GJ) **no está demostrado que sea un intercambio físico**: tampoco está demostrado que no lo sea.
- Ningún veredicto congelado cambia.
