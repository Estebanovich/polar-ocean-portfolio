# PD2 — Validación V2–V4 del fiordo: qué falló y opciones (Task 6, parte B)

Fecha: 2026-11-12 (noche). Material para que Esteban decida en el PD2. **No se ha seguido más allá de la evaluación**: ni hipótesis (Task 7) ni figuras 2–4.

**Archivos:**
- veredictos: `figures/validation.json`, que escribe `analysis/validate.py evaluate` con los umbrales congelados del PD0 §4;
- figura: `figures/fig1_validation.{png,pdf}`;
- diagnósticos *post hoc*: `figures/pd2_v2_diagnostics.json` y `figures/pd2_v4_diagnostics.json`, de `analysis/pd2_v2_diagnostics.py` y `pd2_v4_diagnostics.py`.

**Lo que no se ha tocado:** ningún umbral ni método. Los diagnósticos de abajo se escribieron **después** de ver los fallos, están etiquetados *post hoc* y no cambian ningún veredicto.

**Corridas:** F1, V3 y V4, todas de 280 h y con estado OK en todas las etapas (`require_ok`). La comprobación del hielo de F2 y F4 pasa: el hielo queda a ≥ 146 y ≥ 119 celdas de las fronteras abiertas, en 840 registros (`figures/ice_ob_check.json`).

## Resumen

| Criterio | Umbral (PD0 §4) | Resultado | Veredicto | Clase |
|---|---|---|---|---|
| V2a: pico cercano a la inercial en la plataforma | ≤ 1Δf de f | pico a 1.136 f, **3.11 Δf** (Δf = 1/280 h) | **FALLA** | (c), definitiva (ver "Resultado 1"); la causa del desplazamiento al azul no está identificada |
| V2b: velocidad del frente de Kelvin (F1) | 0.28–0.62 m s⁻¹, ≥ 3 estaciones sucesivas | **1.65 m s⁻¹**, 12 estaciones sin huecos | **FALLA** | (a) |
| V3: presencia (decaimiento) | 1–8 km | 2.50 km | pasa | — |
| V3: sentido ciclónico | hacia el este por la costa sur | hacia el este (dφ/ds = 2.9×10⁻⁴ rad m⁻¹) | pasa | — |
| V3: velocidad de fase | 0.21–0.84 m s⁻¹ | 0.48 m s⁻¹ | pasa | — |
| V4a: flujo baroclínico neto hacia fuera, F1 frente a V4 | ≤ 25 % | F1 −1.17 GJ (entra) y V4 +0.32 GJ (sale): **128 %**, con cambio de signo | **FALLA** | (c) |
| V4b: fracción por banda | ≤ 0.10 | diferencias de 0.43 a 1.60 | **FALLA** | (c) |

V3 reproduce S&F14: c = 0.48 m s⁻¹ frente a sus 0.40–0.45 y los 0.42 de la teoría de dos capas; radio de 2.5 km frente a sus 2–3 km. S&F14 no está en la biblioteca local [verificar].

## V2a — el pico inercial está desplazado al azul en todo el dominio

**Qué se midió.** En la estación de la plataforma (j, i) = (114, 62), a 467.5 km E y 8616.9 km N, con 153 m de profundidad (regla fijada antes del análisis): espectro de energía cinética horizontal en superficie, una ventana Hann sobre 0–280 h. El pico dentro de [0.85, 1.20) f cae a 0.0929 cph, cuando f = 0.0817 cph. Son 3.11 Δf.

**Evidencia *post hoc*** (`pd2_v2_diagnostics.json`; el mismo estimador en cada cuarta columna mojada, 3567 columnas):
- El pico cercano a la inercial está desplazado al azul **en todo el dominio**: la mediana es 1.09 f (percentiles 5–95: 1.00–1.18 f) y solo el 5 % de las columnas tiene el pico a ≤ 1Δf de f. No es un artefacto de la estación elegida.
- En la estación gira en **sentido horario**: la energía horaria cercana a la inercial es 5×10³ veces la antihoraria, con el pico horario también en 1.14 f. Es movimiento inercial, desplazado al azul.

**Clasificación: (c), provisional** (definitiva tras la prueba de columna; ver "Resultado 1").
- *Corrección (2026-11-13, ronda de correcciones 1 de la parte B).* Aquí decía que el dominio "no tiene una plataforma con W ≫ R₁" y atribuía el desplazamiento a los modos de Poincaré de una región acotada. Las dos afirmaciones eran incorrectas:
  - con W = 60 km y R₁ ≈ 2.5–3.6 km, W ≫ R₁ sí se cumple;
  - la relación ω = (f² + c₁²(π/W)²)^½ (Gill 1982 §10.4 [verificar]) da solo ω/f = 1.008–1.018 con W = 60 km y c₁ = 0.35–0.52 m s⁻¹. Llegar a 1.09 f exigiría W = 18–26 km, y a 1.14 f, W = 14–21 km.

  El ancho del dominio, por sí solo, no explica el desplazamiento. Las causas candidatas, sin probar, están en "Resultado 1".
- El propio PD0 §3 anticipaba el desplazamiento al azul de esta banda (Alford et al. 2016 [verificar]). Pero el criterio exigía f en la plataforma, y no se cumple.
- **No está descartada la clase (a)**, un error del modelo en la frecuencia inercial, por ejemplo por el paso de tiempo de Coriolis. Con fΔt = 0.017 el error de fase de AB2 es del orden de 10⁻⁴ f, así que no lo espero, pero no está comprobado.

**Opciones:**
1. **Declarar la falla (c)** con la evidencia de arriba. Costo: 0.
2. **Prueba de columna:** un dominio pequeño y horizontalmente uniforme, con la misma configuración (Δt 120 s escalonado, KPP, perfil de verano) y el mismo pulso. La oscilación inercial debe salir en f a ≤ 1Δf. Separa (a) de (c). Costo: unos minutos de cómputo (np = 1) más ~1 h de preparación. Recomendada.
3. **Un dominio mayor**, con plataforma abierta. No cabe en esta máquina con la memoria del build (PD1 §5) y cambiaría el dominio de todas las F. No recomendada.

## V2b — la métrica del frente no mide una onda que se propaga

**Qué se midió.** Doce estaciones de la costa sur cada 5 km, de 0 a 59 km, sin huecos. Los tiempos de llegada (la mitad del máximo de |ζ| en la picnoclina, a −6.8 m) son 2.0, 2.5, 2.8, 2.9, 3.1, 3.7, 3.1, 3.7, 3.5, 3.1, 3.4 y 9.8 h. El ajuste da c = 1.65 m s⁻¹.

**Evidencia** (*post hoc*):
- 11 de las 12 estaciones cruzan la mitad del máximo **durante el pulso de viento** (centro a 2.5 h, σ = 1 h; PD0 §7). Los tiempos no crecen con la distancia.
- |ζ| llega a 1.4–1.9 m, y su máximo cae entre 2.7 y 25 h.
- Con viento espacialmente uniforme, cada estación responde localmente (Ekman en la costa) al mismo tiempo. El "primer frente" es la respuesta forzada, no una onda de Kelvin que se propaga.
- Como **descriptivo** (*post hoc*, no criterio), el desfase de máxima correlación del desplazamiento filtrado en la banda de Kelvin da 0.64 m s⁻¹. Los desfases son ruidosos (de −1 a 13 h).

**Clasificación: (a), error nuestro en la definición de la métrica prerregistrada.** El tiempo de llegada del primer frente no está definido cuando el forzamiento es simultáneo en todo el fiordo. Es la misma clase que el control del ancho de Mariasundet (Task 1b). El modelo no está en cuestión: V3 mide en el mismo fiordo una onda de Kelvin interna de 0.48 m s⁻¹, dentro del rango de V2b.

**Opciones:**
1. **Declarar la falla (a)** y citar V3 como la medida de la velocidad de Kelvin en el modelo. Costo: 0.
2. **Enmienda fechada** del PD0 con otra métrica de V2b, por ejemplo el desfase de correlación en la banda de Kelvin. Tendría que etiquetarse como *post hoc*, porque ya se vio el resultado (0.64 m s⁻¹, fuera del rango por 0.02). Costo: análisis, ~1 h. Valor como validación: bajo, porque se elegiría después de ver los datos.

## V4 — el flujo neto por los estrechos no converge entre 400 y 200 m

**Qué se midió.** Trabajo de presión baroclínico de Kang y Fringer en las caras de los estrechos, el de H4 tras las Enmiendas 1–2, integrado en 0–280 h.
- **Caja de V4:** a 400 m, `budget.fjord_box` no cierra la caja. La línea de Akselsundet termina al oeste de la punta de Akseløya, que a esa malla se mueve, y la caja se escapa por el flanco oeste de la isla, con extensiones de 1 a 3 celdas.
- **Solución:** la caja de V4 es la de 200 m llevada a 400 m (`validate.coarse_box`: una celda gruesa entra si está mojada y al menos la mitad de sus hijas está en la caja fina). Da 532.2 km² frente a 531.6, con 10 caras abiertas, todas en los dos estrechos.
- Es un método, no del PD0. Lo decidí después del fallo de la caja y **antes de calcular ningún flujo de V4**.

**Evidencia** (*post hoc*; `pd2_v4_diagnostics.json`; en GJ y + hacia fuera):

| | Neto | Bruto ∫\|F\|dt | Neto/bruto | Banda baja (no resuelta) | Suma de las 4 bandas resueltas | Caras | Niveles mojados (máx.) |
|---|---|---|---|---|---|---|---|
| F1 Akselsundet | −0.975 | 1.573 | −0.62 | −0.997 | +0.020 | 6 | 21 |
| F1 Mariasundet | −0.196 | 0.208 | −0.94 | −0.147 | −0.049 | 8 | 11 |
| V4 Akselsundet | +0.318 | 0.381 | +0.84 | +0.165 | +0.153 | 4 | 20 |
| V4 Mariasundet | +0.005 | 0.022 | +0.25 | +0.007 | −0.002 | 6 | 6 |

- A 200 m, el flujo neto es **casi entero de la banda baja**, el intercambio medio de periodo > 140 h: −1.14 de −1.17 GJ. Las cuatro bandas resueltas suman −0.03 GJ. Por eso las fracciones de V4b (banda / suma de las resueltas) están **mal condicionadas**: salen fracciones de −0.60 y +1.45.
- El flujo bruto cae ×4 en Akselsundet y ×10 en Mariasundet al pasar a 400 m. Mariasundet pierde la mitad de sus niveles mojados.
- El PD0 esperaba que Mariasundet fallara sola. Aquí también falla Akselsundet, y el neto total cambia de signo.

**Clasificación: (c), elección numérica (resolución).** A 400 m no se resuelve el intercambio por los estrechos. Además, el neto es un residuo pequeño de un intercambio medio sensible a la geometría.

**Consecuencia prerregistrada (PD0 §4):** "Si falla, se declara la limitación (spec §6)."

**Aviso para la Task 7 (descriptivo, no es una evaluación de H4a):** con estos números, la banda baja lleva más del 25 % del flujo total de F1. Por la regla del PD0 §5, H4a saldría "no concluyente con T_fin". Lo dejo dicho para que la decisión del PD2 lo tenga en cuenta; H4a no se ha evaluado.

**Opciones:**
1. **Declarar la limitación**, como manda el PD0. Costo: 0.
2. **Prueba a 100 m de los estrechos:** no es posible en esta máquina (PD1 §5).
3. **Reportar, como descriptivo y *post hoc*, el flujo de las bandas resueltas por separado y el bruto**, para ver si la partición por bandas, y no el neto, es robusta a la resolución. Costo: ~1 h de análisis; las lecturas ya están hechas en los diagnósticos.

## Qué necesita Esteban decidir

1. V2a: ¿declarar (c) tal cual, o antes correr la prueba de columna (minutos) para descartar (a)?
2. V2b: ¿declarar (a) y apoyarse en V3, o una enmienda *post hoc*?
3. V4: la limitación se declara por regla. ¿Se añade el descriptivo de bandas resueltas y bruto?
4. Para la Task 7: con la banda baja dominando el flujo neto, H4a probablemente salga "no concluyente". ¿Se sigue con H4b–H4d como están?

---

## Decisiones del PD2 (Esteban, 2026-11-13)

1. **V2a:** opción 2, la prueba de columna. La regla se fijó antes de correrla:
   - pico a ≤ 1Δf de f → la frecuencia inercial del modelo es correcta y el desplazamiento al azul de F1 es (c), el dominio acotado;
   - pico a > 1Δf → (a), y se investiga.

   La prueba es *post hoc*: cambia la clase de V2a, no su veredicto (FALLA).
2. **V4:** se declara la limitación (regla del PD0 §4), y se añade un análisis descriptivo *post hoc*: el bruto y cada una de las cuatro bandas resueltas, F1 frente a V4, por estrecho. V4a y V4b siguen en FALLA.
3. **V2b y Task 7:** nada por ahora. Esteban pidió antes una búsqueda bibliográfica, que lleva el controlador. No se toca la clase de V2b ni se escribe ninguna enmienda.
4. **V2b, cierre (Esteban, 2026-11-13, después de la consulta bibliográfica).**
   - **Veredicto:** FALLA, **clase (a)**. La métrica prerregistrada del tiempo de llegada del frente no está definida con un forzamiento espacialmente uniforme. Apoyo (líneas comprobadas en la biblioteca de la tesis, `~/Documents/TESIS_InternalWaves-2/Literatura/`):
     - Gill (2016), §10.10, pp. 399–401, ec. 10.10.21 (`Gill_2016 … .md`, l. 15449–15486; la ecuación está en la l. 15483). Con un viento uniforme de duración finita, la respuesta en la costa es A = 1 − cos(πt), **independiente de x**, hasta que llega la señal del extremo de la costa.
     - Yao (1986), Trinity Bay, doi:10.1080/07055900.1986.9649249 (`Yao_1986 … .md`, l. 553–558). La respuesta forzada tiene alrededor de **la mitad** del gradiente de fase a lo largo de la costa que la onda libre, y "it is not possible to distinguish between a free wave and a forced response from the observed alongshore phase gradients alone".
   - **Descriptivo *post hoc*, no es validación** (`figures/pd2_v2b_free_wave.json`, `analysis/pd2_v2_diagnostics.py --free-wave`).
     - **Método:** el de la literatura para la velocidad del periodo libre: Antenucci et al. (2000); Hodges et al. (2000), l. 798–800, desplazamiento de una isoterma filtrado en una banda; Beletsky et al. (1997), l. 424–428, seguir el máximo y el mínimo de la temperatura alrededor de la cuenca.
     - **Reglas fijadas antes de calcular:**
       - la banda de Kelvin/subinercial del PD0 §3 (el controlador citó §5; las bandas están en §3);
       - solo t ≥ 6 h (centro del pulso 2.5 h + 3σ);
       - las mismas 12 estaciones de la costa sur;
       - el desplazamiento de la picnoclina en k\* (−6.8 m).
     - **Dos detalles de implementación**, fijados con la prueba sintética **antes de calcular F1**:
       - el registro se corta en t ≥ 6 h **antes** de filtrar, porque filtrar el registro entero esparce el pulso forzado más allá de las 6 h (ondulación del filtro);
       - un extremo cuenta si su magnitud es ≥ ½ del máximo de la estación.
     - **Resultado:**
       - se sigue la primera cresta;
       - llegadas en 9.0, 10.3, 12.0, 10.7, 11.7, 12.3, 14.7, 11.7, 13.0, 21.7, 11.0 y 14.0 h a lo largo de 0–59 km;
       - **c = 0.92 ± 0.42 m s⁻¹** (pendiente y error típico del ajuste);
       - como referencia: V3, 0.48 m s⁻¹; teoría, c₁ = 0.35 (H̄) a 0.52 m s⁻¹ (H_máx).
     - **Lectura:** después del pulso, las crestas siguen llegando casi a la vez (9–15 h en 11 de 12 estaciones). Con la incertidumbre del ajuste, el resultado no separa una onda libre de la respuesta forzada, como predicen Gill (2016) y Yao (1986). Es solo descriptivo.
     - La correlación con desfase *post hoc* anterior (0.64 m s⁻¹) abarca el periodo forzado, así que se espera que esté sesgada hacia arriba (Yao 1986).
   - `figures/validation.json` recoge la clase (`V2b.classification`) y un bloque descriptivo aparte (`V2b.descriptive_free_wave`). El campo del veredicto no cambia.

## Resultado 1: prueba de columna (V2a)

**Configuración** (`scripts/column_test.py`, build `config/fjord_vmf_column`):
- **Código:** el de `fjord_vmf`, con las mismas copias modificadas de MITgcm; solo cambia `SIZE.h`: 10 × 10 × 40 celdas y 1 rango. La prueba `test_fjord_configs_share_code_but_size_h` exige que el resto sea idéntico.
- **Namelists:** los de F1, copiados (`data`, `data.free`, `data.kpp`, `data.diagnostics`, `eedata`): Δt 120 s escalonado, KPP y la mezcla del PD0, EOS lineal ajustada, f0 de 78°N, salida cada 1200 s, y etapa forzada de 0–60 h más libre hasta 280 h.
- **Fronteras:** **doblemente periódicas**: sin tierra y con `useOBCS = .FALSE.`, así que MITgcm aplica su periodicidad por omisión en x e y (`output.txt`: "pkg/obcs compiled but not used").
- **Columna:** fondo plano a 153.4 m, la profundidad de la estación de V2a. T y S son la columna inicial de F1 en esa estación. El esfuerzo es el registro de F1 en esa estación, el mismo pulso uniforme, repartido igual en todo el dominio.
- **Pipeline:** no se pasó por `cases.yaml`/`make_case`, porque `make_fjord_case` construye la batimetría IBCAO y no es barato adaptarlo. Es una configuración mínima aparte que **copia las entradas reales de F1**, y las dos etapas corren con `pfcommon.campaign.run_stages`, como en la campaña.

**Corrida:** `V2a_column`, OK en las dos etapas (8400 pasos). Build: 1 min 57 s. Modelo: 13 s + 59 s (76 s de pared en total), np = 1, mientras R3 corría con 4 rangos.

**Resultado** (`figures/pd2_column_test.json`; el mismo estimador que V2a: energía cinética horizontal en superficie, una ventana Hann en 0–280 h; el campo se llama `within_1df`, porque es el resultado de la regla, no un veredicto):
- pico a **1.005 f**, a **0.11 Δf** de f;
- la columna queda horizontalmente uniforme (desviación máxima de U en superficie: 0.0).

**Clase final de V2a: (c)**, por la regla fijada antes de la prueba. El veredicto sigue siendo **FALLA**.

**Qué excluye la prueba y qué no** (redacción corregida el 2026-11-13, ronda de correcciones 1 de la parte B; la versión anterior atribuía la causa al dominio acotado, más de lo que la prueba muestra):
- **Excluye** un error del modelo en la frecuencia inercial: el término de Coriolis y el paso de tiempo. La columna es una losa en plano f sin fronteras.
- **No prueba nada que dependa de las fronteras.** Con `useOBCS = .FALSE.` las copias modificadas de `obcs` y la esponja no actúan, y la viscosidad lateral (Smagorinsky) no tiene gradientes horizontales sobre los que actuar.
- **La causa del desplazamiento al azul de F1 no está identificada:** pico de 1.136 f en la plataforma; en todo el dominio, mediana 1.09 f y percentiles 5–95 de 1.00–1.18 f, en sentido horario. Candidatas, sin elegir ninguna:
  1. la geometría acotada y los modos laterales;
  2. las fronteras abiertas y la esponja;
  3. la vorticidad de fondo o el flujo medio.
- **El ancho del dominio, por sí solo, no basta:** ω/f = (1 + (c₁π/(W f))²)^½ = 1.008–1.018 con W = 60 km y c₁ = 0.35–0.52 m s⁻¹. Para 1.09–1.14 f harían falta W = 14–26 km (Gill 1982 §10.4 [verificar]).

Actualizado en `figures/validation.json` (`V2a.classification`), que ahora se reconstruye con `validate.py pd2-merge`.

*Nota de cita (Task 8, 2026-11-14).* La relación de Poincaré de esta nota cita «Gill 1982 §10.4 [verificar]». Comprobada en la biblioteca de la tesis (`Gill_2016 … .md`, l. 14593–14600), la fuente es **Gill 1982 §10.3, ec. 10.3.5**: ω_n² = f² + n²π²gH/W² en un canal de ancho W, con c₁ en lugar de √(gH) para el modo interno. El §10.4 es el de las ondas de Kelvin. Los números no cambian; la cita se corrigió en el código (`validate.py`, `pd2_v2_diagnostics.py`, `theory_numbers.py`) y en `validation.json`. El texto de arriba y el del PD0 §2 no se editan.

## Resultado 2: V4, limitación declarada y descriptivo *post hoc*

**La limitación se declara** (PD0 §4): a 400 m el flujo baroclínico neto por los estrechos no converge (cambia de signo), así que V4 no puede acotar la incertidumbre de resolución de los cocientes de H4. V4a y V4b quedan en **FALLA**.

**Descriptivo *post hoc*** (`figures/pd2_v4_descriptive.json`, de `analysis/pd2_v4_diagnostics.py --descriptive`). GJ, + hacia fuera. "Peso" es |F_b| / Σ|F_b| sobre las cuatro bandas resueltas.

| Estrecho | | Bruto ∫\|F\|dt | Kelvin/subinercial | Cercana a la inercial | Superinercial | Seiche | Banda baja (no resuelta) |
|---|---|---|---|---|---|---|---|
| Akselsundet | F1 (200 m) | 1.573 | +0.0520 (peso 0.49) | −0.0340 (0.32) | +0.0108 (0.10) | −0.0093 (0.09) | −0.997 |
| | V4 (400 m) | 0.381 (×0.24) | +0.1534 (0.85) | −0.0121 (0.07) | +0.0131 (0.07) | −0.0012 (0.01) | +0.165 |
| Mariasundet | F1 (200 m) | 0.208 | −0.0341 (0.69) | −0.0091 (0.19) | −0.0006 (0.01) | −0.0053 (0.11) | −0.147 |
| | V4 (400 m) | 0.022 (×0.11) | −0.0021 (0.60) | −0.0006 (0.17) | +0.0003 (0.08) | +0.0005 (0.15) | +0.007 |

**Lectura** (descriptiva, sin umbral):
- **Akselsundet:** el **signo** de las cuatro bandas resueltas es el mismo a 200 y a 400 m, pero los **pesos no**:
  - la banda cercana a la inercial pasa de 0.32 a 0.07;
  - la de Kelvin, de 0.49 a 0.85;
  - diferencia máxima de peso 0.36.
- **Mariasundet:** los pesos se parecen (diferencia máxima 0.10), con Kelvin primero y la cercana a la inercial después en las dos mallas. Pero todo el flujo es 10 veces menor a 400 m, y las bandas superinercial y de seiche cambian de signo con valores minúsculos.
- **En ningún estrecho la partición es robusta**:
  - en Akselsundet, que lleva el 88 % del bruto de F1, cambia el reparto entre bandas;
  - en Mariasundet, el estrecho apenas existe a 400 m (6 niveles mojados frente a 11).
- La banda baja domina el neto en F1 (−1.14 de −1.17 GJ) y cambia de signo a 400 m.


## Cómo se reconstruyen los bloques del PD2 en `validation.json`

`validate.py evaluate` solo escribe los veredictos. Los bloques del PD2 (`classification`, `pd2_decision`, `column_test`, `descriptive_free_wave`, `note_on_lag_correlation`) los reconstruye **`validate.py pd2-merge`** a partir de:
- `pd2_column_test.json`;
- `pd2_v2_diagnostics.json`;
- `pd2_v2b_free_wave.json`;
- `pd2_v4_diagnostics.json` y `pd2_v4_descriptive.json`;
- los textos fijos de clasificación, que están en el código.

No cambia ningún veredicto (prueba `test_pd2_merge_rebuilds_the_pd2_blocks`). El orden documentado es `evaluate` y después `pd2-merge`.

*Nota de cita (Task 8, ronda de correcciones 1, 2026-11-14).* La decisión del PD2 cita «Gill (2016)». El archivo `Gill_2016 … .md` de la biblioteca de la tesis es una reimpresión de Academic Press de la edición de 1982 (International Geophysics Series 30; copyright 1982, l. 15–25 del archivo), con las mismas secciones y páginas. La cita correcta es **Gill (1982)**: §10.10, pp. 399–401, ec. 10.10.21 para V2b, y §10.3, ec. 10.3.5 para la relación de Poincaré. Así se cita en el README, en `validation.json` y en el código. El texto de arriba no se edita.
