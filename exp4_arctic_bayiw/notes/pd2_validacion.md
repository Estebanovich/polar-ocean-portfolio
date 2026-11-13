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
| V2a: pico cercano a la inercial en la plataforma | ≤ 1Δf de f | pico a 1.136 f, **3.11 Δf** (Δf = 1/280 h) | **FALLA** | (c), provisional |
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

**Clasificación: (c), provisional.**
- El dominio (caja de 100 × 60 km, R5, con fronteras abiertas y esponja; la estación queda a ≥ 20 celdas de ellas) no tiene una plataforma con W ≫ R₁. La respuesta cercana a la inercial de una región acotada y estratificada se desplaza al azul, por los modos de Poincaré: ω = (f² + c_n²(π/W)²)^½ (Gill 1982 §10.4 [verificar]).
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
