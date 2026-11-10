# PD0 — Preregistro del Exp. 2 (fijado antes de cualquier corrida)

Fecha: 2026-11-08. No hay ninguna salida de MITgcm del Exp. 2 en el momento de escribir esto.

**Fuentes**
- Asay-Davis et al. 2016 (AD16), *GMD* 9:2471.
- Yung et al. 2026, versión final (Y26), *The Cryosphere* 20:2053–2088; en `literatura/Yung_2026*.md`. Página de la revista = página del PDF + 2052. "L" = línea del `.md` final.
- Suplemento de Y26 (Y26S).
- Holland et al. 2008 (H08).

Los valores marcados **(leído de figura, aprox.)** salen de la imagen y no de un número impreso. El README los presenta con esa etiqueta.

## Definiciones comunes

- **Calendario:** 365 días; mes = 2 628 000 s.
- **"Últimos 6 meses"** de una corrida de 2 años: los registros mensuales 19–24 (índices 18–23 en base 0). **"Año 20"** de Ocean1: los registros 229–240.
- **Melt:** m a⁻¹ de agua dulce (ρfw = 1000 kg m⁻³), positivo = fusión. El signo y las unidades de `SHIfwFlx` se fijan en `pd1_config.md` (Task 2) leyendo `shelfice_thermodynamics.F`.
- **Plataforma:** celdas mojadas con base del hielo < 0 tras el vaciado (AD16 §3.1.2).
- **Melt profundo:** celdas de plataforma con base < −300 m (AD16 Eq. 37).
- **Medias:** ponderadas por área.
- **ψ barotrópica:** ψ(x, y) = −∫₀^y U dy′, con U el transporte zonal integrado en la vertical; en Sv. La overturning ψ(x, z) se integra en y.
- **Fuerza de la circulación** (Y26 Fig. 13, L1365–1369): máx − mín de ψ para x < 640 km, por mes.

## K1 — Calibración (no es validación)

- **Regla:** ⟨m⟩ profundo de Ocean0, en los últimos 6 meses, = 30 ± 2 m a⁻¹ (AD16 Eq. 37), por separado a 2 km y a 4 km.
- **Procedimiento:**
  - se muestrea ΓT ∈ {0.015, 0.03, 0.06}, con ΓS = ΓT/35;
  - se interpola en log–log (ln ⟨m⟩ contra ln ΓT, lineal por tramos) hacia 30 m a⁻¹;
  - se verifica con una corrida;
  - si el valor cae fuera de 28–32, se hace una sola iteración más con los 4 puntos.

  Si aún falla, se va a PD2 y se reporta.
- ΓT_cal,2km sirve para Ocean0 y Ocean1 a 2 km; ΓT_cal,4km, para el ensamble. No se vuelve a tocar.

## V1 — ΓT dentro del rango de ISOMIP+

- **Regla:** pasa si 0.011 ≤ ΓT_cal,2km ≤ 0.2. Es el rango COM de Y26 Tabla 2 (L472–600), idéntico al del preprint.
- **Referencia informativa, no criterio:** los modelos z de Y26 están en 0.0135–0.0325; MITgcm-JPL (Nakayama) usa 0.0325, MITgcm-BAS 0.011 y MITgcm-BAS-Coupled 0.0135.

## V2 — Ocean1 a 2 km, año 20, contra los valores publicados

Cada métrica se reporta por separado. "±" no se define en Y26: se toma el intervalo tal como está escrito.

| Métrica | Rango de paso | Fuente |
|---|---|---|
| V2a: melt medio de toda la plataforma | 4–12 m a⁻¹ | Y26 §4.3, p. 2068, L1260–1261 ("4–12 m yr−1") |
| V2b: máx. de la overturning | 0.20–0.25 Sv | Y26 §4.2, p. 2066, L1100–1102 ("most models"); COCO 0.39, y ROMS y FVCOM más débiles, quedan fuera de ese rango |
| V2c: máx. de ψ barotrópica en la cavidad (y > 50 km, x < 640 km) | 0.26 ± 0.17 → 0.09–0.43 Sv | Y26 p. 2066, L1118 |
| V2d: máx. de ψ barotrópica cerca del frente (630 < x < 650 km) | 0.41 ± 0.27 → 0.14–0.68 Sv | Y26 p. 2067, L1155 |

Para V2b, "máx." es el máximo de |ψ(x, z)| en la cavidad (x < 640 km), igual que en la Fig. 8 de Y26.

## V3 — Relación melt–circulación (Y26 Fig. 13; leído de figura, aprox.)

Se usan los puntos mensuales de Ocean1 a 2 km del año 20: el promedio de los 12 meses, (m̄, C̄), con C = fuerza de la circulación.

**Regla:** pasa si C̄ cae dentro de la banda de los modelos COM sin COCO, interpolada linealmente en m̄ entre los anclajes. Si m̄ cae fuera de [4, 12], V3 = "no evaluable" (fuera del rango de datos de la figura).

| Melt (m a⁻¹) | Barotrópica (Sv) | Overturning (Sv) |
|---|---|---|
| 5 | 0.17–0.30 | 0.10–0.18 |
| 10 | 0.25–0.45 | 0.19–0.28 |

**Informativo, sin umbral:** pendiente y ordenada de nuestro ajuste con los 240 meses, contra el ajuste multi-modelo impreso en la figura. En Ocean1 son y = 0.024x + 0.071 (barotrópica, r = 0.89) e y = 0.02x + 0.019 (overturning, r = 0.97).

## N1–N3 — Numéricos, por corrida

- **N1:** |⟨m⟩(meses 19–24) − ⟨m⟩(meses 13–18)| / ⟨m⟩(meses 19–24) < 3 %, con ⟨m⟩ = melt de toda la plataforma. Para Ocean1 es |⟨m⟩(año 20) − ⟨m⟩(año 19)| / ⟨m⟩(año 20) < 3 %. Y26 L1260 dice que todos los modelos salvo ROMS son estables desde el año 14.
- **N2:** |ΔInventario − (∫fuente − ∫sumidero rbcs)| / ∫fuente < 1 % sobre toda la corrida.
- **N3:** `require_ok`: OK, nTimeSteps completos y sin NaN.
- **N4:** la rampa con z_c = 360 y d = 720 reproduce Ocean0 con error < 1e-12 (prueba unitaria).

## H1 — Forma cuadrática (H08 §4, §6)

En el núcleo a 4 km (s = 1), para cada z_c ∈ {250, 400, 550}, con las 5 corridas de Tbot ∈ {−0.5, …, 1.5}:

1. **Ajuste cuadrático:** ⟨m⟩ = a·Tbot² + b·Tbot + c, por mínimos cuadrados. Pasa si R² ≥ 0.98.
2. **Ley de potencia:** ⟨m⟩ = k·(Tbot − T_f,GL)ⁿ, ajustada en log–log. Pasa si n ∈ [1.3, 2.0]. T_f,GL es el punto de fusión (AD16 Tabla 4: λ1–λ3) con Sbot = 34.7 y la presión de la base más profunda tras el vaciado.

H1 pasa en un z_c si pasan 1 y 2. Se reporta por z_c.

**Informativo:** el suplemento de Y26 (Fig. S1a) da exponentes de 0.99–1.70 en el transitorio de Ocean1, con otra definición de temperatura (T + 2.1 °C bajo 300 m a x = 520 km). Se cita como contexto, no como umbral.

## H2 — Profundidad de salida del agua de deshielo (H08 §4)

- **Definición:** z_out = la profundidad media del transporte saliente de trazador en x = x_front, ponderada por ese transporte, en los últimos 6 meses. x_front es la columna de celdas a mar abierto inmediatamente al norte del frente de hielo, con x máximo de la plataforma + 1 celda, en cada y. Solo cuenta el transporte con u·C > 0.
- **Regla:** D_out = −z_out (profundidad positiva, m). En el núcleo a 4 km (s = 1), para cada z_c, pasa si Spearman ρ(Tbot, D_out) ≤ −0.9: el agua de deshielo sale menos profunda cuando aumenta Tbot. Se reporta por z_c.

## Prueba de resolución (descriptiva)

Ocean0 a 2 km contra 4 km, cada uno con su ΓT_cal: se reportan las diferencias relativas del melt de toda la plataforma, de la máx. overturning y de la máx. ψ barotrópica en la cavidad, y la correlación espacial de los mapas de melt (el de 2 km promediado a 4 km). No hay umbral. Si la diferencia en el melt de toda la plataforma supera el 25 %, el README lo destaca como limitación del ensamble.

## Lo que queda fuera

- **La correlación de patrón contra la media multi-modelo,** que era V3 en la primera versión del plan. Los datos de los modelos se publicaron en OSF (Y26 L2102–2107, Tabla 3; MITgcm-JPL de Nakayama: doi:10.17605/OSF.IO/QNA89), pero no se han descargado. Si se descargan después, cualquier comparación con ellos se etiqueta **post hoc**.

## Aclaración previa a los datos (R22, 2026-11-09 16:39 CST)

Decisión de Esteban, tomada **antes de evaluar cualquier salida del ensamble** (la cola de producción empezó a las 16:07 CST; no se había leído ni evaluado ninguna salida de las corridas `ens_*`, solo su estado de ejecución). Es un agregado: el texto de arriba no cambia.

- **⟨m⟩ de H1 = melt medio de toda la plataforma** (ponderado por área, "Plataforma" de las definiciones comunes; `mean_melt_shelf`). Razones: es la métrica de H08, el plan dice "derretimiento medio" y es la definición de ⟨m⟩ de N1.
- **El ajuste con el melt profundo** (base < −300 m) se calcula y se reporta **como secundario**, con esa etiqueta. No es criterio.
- **Corridas del núcleo que fallan N1 o N2:** se quedan en el ajuste de H1 y H2 y se señalan en el veredicto. Solo una corrida que falla N3 (o que falta) es un hueco; con un hueco, ese z_c queda "no evaluable", sin imputar.
- **Figura 3(c):** el mapa del trazador integrado en la vertical es el del caso Tbot = 1.5 °C, z_c = 400 m.
