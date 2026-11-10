# Comparación post hoc con MITgcm-JPL (Nakayama y Menemenlis), Ocean0 a 2 km

Fecha: 2026-11-09. **Post hoc**, a pedido de Esteban; no es un criterio de pd0 y no cambia ningún veredicto. Script: `analysis/jpl_compare.py` → `figures/jpl_comparison.json`. Datos de JPL: `data_external/osf_jpl/` (OSF, DOI 10.17605/OSF.IO/QNA89; ver `data_external/SOURCES.md`).

## Qué se compara

- **Nuestro:** `ocean0_2km_r21` (ΓT = 0.030959, ley de u\* de R21), meses 19–24 de 2 años (ventana `last6` de pd0).
- **JPL:** `Ocean0_yosi_COM_MITgcm_JPL_new.nc`, un año de medias mensuales (12 registros); se toman los registros 7–12 como «últimos 6 meses». ΓT = 0.0325, ΓS = 0.0009286, C_D,top = 0.0025 (resumen PDF).
- **Regiones,** cada modelo con su propio calado y su propia máscara: plataforma completa, y región profunda = plataforma con base del hielo < −300 m (AD16 Eq. 37). La plataforma de JPL son las celdas con `meltRate` no enmascarado y `iceDraft` < 0. Medias ponderadas por área (las dos mallas son uniformes de 2 km).
- **Malla común:** los centros coinciden exactamente (x, y de JPL = XC, YC nuestros; comprobado en el script). Nuestra fila y = 81 km (pared cerrada) se descarta.

## Resultados

| Cantidad | Región | Nuestro | JPL | Nuestro/JPL − 1 |
|---|---|---|---|---|
| Melt medio (m a⁻¹) | profunda | 30.00 | 31.35 | −4.3 % |
| Melt medio (m a⁻¹) | plataforma | 10.19 | 10.33 | −1.3 % |
| u\* medio (m s⁻¹) | profunda | 0.00365 | 0.00397 | −7.9 % |
| u\* medio (m s⁻¹) | plataforma | 0.00254 | 0.00282 | −10.0 % |
| Thermal driving (°C) | profunda | 0.799 | 0.703 | +13.6 % |
| Thermal driving (°C) | plataforma | 0.326 | 0.259 | +25.8 % |
| ΓT | — | 0.030959 | 0.0325 | −4.7 % |
| Celdas de plataforma | — | 2689 | 3170 | — |
| Celdas profundas | — | 618 | 687 | — |

**Correlación espacial de los mapas de melt** (medias temporales, malla común de 2 km):
- r = 0.911 en las 3170 celdas donde alguno de los dos tiene plataforma (melt = 0 fuera de ella);
- r = 0.938 en las 2689 celdas donde los dos tienen plataforma.

## Lectura (post hoc, sin veredicto)

- Con la ley de u\* corregida, nuestro u\* profundo (0.00365) queda a 8 % del de JPL (0.00397). Con la ley anulada era 0.0109 (corrida anulada `ocean0_2km`, meses 19–24), 2.7 veces el de JPL. El piso de JPL es 0.0005 m s⁻¹ (mínimo de `frictionVelocity`), el mismo √C_D·u_tidal de AD16 Eq. 27.
- Nuestro ΓT_cal es 4.7 % menor que el de JPL, y aun así tenemos menos u\* y más thermal driving. Los dos modelos llegan a ~30 m a⁻¹ en la región profunda con un balance distinto entre velocidad y temperatura. No se atribuye causa aquí.
- **Máscaras.** JPL tiene 481 celdas de plataforma que nosotros no tenemos (nosotros no tenemos ninguna que ellos no tengan):
  - 205 son secas en nuestro modelo (x = 459–567 km, columna de agua de JPL ≈ 100 m): el lado de la línea de apoyo. Candidatas, no verificadas: nuestra regla de mayoría de subceldas flotantes de 1 km (R9) y la columna mínima de 2 celdas, frente a la interpolación lineal de JPL a 2 km. JPL tiene allí 7.8 m a⁻¹ de media;
  - 276 son mojadas pero sin hielo en nuestro modelo (calado 0), junto a las paredes laterales (y ≤ 7 km o ≥ 73 km, x = 515–639 km). JPL tiene allí hielo fino y melt bajo (0.6 m a⁻¹ de media).
  
  Esto entra en las medias de plataforma y explica parte de la diferencia de celdas profundas (618 contra 687). No se clasifica aquí.

## Diferencias de definición y de configuración

- **Thermal driving.** La de JPL es la variable `thermalDriving` de AD16: T_w − T_zd, con T_w promediada sobre una capa de 20 m bajo el calado. No tenemos diagnóstico de T_w ni de T_b, así que la nuestra se invierte del balance de calor AD16 Eq. 24 con κ_i = 0 (exacto en nuestro código): T_w − T_b = −⟨SHIfwFlx⟩·L / (ρ_ref·c_w·⟨SHIgammT⟩), con medias mensuales. Es la media de T_w − T_b ponderada por γT, no la media simple. La diferencia entre las dos definiciones depende de la covarianza temporal y espacial de γT y T_w − T_b, cuyo signo no conocemos; no se cuantifica aquí, y la comparación de thermal driving vale solo como orden de magnitud.
- **Velocidad de u\*.** JPL promedia u sobre una capa de 20 m y los 4 vecinos horizontales, y excluye las velocidades sin flujo normal en las caras verticales del hielo (toma los 20 m desde el fondo de la capa límite). Nosotros usamos `SHELFICEboundaryLayer`: la celda superior más la parte de la siguiente que completa el espesor de una celda (20 m), con u y v llevados al centro de la celda (`metrics.boundary_layer_velocity_sq`).
- **Configuración** (resumen PDF de JPL contra la nuestra):
  - advección: JPL usa DST de 3.er orden con limitador (33); nosotros, el limitador de 2.º orden (77);
  - convección: JPL usa ajuste convectivo; nosotros, `ivdc_kappa` = 0.1 (AD16 Tabla 4);
  - celdas parciales: JPL usa hFacMin = 0.25; nosotros, 0.2;
  - JPL: MITgcm checkpoint 66o, con el código de la plataforma modificado y no publicado; flujos de sal virtuales; topografía interpolada linealmente a 2 km.
- **Duración.** JPL entregó un año de Ocean0; nosotros corremos 2 años (pd0). Los meses 7–12 de JPL y los 19–24 nuestros no son la misma fase de ajuste.
- **Meses.** Los registros de JPL son meses de calendario (`time` = 0, 31, 59, … días: meses de 28–31 días), así que sus registros 7–12 cubren 184 días (julio–diciembre). Los nuestros son meses fijos de 2 628 000 s (30.42 días, pd0), así que nuestros 6 meses cubren 182.5 días. Las medias de JPL sobre la ventana se toman sin ponderar por la duración de cada mes (los meses de 30 y 31 días difieren en ~3 %).

(La serie R20, que probó 33, hFacMin = 0.25 y ajuste convectivo de uno en uno, está anulada por R21 y no se usa aquí.)

# Ocean1 (post hoc, 2026-11-10)

Mismo procedimiento que arriba, a pedido de Esteban; **post hoc**, no cambia ningún veredicto (los de pd0 están en `pd1_ocean1.md`). `figures/jpl_comparison.json` tiene ahora una entrada por caso (`ocean0`, `ocean1`); los números de Ocean0 de arriba se reprodujeron sin cambios.

- **Nuestro:** `ocean1_2km_r21`, año 20 (registros 229–240).
- **JPL:** `Ocean1_yosi_COM_MITgcm_JPL_new.nc` tiene 240 registros mensuales de calendario (`time` de 0 a 19.92 años). Se toman los últimos 12 (`time` ≥ 19 años) como su año 20.
- Mismas regiones, máscaras y definiciones que en Ocean0, incluida la advertencia sobre la thermal driving: la nuestra está ponderada por γT, la de JPL es la media simple, y el signo de la diferencia es desconocido.

| Cantidad | Región | Nuestro | JPL | Nuestro/JPL − 1 |
|---|---|---|---|---|
| Melt medio (m a⁻¹) | profunda | 30.37 | 31.55 | −3.7 % |
| Melt medio (m a⁻¹) | plataforma | 10.34 | 10.46 | −1.1 % |
| u\* medio (m s⁻¹) | profunda | 0.00371 | 0.00396 | −6.2 % |
| u\* medio (m s⁻¹) | plataforma | 0.00256 | 0.00282 | −9.3 % |
| Thermal driving (°C) | profunda | 0.801 | 0.709 | +13.0 % |
| Thermal driving (°C) | plataforma | 0.330 | 0.262 | +26.2 % |
| máx. \|ψ_bt\|, cavidad (y > 50 km, x < 640 km) (Sv) | — | 0.258 | 0.246 | +4.6 % |
| máx. \|ψ_bt\|, frente (630 < x < 650 km) (Sv) | — | 0.492 | 0.276 | +78 % |
| máx. \|ψ(x, z)\|, x < 640 km (Sv) | — | 0.185 | 0.203 | −8.8 % |
| Fuerza barotrópica (máx − mín, x < 640 km) (Sv) | — | 0.261 | 0.271 | −3.8 % |
| Fuerza de la overturning (máx − mín, x < 640 km) (Sv) | — | 0.186 | 0.222 | −16 % |
| ΓT | — | 0.030959 | 0.0325 | −4.7 % |

**Correlación de los mapas de melt del año 20:** r = 0.913 en las 3170 celdas donde alguno tiene plataforma, y r = 0.942 en las 2689 donde los dos la tienen. Las diferencias de máscara son las mismas 481 celdas que en Ocean0.

**Streamfunctions.** Las extremas se toman sobre la media del año 20 en los dos modelos. Las nuestras están en las esquinas de celda; las de JPL son las que entregaron, interpoladas a la malla estándar de AD16 (centros en x, y; niveles de 5 m). Por eso no se evalúan en los mismos puntos.
- La fuerza aquí es máx − mín de la ψ media; la de V3 es la media de los máx − mín mensuales. Para nosotros dan casi lo mismo (0.261 contra 0.261 en la barotrópica, 0.186 contra 0.186 en la overturning).
- El máximo de \|ψ_bt\| cerca del frente es 78 % mayor en nuestro modelo. Los dos valores caen dentro de V2d.

**Lectura (sin atribución).** Con la misma meta de melt, JPL tiene una overturning ~10 % más fuerte (0.203 contra 0.185 Sv; la suya cae dentro de V2b) y algo más de u\*; nosotros, más thermal driving. Es coherente con lo visto en Ocean0. Las diferencias de configuración están listadas arriba; la serie R20, que las probaba, está anulada. No se clasifica aquí: es material para PD2 (V2b y V3 overturning fallan).
