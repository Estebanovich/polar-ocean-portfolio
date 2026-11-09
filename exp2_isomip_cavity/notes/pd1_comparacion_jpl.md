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

- **Thermal driving.** La de JPL es la variable `thermalDriving` de AD16: T_w − T_zd, con T_w promediada sobre una capa de 20 m bajo el calado. No tenemos diagnóstico de T_w ni de T_b, así que la nuestra se invierte del balance de calor AD16 Eq. 24 con κ_i = 0 (exacto en nuestro código): T_w − T_b = −⟨SHIfwFlx⟩·L / (ρ_ref·c_w·⟨SHIgammT⟩), con medias mensuales. Es la media de T_w − T_b ponderada por γT, no la media simple. Como el melt crece donde γT es grande, esta ponderación puede sesgar nuestra cifra hacia arriba frente a la de JPL; no se cuantifica aquí.
- **Velocidad de u\*.** JPL promedia u sobre una capa de 20 m y los 4 vecinos horizontales, y excluye las velocidades sin flujo normal en las caras verticales del hielo (toma los 20 m desde el fondo de la capa límite). Nosotros usamos `SHELFICEboundaryLayer`: la celda superior más la parte de la siguiente que completa el espesor de una celda (20 m), con u y v llevados al centro de la celda (`metrics.boundary_layer_velocity_sq`).
- **Configuración** (resumen PDF de JPL contra la nuestra):
  - advección: JPL usa DST de 3.er orden con limitador (33); nosotros, el limitador de 2.º orden (77);
  - convección: JPL usa ajuste convectivo; nosotros, `ivdc_kappa` = 0.1 (AD16 Tabla 4);
  - celdas parciales: JPL usa hFacMin = 0.25; nosotros, 0.2;
  - JPL: MITgcm checkpoint 66o, con el código de la plataforma modificado y no publicado; flujos de sal virtuales; topografía interpolada linealmente a 2 km.
- **Duración.** JPL entregó un año de Ocean0; nosotros corremos 2 años (pd0). Los meses 7–12 de JPL y los 19–24 nuestros no son la misma fase de ajuste.

(La serie R20, que probó 33, hFacMin = 0.25 y ajuste convectivo de uno en uno, está anulada por R21 y no se usa aquí.)
