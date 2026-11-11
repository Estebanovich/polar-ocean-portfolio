# PD1 — Batimetría de Van Mijenfjorden para MITgcm (Task 1b)

Material para el PD1. La batimetría se aprueba en el PD1 junto con la configuración. Si el PD1 cambia la caja o la malla, se edita `cases.yaml` (`configs: fjord_vmf:`) y se vuelve a correr la Task 1b con sus pruebas:

```
$PY exp4_arctic_bayiw/analysis/bathy_report.py      # JSON, figura y bathyFile (no versionado)
$PY -m pytest exp4_arctic_bayiw/tests/test_bathy.py
```

Código: `gendata/bathy.py` (procedimiento), `gendata/vgrid.py` (malla vertical), `analysis/bathy_report.py` (controles y figura). Resultados: `figures/bathy_report.json` y `figures/bathy_report.{png,pdf}`. El `bathyFile` (500 × 300, float64 big-endian, 1.2 MB) va a `data_external/derived/bathy_fjord_vmf.bin`, que git ignora; lo regenera `make_case` (Tasks 3–4).

## 1. Referencia: Støylen & Fer (2014)

Støylen, E. y Fer, I. (2014), *Tidally induced internal motion in an Arctic fjord*, Nonlin. Processes Geophys. 21, 87–100, doi:10.5194/npg-21-87-2014. No está en `literatura/` (`litsearch.py "Van Mijen"` no encontró nada). Se leyó completo en la **versión en línea** (PDF de acceso abierto de Copernicus, https://npg.copernicus.org/articles/21/87/2014/). No se guardó en el repositorio: Esteban lo descarga.

Lo que reportan sobre la geometría (página y figura de la versión en línea):

| Dato | Valor | Dónde |
|---|---|---|
| Estrechos de la boca | Akseløya divide la entrada en dos: Akselsundet al norte y Mariasundet al sur | p. 89, Sect. 3; Fig. 1 |
| Ancho de Mariasundet | "relatively narrow (0.7 km wide)" | p. 96, Sect. 5 |
| Akselsundet | No dan ancho ni profundidad. La caja idealizada usa un estrecho de 20 m × 1300 m "giving a cross-sectional area similar to Akselsundet" (≈ 2.6 × 10⁴ m²) | p. 93, Sect. 4.2 |
| Corrientes en Akselsundet | pueden pasar de 2 m s⁻¹ (Bergh 2004) | p. 89 |
| Profundidad del umbral | **No la reportan.** Solo dicen que el agua de fuera "is lifted over the sill" (sección D) | p. 90; Fig. 3a |
| Profundidades de la cuenca | Profundidad mín.–máx. de las secciones CTD: A-1 61–117 m, C-1 41–110 m, C-2 42–114 m, D 33–113 m, B 34–88 m; TS2 (justo dentro de Akselsundet) 60 m | Tabla 1, p. 89 |
| Batimetría del modelo | La escala de color de la topografía llega a −120 m. Omiten "the shallow, secondary inner basin on the east" | Fig. 15, p. 96; Sect. 4.3, p. 95 |
| Fuente y resolución de su batimetría | No dicen de dónde sale la batimetría. MITgcm hidrostático, celdas recortadas (*shaved cells*), Δx = 100 m, 32 niveles (0.75 m donde el gradiente de densidad es máximo, hasta 14 m por debajo de 100 m), esponja de 64 puntos y frontera abierta al oeste en la caja idealizada | p. 92, Sect. 4.1; p. 93 |
| Advertencia | Mariasundet es estrecho frente a su malla de 100 m y "the representation of this sound in the topography matrix is coarse" | p. 96 |

Con esto se fijaron los controles **antes** de correr el procedimiento completo (en `cases.yaml`, `bathy: controls:`):

- el ancho de Mariasundet frente a 0.7 km, con ±200 m (una celda);
- el área de Akselsundet frente a 2.6 × 10⁴ m², con ±50 % (es un área "parecida", no medida);
- la profundidad máxima de la cuenca interior frente a 117 m, con ±10 %.

Como no dan la profundidad del umbral, el umbral se comparó con IBCAO a 100 m (±20 %), y la profundidad máxima y el volumen también se comparan contra IBCAO (±5 %).

## 2. Datos y sistemas de referencia

- **IBCAO v5.2** (Jakobsson et al. 2024, doi:10.1038/s41597-024-04278-w; malla doi:10.17043/ibcao-5.2). Es una ventana de 100 m en EPSG:3996 con elevación, `sid` y `tid` (`data_external/SOURCES.md`). Los GeoTIFF se leen con PIL. La georreferencia sale de las etiquetas GeoTIFF (escala de píxel, punto de amarre y `PixelIsArea`), así que los centros de celda son esquina + 50 m. Prueba con datos reales: en la boca de Isfjorden (78.1°N, 13.6°E) la profundidad es −291 m con TID 11.
- **Códigos TID** (GEBCO, https://www.gebco.net/gebco-tid-grid): 10–17 son medición directa (11 multihaz, 17 combinación de métodos directos); 40–46, indirecta (41 interpolado por algoritmo); 70–72, desconocida. IBCAO usa además 47 (2 celdas en la ventana) y **100**. El 99.5 % de las celdas con TID 100 está sobre el nivel del mar, así que se trata como tierra [verificar en la documentación de IBCAO 5].
- **Costa: NPI S100 Kartdata** (Norwegian Polar Institute, doi:10.21334/npolar.2014.645336c7, CC BY 4.0), polígonos `S100_Land_f` (6482 polígonos con `OBJTYPE = Landareal`). Se leen con pyshp y se unen con shapely. El dato está en ETRS89 / UTM 33N (EPSG:25833). Se usa como EPSG:32633 (WGS 84), porque la diferencia es menor de 1 m en Svalbard [verificar]. La licencia pide acreditar al NPI; la figura lo hace en su pie. El texto `_Lisens_License.txt` debe acompañar cualquier copia de datos derivados, lo que importa si algún día se publica el `bathyFile`.

## 3. Caja y mallas (decisión de la Task 1b; ruling R2)

- **Caja:** UTM 33N sin rotar, x = 455–555 km e y = 8594–8654 km. Son **100 × 60 km** a 200 m, 500 × 300 celdas, de las que 56 531 son de mar. El plan pedía unos 45 km de norte a sur, pero no alcanzan:
  - Van Mijenfjorden con Vallunden llega a y ≈ 8651.5 km.
  - Van Keulenfjorden y Recherchefjorden, los fiordos laterales de Bellsund, bajan hasta y ≈ 8599.5–8600 km en la costa del NPI (8598 km con la región R7 de §4).
  - Con 45 km, la frontera sur cruzaría un fiordo, con una frontera abierta y una esponja dentro de él, o habría que cerrar a mano un fiordo entero.

  Con 60 km ambos quedan completos y el costo sube un 33 % frente a 100 × 45. **Alternativa para el PD1:** 100 × 45 con Van Keulenfjorden y Recherchefjorden cerrados como tierra en su boca.
- **Fronteras abiertas** O, N y S. Todo el borde oeste es mar. El borde norte es mar para x < 467.8 km y el sur para x < 474.4 km; el resto de esos bordes es tierra. El borde este es tierra (la cabeza del fiordo está en x ≈ 551.5 km). Las esquinas NO y SO son mar, y `open_boundaries` lo exige: si cae tierra en una esquina abierta, da error.
- **Malla vertical** (`vgrid.stretched_delR`): dz_k = dz_top + (dz_bot − dz_top)(k/(nz−1))^p, con nz = 40, dz_top = 1.5 m y p = 3.5. dz_bot se despeja de Σdz = H, y se redondea a 1 mm metiendo el residuo en el último nivel.
  - H = 240 m: es la profundidad máxima del mar final del modelo (238.8 m) redondeada a la decena.
  - Las celdas de IBCAO de hasta 261.7 m de la esquina SO caen en la esponja sur, y la copia de la esponja las reemplaza (99 celdas se recortaron primero al fondo).
  - Quedan 25 niveles en los primeros 60 m, de 1.5 a 5.1 m, casi uniformes en la picnoclina para KPP; 21.1 m en el fondo; la razón entre niveles vecinos es ≤ 1.11.
  - Comparación: Støylen & Fer usaron 32 niveles, de 0.75 m a 14 m.
  - Las pruebas verifican que Σdz = H, que dz no decrece, que la lista de `cases.yaml` es la que da la receta, y que H es la decena de la profundidad máxima final.
- **Celdas parciales:** hFacMin = 0.2 y hFacMinDr = 2.0 m (ruling R1), en `cases.yaml`. Como la celda de arriba mide 1.5 m, hFacMnSz = 1 en los niveles 1–2 (`model/src/ini_masks_etc.F:108`), así que toda columna mojada tiene al menos dos celdas completas.

## 4. Procedimiento (`bathy.build`)

1. **Reproyección por promedio de área** (`to_model_grid`). Cada punto de IBCAO de 100 m se transforma a UTM con pyproj y se asigna a la celda del modelo que lo contiene. La celda toma la media de sus puntos (unos 4 por celda), con la tierra contada como profundidad 0. El volumen se conserva solo en conjunto: cada celda recibe 3–5 puntos de área UTM algo distinta, así que profundidad × dx² no es exactamente el volumen de mar de esa celda, pero la suma sobre una región coincide con el volumen de los puntos dentro de mucho menos del 1 %.
   - El área de cada punto en UTM sale del jacobiano de la transformación; se usa en los volúmenes de IBCAO.
   - Pruebas sintéticas:
     - una rampa lineal da exactamente la media de los puntos (error < 1e-6);
     - una cuenca gaussiana conserva el volumen dentro del 1 %;
     - un punto en 77.8°N, 15.5°E cae en la celda UTM correcta, en una caja asimétrica para que un eje invertido falle.
2. **Máscara** (`coast_mask` y `make_wet`). Una celda es tierra si más de la mitad de su área es tierra en el polígono del NPI. El área es exacta, por intersección con shapely en bloques de 25 × 25 celdas. Donde IBCAO y el NPI difieren, manda el NPI (regla `npi`). IBCAO cuenta como mar si más de la mitad de sus puntos están bajo el nivel del mar.
   - **Ediciones del PD1 desde `cases.yaml`** (`bathy: force_land / force_sea`). Cada región es un polígono en EPSG:32633, un punto semilla que rellena una región conectada por 4 vecinos, o ambos. `force_sea` puede llevar una condición: dato medido de IBCAO (TID 10–17) y profundidad mayor que 2 niveles (3.0 m). Se aplican después de la máscara del NPI y antes de la limpieza (`apply_force_regions`), y cada región se cuenta en `bathy_report.json` (`force_regions`). Valores actuales:
     - **Vallunden cerrada** (ruling R6). `force_land` con semilla en (543.5, 8648.0 km) y un polígono que toma la laguna al norte de su boca (y ≥ 8646.6 km) y al este de la barra (x ≥ 542.3 km, y ≥ 8644.2 km). Pasan a tierra **453 celdas (18.1 km²)**.
     - **Cabeza de Van Keulenfjorden** (ruling R7). `force_sea` con semilla en (528.4, 8600.0 km), solo con dato medido y más de 3.0 m. Pasan a mar **455 celdas (18.2 km²)**, todas con TID medido (prueba con datos reales).
       - El polígono se detiene en y = 8598 km, 2 km al norte de la esponja sur. El mar de IBCAO sigue hacia el sur fuera de la caja, y sin ese límite la cabeza del fiordo tocaría la frontera abierta sur.
       - Sin el límite y con la misma condición, se habrían convertido 843 de las 1058 celdas del grupo.
3. **Limpieza** (`clean`), que se repite hasta que nada cambia:
   - las celdas mojadas que no se conectan por 4 vecinos (caras de la malla C) con una frontera abierta pasan a tierra;
   - las ensenadas de una celda (tres caras de tierra) se rellenan;
   - las celdas menos profundas que max(`min_depth_m` = 3 m, fondo del nivel 2 = 3.0 m) se profundizan a 3.0 m. Eso es más de tres veces la amplitud de marea fuera de Akseløya, de 0.3–0.8 m (Støylen & Fer 2014, p. 89);
   - las profundidades se ajustan a lo que MITgcm representa con celdas parciales (`ini_masks_etc.F:108-121`), de modo que el archivo es idempotente.
4. **Fronteras abiertas** (`open_boundaries`). En la esponja de 10 celdas de cada frontera abierta, la profundidad (y con ella la máscara) se copia de la primera fila o columna interior, así que ∂h/∂n = 0 exactamente. Después se repite la limpieza: no cambia nada. `build` comprueba al final que la esponja sigue con ∂h/∂n = 0 (`check_sponge`), y hay una prueba con datos reales.
5. **Suavizado** (`smooth_if_needed`): apagado (`max_step_cells: null`). Ver §6.

## 5. Resultados (`figures/bathy_report.json`)

### Controles

| Control | Modelo | Referencia | Tolerancia | Resultado |
|---|---|---|---|---|
| Ancho de Mariasundet, NPI | 1256 m de agua | 700 m (S&F, p. 96) | ±200 m | **FALLA**, clasificada (a) *post hoc* |
| Ancho de Mariasundet, malla | 1329 m de celdas mojadas | 700 m | ±200 m | **FALLA**, clasificada (a) *post hoc* |
| Área de Akselsundet | 3.24 × 10⁴ m² | 2.6 × 10⁴ m² (S&F, p. 93) | ±50 % | pasa (+25 %) |
| Umbral de Akselsundet | 22.4 m | 21.0 m (IBCAO 100 m) | ±20 % | pasa |
| Umbral de Mariasundet | 10.6 m | 11.2 m (IBCAO 100 m) | ±20 % | pasa |
| Profundidad máx. del fiordo interior | 119.0 m | 117 m (S&F, Tabla 1) | ±10 % | pasa |
| Profundidad máx. del fiordo interior | 119.0 m | 119.2 m (IBCAO) | ±5 % | pasa |
| Volumen del fiordo interior | 32.95 km³ | 32.93 km³ (IBCAO) | ±5 % | pasa (+0.06 %) |

El **umbral** se define como la profundidad de cuello de botella entre Bellsund y el fiordo interior con el otro estrecho cerrado: la mayor h tal que las celdas con profundidad ≥ h todavía unen ambos lados. En IBCAO se calcula igual, sobre su malla de 100 m.

**Clasificación de la falla de Mariasundet.** Se escribió después de correr, como etiqueta *post hoc*; el control sigue marcado como fallado en el JSON. La clasificación es **(a), definición de la métrica** (ruling R8). El control mide el agua en el segmento más corto entre Akseløya y la costa sur, que cruza un islote de 0.35 km² y tiene 1.26 km de agua. No mide los dos pasos que el NPI S100 separa con ese islote: **0.50 km** (Akseløya–islote) y **0.48 km** (islote–costa), 0.98 km en total.

- El paso mayor, de 0.50 km, queda 0.20 km por debajo de los 0.7 km de Støylen & Fer, justo en el borde de la tolerancia de ±200 m. Es una comparación *post hoc* y no reemplaza el control.
- La malla reproduce la sección del NPI con 73 m de diferencia, menos de una celda.

### Estrechos (revisión 1 del plan)

| | Akselsundet | Mariasundet |
|---|---|---|
| Sección | Akseløya (489.84, 8630.27 km) → costa norte (489.85, 8631.34) | Akseløya (494.97, 8623.65) → costa sur (496.16, 8622.32) |
| Ancho de agua NPI | 1073 m | 1256 m (dos pasos de 502 y 481 m) |
| Celdas mojadas que cruza la sección | 6 | 10 |
| Profundidad máx. en la sección, malla / IBCAO | 48.8 / 50.5 m | 26.7 / 29.0 m |
| Área de la sección, malla / IBCAO | 3.24 / 3.04 × 10⁴ m² | 1.94 / 2.39 × 10⁴ m² |
| Umbral, malla / IBCAO | 22.4 / 21.0 m | 10.6 / 11.2 m |
| Abierto con ≥ 3 celdas | sí | sí |

"Abierto con ≥ 3 celdas" es una prueba topológica. Con el otro estrecho cerrado, se erosiona el mar una celda con una cruz (caras de la malla C). Si Bellsund y el fiordo interior siguen conectados, el paso mide al menos 3 celdas en todas partes. Los dos estrechos la pasan, y la prueba con datos reales lo verifica.

### Volumen e hipsometría (fiordo interior)

El fiordo interior es el mar del modelo conectado con un punto del fiordo (510, 8632 km) con los dos estrechos cerrados.

**La región de referencia de IBCAO es construida, y eso explica en parte el acuerdo de volumen.** Para IBCAO se toman los puntos de las celdas de mar del fiordo interior más las celdas costeras de tierra a menos de 3 celdas, 600 m (`inner_buffer_cells: 3`), cuyo mar más cercano es el fiordo interior. Con eso entra el mar de IBCAO que el NPI llama tierra a lo largo de la costa, que es justo lo que compensa a la malla. Con otro margen el número cambia; el +0.06 % no es una prueba independiente de conservación. La prueba independiente es la sintética (cuenca gaussiana, < 1 %).

Valores con Vallunden cerrada (R6) y la cabeza de Van Keulenfjorden abierta (R7):

- **Área y volumen:** el área es 526.8 km² en la malla y 519.3 km² en IBCAO (+1.4 %). El volumen sube un 0.06 %, de 32.93 a 32.95 km³. Antes de cerrar Vallunden eran 545.0 frente a 525.0 km² y +0.18 %.
- **Hipsometría por bandas de 10 m:** la diferencia está en la banda de 0–10 m, con **49.5 km² en la malla frente a 41.8 km² en IBCAO** (+18 %). Con Vallunden abierta eran 66.5 frente a 46.2 km². Son las celdas que el NPI llama mar y que IBCAO llama tierra o casi tierra, profundizadas a 3 m (planicies de marea de la cabeza del fiordo y la costa norte). De 10 m hacia abajo las bandas coinciden dentro de ±5 % (por ejemplo, 60.1 frente a 60.3 km² en 80–90 m).
- **Área bajo el nivel de 58.5 m:** 299.6 km² en la malla frente a 303.8 km² en IBCAO (−1.4 %).
- **Volumen bajo 58.5 m:** 8.862 frente a 8.872 km³ (−0.1 %).

### Calidad del dato (TID) en el fiordo

- **Fiordo interior:** el **48.6 %** de las celdas de mar del modelo tiene como TID más frecuente un dato medido (10–17). El **50.0 %** es interpolado (40–47) y el 1.4 % es código de tierra.
- **Toda la caja:** el 74 % es medido: Bellsund y la plataforma tienen sobre todo dato medido.
- **Dónde está lo interpolado** (figura, panel c): en las franjas costeras del fiordo interior y en su mitad este, mientras que el eje central del fiordo es medido.

**Como la mayor parte del fiordo interior es dato interpolado, esto se declara como limitación en el README del Exp. 4.** Ese README todavía no existe; queda anotado para cuando se escriba.

### Limpieza (`clean`), primera pasada

- Se quitaron 2 componentes desconectados (17 celdas), lagunas de profundidad 0 dentro de la barra de arena de la boca de Vallunden.
- Se rellenaron 29 ensenadas de una celda.
- Se profundizaron 2591 celdas a 3.0 m, en promedio 2.2 m cada una.
- Se recortaron 99 celdas al fondo de la malla, todas en la esponja sur, que la copia de la esponja reemplaza.
- 22 181 celdas cambiaron por las celdas parciales: un máximo de 2.09 m y una media de +0.002 m.
- La segunda pasada, después de la esponja, no cambia nada.
- **Esponja:** cambiaron 4144 celdas; 31 pasaron de tierra a mar, donde la costa cruza oblicua el borde N o S, y ninguna de mar a tierra.

### Desacuerdos entre IBCAO y la costa del NPI

En total hay 3320 celdas en desacuerdo, 1017 de ellas en el fiordo interior:

- **IBCAO mar y NPI tierra:** 1420 celdas (192 en el fiordo interior). El grupo mayor tiene **1058 celdas (42 km²) en la cabeza de Van Keulenfjorden** (528.4, 8598.4 km), donde IBCAO tiene mar de hasta **76 m**. Lo más probable es el retroceso de Nathorstbreen después de su *surge* de 2008–2009 [verificar]: la batimetría de IBCAO es más nueva que el frente de glaciar del S100. Con la regla `npi` sola, ese mar se perdía. Con R7 se recuperan 455 celdas (18.2 km²) con dato medido al norte de y = 8598 km (ver §4). El resto es la parte al sur de 8598 km (esponja y su margen) y las celdas interpoladas o de menos de 3 m.
  - Los demás grupos son pequeños: 61 celdas en la costa norte del fiordo (≤ 6 m), 47 en Svea (≤ 6 m) y 24 en la bahía al norte de Akselsundet (≤ 22 m).
- **IBCAO tierra y NPI mar:** 1900 celdas (825 en el fiordo interior). Los grupos mayores son:
  - Vallunden, 349 celdas (IBCAO ≤ 2 m);
  - la costa oeste de la plataforma en (467.9, 8647.2), 274 celdas;
  - la costa norte del fiordo en (529.1, 8635.0), 99 celdas;
  - la cabeza de Recherchefjorden en (493.3, 8600.4), 96 celdas.

  En todos esos lugares IBCAO no tiene profundidad (0–2 m) y el modelo pone 3 m.

### Vallunden

**Cerrada por omisión** (ruling R6, región `force_land: Vallunden`): 453 celdas (18.1 km²) pasan a tierra, y la prueba con datos reales verifica que el punto (543.5, 8648.0 km) es tierra.

Sin esa región no habría quedado como lago. En la costa del NPI, Vallunden se abre al fiordo por un paso al oeste de una barra de arena que a 200 m tiene ≥ 3 celdas, así que `clean` no la quitaba. Además, IBCAO no tiene profundidades en la mayor parte de la laguna (solo el brazo oeste tiene dato, de 15–23 m), y el modelo la habría convertido en una laguna de 3 m. Si el PD1 quiere abrirla, basta borrar la región de `cases.yaml`.

## 6. ¿Hace falta suavizar?

Por omisión, **no**. En coordenadas z con celdas parciales, MITgcm no tiene el error de gradiente de presión de las coordenadas σ, que es lo que motiva limitar rx0. Aquí rx0 llega a 0.85 en las costas empinadas (p99 = 0.44), y eso no es un criterio para z.

El salto de fondo entre vecinos es de hasta 20 niveles (p99 = 7, p99.9 = 13). Está en las paredes de Bellsund y del fiordo y es real: la pendiente costera del fiordo es empinada (S&F, p. 89).

Si el cronometraje del PD1 muestra inestabilidad en el fondo, `smooth_if_needed` limita el salto a `max_step_cells` niveles, solo haciendo menos profundas las celdas que lo exceden, y reporta el cambio máximo y el medio.

## 7. Pendiente para el PD1

- Aprobar la caja de 100 × 60 km o volver a 100 × 45 con Van Keulenfjorden y Recherchefjorden cerrados.
- Confirmar las ediciones de máscara: Vallunden cerrada (R6) y la cabeza de Van Keulenfjorden con dato medido abierta hasta y = 8598 km (R7, 18.2 km²).
- Ancho de Mariasundet: falla clasificada (a), definición de la métrica, *post hoc*. Pasos del NPI: 0.50 y 0.48 km.
- **Área de la sección de Mariasundet:** 1.94 × 10⁴ m² en la malla frente a 2.39 × 10⁴ m² en IBCAO (−19 %). El estrecho pierde área transversal a 200 m, lo que importa para el flujo por los estrechos.
- Aprobar la malla vertical (40 niveles, 1.5 m arriba, H = 240 m) después del cronometraje.
- Agregar al README la limitación de dato interpolado (50 % del fiordo interior, frente a 48.6 % medido).
