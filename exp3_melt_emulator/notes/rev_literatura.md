# Revisión de literatura del Exp. 3: emulador del melt basal

Fecha: 2026-11-10. Se hizo antes de diseñar el experimento y antes de entrenar nada.

**Pregunta del Exp. 3** (spec §5): ¿puede una red pequeña, entrenada con el ensamble del Exp. 2, predecir el mapa de melt mejor que una parametrización cuadrática simple, y cuándo falla al extrapolar a CDW más cálida?

**Fuente.** Diez búsquedas en Consensus, más la búsqueda inicial de reconocimiento (detalle en el registro de búsquedas, al final). En la biblioteca local había pocos artículos del tema: solo Holland 2008, Favier 2019 y Burgard 2022. Los artículos marcados con ⬜ todavía no están en `literatura/`. Lo que dice este documento de ellos sale **solo de su resumen**, y hay que leerlos completos antes de citarlos en el README. Todos están en `readinglist.md`, con su DOI.

## 1. Lo esencial en cinco puntos

1. **El emulador con red neuronal ya existe, con tres variantes.** MELTNET (Rosier et al., 2023) usa una CNN de segmentación entrenada con NEMO en geometrías sintéticas. Burgard et al. (2023) usan un MLP por celda, entrenado con NEMO circumantártico. DeepMelt-GL (Ockenden et al., 2026) es un MLP para cavidades pequeñas y zonas cercanas a la línea de apoyo. Ninguno de estos es "el primero", y el README no lo dirá.
2. **Lo que el spec planeaba como aporte ya está publicado.** El plan era usar la dispersión de un ensamble de redes como indicador de confianza y probar la extrapolación a CDW más cálida. DeepMelt-GL ya muestra las dos cosas: la red es "particularly sensitive to temperatures and ice-draft slopes outside those used in the training dataset", y la dispersión del ensamble de redes "is a good indicator of the reliability of the emulator". El Exp. 3 tiene que **reproducir** ese resultado en un caso controlado y citarlo, no presentarlo como nuevo.
3. **El antecedente más cercano es Bouissou, Burgard & Jourdain (2022).** Entrenan una red profunda en una **geometría idealizada**, con ocean simulations de alta resolución, y concluyen tres cosas: las redes grandes rinden mejor, el formato de entrada de T y S es lo que más importa, y la red puede aplicarse fuera del rango de entrenamiento "if trained appropriately". Es un artículo de congreso (ECCOMAS) y hay que leerlo completo.
4. **Con pocos datos, las líneas base simples a menudo ganan.**
   - En ClimateBench, una regresión lineal tipo *pattern scaling* superó a un modelo de 100 M de parámetros en 3 de 4 variables (Lütjens et al., 2025).
   - En emulación de simuladores, los procesos gaussianos fueron más precisos y estuvieron mejor calibrados que las redes (Myren & Lawrence, 2021).
   - Para la contribución antártica al nivel del mar, en cambio, una LSTM superó a un proceso gaussiano (Van Katwyk et al., 2023).

   Con 26 corridas de entrenamiento, el Exp. 3 tiene que comparar contra una regresión simple y contra un proceso gaussiano, no solo contra la cuadrática.
5. **Las parametrizaciones divergen justo donde el Exp. 3 va a extrapolar.**
   - Con un calentamiento idealizado bajo la termoclina, cinco modelos de melt dan entre +67 % y +240 % de melt total con +1 °C, y entre +141 % y +680 % con +2 °C. Calibrar todos al presente "does not guarantee consistent melt sensitivities" (Lambert & Burgard, 2025).
   - En condiciones mucho más cálidas, tanto las parametrizaciones como la red "struggle" (Burgard et al., 2023).
   - En Favier et al. (2019), la de plumas falla en escenarios cálidos y subestima la contribución al nivel del mar en decenas de porcentaje.

## 2. Orden de lectura sugerido

1. **Burgard et al. (2022)**, *An assessment of basal melting parameterisations for Antarctic ice shelves* (✅ local). Es la mejor visión de conjunto de las líneas base: cuadrática local y no local, plumas y PICO, evaluadas contra NEMO. Leer cómo calibran cada parametrización y qué métricas usan.
2. **Holland et al. (2008)** (✅ local), el fundamento de la ley cuadrática y del H1 del Exp. 2.
3. **Burgard et al. (2023)**, *Emulating present and future simulations of melt rates… with neural networks* (⬜), el diseño de referencia de un MLP por celda contra parametrizaciones. Fijarse en la prueba sobre un ensamble independiente y en el caso "much warmer".
4. **Ockenden et al. (2026), DeepMelt-GL** (⬜). Es lo más reciente y define qué ya no es nuevo. Revisar cómo construyen el ensamble de redes y cómo relacionan la dispersión con el error.
5. **Bouissou, Burgard & Jourdain (2022)** (⬜), el caso idealizado. Revisar qué entradas usan para T/S y qué significa "trained appropriately" para extrapolar.
6. **Lambert & Burgard (2025)** (⬜). Muestra por qué extrapolar es la pregunta abierta: es la incertidumbre que el Exp. 3 puede ilustrar.
7. **Lütjens et al. (2025)** (⬜). Advierte sobre comparar contra líneas base débiles y sobre el sobreajuste al ruido.

## 3. Cómo llegó el campo aquí

| Año | Hito | Por qué importa |
|---|---|---|
| 2003 | Beckmann & Goosse: melt ∝ ΔT con un área efectiva | Primera parametrización simple para modelos de clima |
| 2008 | Holland, Jenkins & Holland: melt cuadrático en ΔT | Base de las parametrizaciones cuadráticas (ISMIP6) |
| 2009 | Little et al.: la eficiencia del flujo de calor depende de la forma y de T | Las simulaciones se apartan de la ley cuadrática simple |
| 2010 | Olbers & Hellmer: modelo de caja | Melt cuadrático para T pequeña, lineal para T grande |
| 2018 | PICO (Reese et al.) y plumas (Lazeroms et al.) | Parametrizaciones con circulación en la cavidad |
| 2019–2020 | Favier et al. (evaluación con NEMO–Elmer/Ice); protocolo ISMIP6 (Jourdain et al., 2020) | Cuadrática local o no local como estándar de proyecciones |
| 2022–2023 | Bouissou et al.; MELTNET (Rosier et al.); Burgard et al. (2023) | Primeros emuladores de red neuronal |
| 2025–2026 | Lambert & Burgard (dispersión entre modelos de melt); DeepMelt-GL (ensamble de redes con indicador de confianza); LADDIE v2 | La extrapolación y la incertidumbre pasan al centro |

**Términos.** Antes de 2015 se habla de *parameterization* y *melt rate sensitivity*; desde 2022, de *emulator*, *neural network*, *multilayer perceptron* y *out-of-distribution*. La forma cuadrática aparece como *local / non-local thermal forcing* (ISMIP6). Hay que buscar con los dos conjuntos de términos.

## 4. Sub-áreas

### 4.1 Emuladores de red neuronal del melt basal

**Qué muestra la literatura.**
- Un MLP por celda o una CNN reproducen bien el melt de NEMO dentro del rango de entrenamiento (**Burgard et al., 2023; Rosier et al., 2023; Ockenden et al., 2026**).
- Fuera de ese rango se degradan (**Ockenden et al., 2026**), aunque en condiciones cálidas quedan más cerca de la referencia que la mayoría de las parametrizaciones (**Burgard et al., 2023**).
- En geometría idealizada, el formato de las entradas de T/S es lo que más pesa (**Bouissou et al., 2022**).

**Artículos clave:**
- Rosier et al. (2023), *The Cryosphere* 17:499, 12 citas. MELTNET, una CNN sobre geometrías sintéticas. Su preprint de 2022 tiene revisiones públicas (T. Smith, J. Bolibar, G. Jouvet) que discuten las debilidades del método.
- Burgard et al. (2023), *JAMES*, 10 citas. MLP por celda y prueba en un ensamble independiente. **Apareció en 4 búsquedas.**
- Ockenden et al. (2026), *GMD* 19:7979. Ensamble de redes con indicador de confianza; sensibilidad fuera de la distribución de entrenamiento.
- Bouissou, Burgard & Jourdain (2022), ECCOMAS. Caso idealizado. **Apareció en 4 búsquedas.**

**Términos:** basal melt emulator, neural network parameterisation, multilayer perceptron, MELTNET, DeepMelt, sub-shelf melt, ice-shelf cavity, deep learning.

```
("basal melt" OR "sub-shelf melt" OR "ice shelf melt") AND ("neural network" OR "deep learning" OR emulator)
("melt rate" AND "ice shelf") AND ("multilayer perceptron" OR "convolutional") AND (idealised OR synthetic geometry)
```

### 4.2 Parametrizaciones de referencia (lo que el emulador debe superar)

**Qué muestra la literatura.**
- La cuadrática de ISMIP6, local o no local, con dos calibraciones (MeanAnt y PIGL), es el estándar (**Jourdain et al., 2020**).
- Las de plumas van bien en frío y fallan en cálido. PICO con cinco cajas "compares fairly well" en casi todos los escenarios. Entre las escalas simples, hace falta una dependencia cuadrática y no lineal, y la cuadrática mejora si combina T local y no local (**Favier et al., 2019**).
- Desde 2009 se sabe que la respuesta simulada se aparta de la ley cuadrática simple, porque la eficiencia del flujo de calor depende de la forma del hielo y de T (**Little et al., 2009**).
- Calibradas igual al presente, las parametrizaciones discrepan mucho en su sensibilidad (**Lambert & Burgard, 2025**).

**Artículos clave:**
- Favier et al. (2019), *GMD* 12:2255, 138 citas (✅ local). Evalúa parametrizaciones contra un modelo acoplado en una geometría idealizada tipo Pine Island. **Apareció en 4 búsquedas.**
- Jourdain et al. (2020), *The Cryosphere* 14:3111, 160 citas. La fórmula cuadrática exacta que sirve de línea base.
- Burgard et al. (2022) (✅ local). La evaluación sistemática.
- Lambert & Burgard (2025), *The Cryosphere* 19:2495. Dispersión de la sensibilidad entre modelos de melt.
- Reese et al. (2018) (PICO), Lazeroms et al. (2018) (plumas) y Lambert et al. (2023) (LADDIE). Líneas base más complejas, opcionales para el Exp. 3.

**Términos:** quadratic melt parameterisation, local/non-local thermal forcing, ISMIP6 melt, PICO box model, plume parameterisation, LADDIE, melt sensitivity to ocean warming.

```
("basal melt parameterisation" OR "melt parameterization") AND (quadratic OR "thermal forcing") AND (evaluation OR assessment)
("ice shelf" AND "basal melt") AND ("plume model" OR PICO OR "box model") AND ("ocean warming" OR sensitivity)
```

### 4.3 Extrapolación e incertidumbre

**Qué muestra la literatura.**
- Los ensambles de redes (*deep ensembles*) son una forma simple de estimar la incertidumbre (**Lakshminarayanan et al., 2017**).
- Pero pueden quedar sobreconfiados fuera de la distribución (de Mathelin et al., 2023; preprint).
- En emulación de melt, la dispersión entre redes sí indicó la fiabilidad (**Ockenden et al., 2026**).
- En downscaling climático, la extrapolación es una limitación conocida del aprendizaje profundo (Hernanz et al., 2023).

**Artículos clave:**
- Lakshminarayanan, Pritzel & Blundell (2017), 9103 citas. El método.
- Ockenden et al. (2026). La aplicación al melt.
- Gawlikowski et al. (2021), *AI Review*, 2128 citas. Revisión de la incertidumbre en redes.

**Términos:** deep ensembles, predictive uncertainty, out-of-distribution, extrapolation, epistemic uncertainty, ensemble spread, calibration of uncertainty.

```
("deep ensemble" OR "ensemble of neural networks") AND ("out-of-distribution" OR extrapolation) AND (emulator OR surrogate)
("neural network emulator") AND ("uncertainty") AND (climate OR ocean OR "ice sheet")
```

### 4.4 Pocos datos: procesos gaussianos, regresiones y redes

**Qué muestra la literatura.**
- No hay un ganador universal. Con simuladores costosos y pocos datos, el proceso gaussiano fue más preciso y estuvo mejor calibrado (**Myren & Lawrence, 2021**). En ISMIP6, una LSTM superó al proceso gaussiano (**Van Katwyk et al., 2023**).
- Un emulador lineal le ganó a un modelo profundo en 3 de 4 variables, porque las redes sobreajustan el ruido de la variabilidad interna (**Lütjens et al., 2025**).
- Para un ensamble perturbado de melt, Berdahl et al. (2020) usaron un emulador estadístico.

**Artículos clave:**
- Lütjens et al. (2025), *JAMES*. Las líneas base simples ganan.
- Myren & Lawrence (2021). Proceso gaussiano contra red con pocos datos.
- Van Katwyk et al. (2023), *JAMES*. El resultado contrario, en hielo antártico.
- Berdahl et al. (2020), *The Cryosphere* (preprint tc-2020-178). Emulación de un ensamble de melt.

**Términos:** Gaussian process emulator, computer model emulation, perturbed parameter ensemble, small training set, pattern scaling, linear emulator, benchmark.

```
("Gaussian process" OR "kriging") AND ("neural network") AND (emulator OR surrogate) AND ("small" OR "limited training")
("emulator" AND benchmark) AND ("linear regression" OR "pattern scaling") AND ("deep learning")
```

### 4.5 Validación y restricciones físicas

**Qué muestra la literatura.**
- La validación usual es probar con corridas no vistas o con un ensamble independiente (**Burgard et al., 2023**).
- Los *benchmarks* deben usar promedios de ensamble para no premiar el sobreajuste al ruido (**Lütjens et al., 2025**).
- En parametrizaciones de submalla del océano en un modelo idealizado, las redes con restricciones físicas fueron estables en línea pero "generalize poorly to new regimes"; una regresión lineal con programación genética rindió igual en el dominio de entrenamiento (Ross et al., 2023).
- No encontré trabajos de restricciones físicas (por ejemplo, forzar la dependencia del melt con u*·ΔT) aplicadas a emuladores del melt basal. Esta búsqueda devolvió sobre todo flujos turbulentos aire–hielo (Cummins et al., 2023), así que la cobertura de este punto es incompleta.

**Términos:** leave-one-out, cross-validation, held-out simulations, physics-informed, physically constrained, benchmark, online vs offline.

```
("emulator" OR "neural network") AND ("leave-one-out" OR "cross-validation") AND ("ocean model" OR "ice shelf")
("physics-informed" OR "physically constrained") AND ("basal melt" OR "ice-ocean boundary layer")
```

## 5. Grupos que dominan el tema

- **Burgard, Jourdain y Mathiot** (IGE Grenoble / LOCEAN). Aparecen en las sub-áreas 4.1 y 4.2 y en la mayoría de las búsquedas: Bouissou et al. 2022, Burgard et al. 2022 y 2023, Jourdain et al. 2020, Favier et al. 2019, Ockenden et al. 2026, Lambert & Burgard 2025. Están en tu tracker (#21–22).
- **Rosier y Gudmundsson** (Northumbria), MELTNET (tracker #15). **De Rydt**, en el mismo grupo, aparece por la amplificación geométrica (De Rydt & Naughten, 2024).
- **Lambert** (IMAU/KNMI), LADDIE v1 y v2, y la comparación de sensibilidades con Burgard (tracker #30).
- **Lütjens y Ferrari** (MIT), el *benchmarking* de emuladores contra líneas base lineales. No está en el tracker; es relevante para MIT-WHOI.

## 6. Huecos y lo que significan para el Exp. 3

**Lo que el Exp. 3 puede aportar sin decir que es el primero:**

1. **Un caso donde la ley cuadrática falla de forma conocida.** El Exp. 2 mostró que el melt es casi exactamente cuadrático en Tbot, pero con exponente 2.4–2.8, y que depende de z_c. La cuadrática no local de ISMIP6 usa una sola T lejana, y la local usa la T a la profundidad del hielo. El Exp. 3 puede medir cuánto de la dependencia en z_c captura cada una, porque la local "ve" la termoclina y la no local no. Esto es una pregunta física con respuesta medible, y no encontré que nadie la haya probado en un ensamble controlado.
2. **Una comparación honesta con 26 corridas.** Hay que comparar el ensamble de redes contra la cuadrática local, la cuadrática no local, una regresión lineal sobre los parámetros (Tbot, z_c, s) y un proceso gaussiano. La literatura no da un ganador para este tamaño de datos (Myren & Lawrence 2021 frente a Van Katwyk et al. 2023), y Lütjens et al. (2025) muestran que omitir la línea base simple puede inflar el mérito de la red.
3. **La extrapolación a Tbot = 2.0 °C como réplica controlada de DeepMelt-GL.** El objetivo es comprobar si la dispersión del ensamble de redes crece donde crece el error en las tres corridas de prueba, y citar a Ockenden et al. (2026) como el resultado que se reproduce.

**Huecos metodológicos:**
- Ningún emulador revisado se entrena con un ensamble tan pequeño (decenas de corridas). MELTNET usa muchas geometrías sintéticas y Burgard et al. (2023), NEMO circumantártico. No sabemos si sus conclusiones valen con 26 corridas.
- Las restricciones físicas en emuladores del melt basal están casi ausentes de los resultados (sub-área 4.5).

**Huecos conceptuales:**
- La forma de la respuesta del melt al calentamiento sigue abierta: cuadrática (Holland et al., 2008), cuadrática a lineal (Olbers & Hellmer, 2010) o apartada de la cuadrática por la eficiencia (Little et al., 2009). El Exp. 2 encontró una respuesta más empinada que la cuadrática. Un emulador que no supone ninguna forma puede ayudar a describirla, pero no a explicarla.

**Riesgos de diseño (para el plan):**
- Con 26 mapas de 21 × 120 celdas, una red por celda tiene muchas muestras, pero solo 26 condiciones de forzamiento independientes. El error de generalización depende de esas 26, no del número de celdas. La validación tiene que dejar fuera **corridas** completas, no celdas.
- Hay 5 corridas que fallan N1 y están marcadas en el dataset. Hay que decidir de antemano si entrenan, y reportarlo.

## 7. Registro de búsquedas

Consensus devolvió hasta 20 artículos por búsqueda. El plan permite hasta 300 por búsqueda con `page_size`, que no se usó.

| # | Búsqueda | Filtros | Artículos | Estado |
|---|---|---|---|---|
| 1 | machine learning emulator of ice shelf basal melt rates compared with melt parameterizations | — | 20 | OK |
| 2 | neural network emulation of sub-ice-shelf melt rates from ocean temperature and cavity geometry | — | 20 | OK |
| 3 | basal melt parameterisations quadratic thermal forcing plume model PICO evaluation against ocean model | — | 20 | OK |
| 4 | out-of-distribution extrapolation and uncertainty of neural network emulators deep ensembles climate model | — | 20 | OK |
| 5 | Gaussian process emulator small ensemble ocean or ice sheet model comparison with neural network | — | 20 | OK |
| 6 | physics-constrained machine learning parameterization of ice-ocean melt boundary layer heat flux | — | 19 | OK (fuera de tema: casi todo flujos aire–hielo) |
| 7 | review of machine learning in glaciology and ice sheet ocean interaction | — | 20 | OK (casi todo teledetección) |
| 8 | review ice shelf basal melt parameterizations for ice sheet models sensitivity to ocean warming | — | 20 | OK |
| 9 | parameterization of ice shelf basal melting as function of far-field ocean temperature | year_max 2015 | 20 | OK |
| 10 | (igual que 9) | year_min 2021 | 18 | OK |
| 11 | evaluating melt parameterizations and emulators against cavity-resolving ocean simulations non-local quadratic sensitivity | year_min 2020 | 20 | OK |

- **Búsquedas:** 11 ejecutadas: la inicial más 10 del plan estándar (5 por sub-área, 2 de revisiones, 2 por época y 1 de seguimiento). 11 con resultados y 0 fallidas.
- **Artículos únicos recibidos:** unos 153 (134 en las búsquedas 1–8, 10 y 11, más 19 nuevos en la 9).
- **Citados aquí:** 30.
- **Cobertura débil:** la sub-área 4.5 (restricciones físicas) y la búsqueda de revisiones de ML en glaciología devolvieron poco material pertinente. Si el plan necesita restricciones físicas, hay que buscar a mano.
- **Biblioteca local:** antes de Consensus solo estaban Holland 2008, Favier 2019 y Burgard 2022. `litsearch.py` no se corrió por separado; se revisó la lista de `literatura/`.

## Para descargar (prioridad)

Burgard et al. 2023 · Ockenden et al. 2026 · Bouissou et al. 2022 · Jourdain et al. 2020 · Lambert & Burgard 2025 · Lütjens et al. 2025 · Rosier et al. 2023 · Little et al. 2009 · Myren & Lawrence 2021. Los DOI están en `readinglist.md`, sección Exp. 3.

## 8. Lectura completa (2026-11-10, después de la descarga)

Esteban descargó 21 artículos, que están en `literatura/` con su versión `.md`. Falta Bouissou et al. 2022. Los grupos clave se leyeron completos. Lo que cambia o precisa respecto a las secciones 1–7:

**Emuladores (Burgard 2023, Ockenden 2026, Rosier 2023).**
- **Burgard 2023.** Usa 4 corridas NEMO circumantárticas con la geometría fija; no son experimentos de forzamiento controlados. El MLP por celda tiene 14 entradas, y la T/S se toma del perfil a 50 km del frente, llevado a la profundidad local. Hay un ensamble de 10 redes, pero no se analiza su dispersión.
  - En la validación cruzada por plataforma, RMSE_int baja 31–61 % frente a las parametrizaciones.
  - En 4xCO2, RMSE_GL ≈ la señal misma (sin destreza), y la red da ~3× menos melt en Pine Island.
  - Por importancia de variables, la posición y la altura de la columna de agua pesan más que T dentro del régimen.
- **Ockenden 2026.** Un MLP de 32/64/32 con 10 semillas. No compara contra ninguna parametrización.
  - Extrapolar en T falla: de ANT2000 a ANT2300, RMSE 2220 Gt/a.
  - "We do not see any examples of a low ensemble spread and high RMSE"; el umbral provisional es una dispersión mayor que ~10 % del melt integrado.
  - No analiza la termoclina.
- **Rosier 2023.** Unas 11 000 corridas NEMO idealizadas (ISOMIP+ con geometrías aleatorias), con perfiles de T/S generados por una GAN a partir de WOA. No mide incertidumbre ni prueba extrapolación.
  - En una prueba de una losa, el melt es cuadrático en T uniforme, sin exponente reportado.
  - La termoclina se menciona como motivación, pero no se prueba.
- **Lo que ninguno hace:** variar Tbot, z_c y la geometría de forma controlada, uno a la vez; probar la profundidad de la termoclina; imponer restricciones físicas (las tres lo proponen como trabajo futuro).

**Parametrizaciones (Favier 2019, Jourdain 2020, Burgard 2022, Lambert 2025, Little 2009, Olbers 2010, Reese 2018, Lazeroms 2018).**
- **Fórmulas implementables con solo el calado y T(z), S(z):**
  - lineal local (Favier, Ec. 2);
  - cuadrática local (Favier, Ec. 4; Jourdain, Ec. 2; Burgard, Ec. 14, con constante K);
  - no local / M+ (Favier, Ec. 5; Jourdain, Ec. 1) y semilocal (Burgard, Ec. 18).

  Las tres tienen un solo parámetro libre. PICO tiene dos (γ_T*, C) y complejidad moderada. La pluma de Lazeroms es la más difícil, por el origen efectivo de la pluma.
- **Elegir qué T usar es elegir la física.** Con T(z_draft), la parametrización "ve" la termoclina. Con T a profundidad fija (500 o 700 m), supone agua que entra por el fondo. Favier (Warm3) encontró que subir la termoclina 200 m no cambia el melt de NEMO acoplado, y que la versión de profundidad fija lo reproduce y la local no.
- **Forma de la respuesta al calentamiento.** Ningún artículo predice un exponente > 2.
  - Little 2009 da exponentes de 1.55–1.81 y desviaciones negativas respecto a la cuadrática, porque la eficiencia hielo–océano baja con T.
  - Olbers & Hellmer 2010 y PICO dan cuadrática a lineal.
  - En Lambert & Burgard 2025, la cuadrática y LADDIE son los más sensibles, con amplificación profunda de 6.7 y 4.2 contra 3.3 observada; PICO da 2.1 y la pluma ~1.
- **Inferencia propia, a verificar en el Exp. 3.** Un polinomio cuadrático en Tbot que no se anula en T_f da necesariamente un exponente ≠ 2 respecto a TF. Hay que comparar ajustes contra Tbot y contra TF antes de interpretar el exponente de 2.4–2.8 del H1.
- **Métricas recomendadas:**
  - RMSE del melt integrado (Burgard, Ec. 33);
  - RMSE cerca de la línea de apoyo, en la primera caja (~10 % del área; Ec. 35);
  - sensibilidad dM/dT, no linealidad y amplificación profunda (Lambert, Ec. 1–3);
  - patrón normalizado por el percentil 90.

**Emulación con pocos datos (Myren 2021, Lütjens 2025, Van Katwyk 2023, Berdahl 2021, Lakshminarayanan 2017, Ross 2023).**
- **Myren 2021.** El caso de cosmología (37 corridas, 5 entradas) es el más parecido al nuestro. Un proceso gaussiano sobre los pesos de PCA (99 % de la varianza), con validación dejando una fuera, gana en todos los casos. Las redes requirieron probar "hundreds of architectures".
- **Lütjens 2025.** La línea base lineal por celda gana por el compromiso sesgo–varianza: la red sobreajusta el ruido. Recomiendan regularizar con supuestos físicos sobre la forma funcional, y advierten que la variación entre semillas puede invertir el orden de los métodos.
- **Van Katwyk 2023.** La red gana con muchos datos, pero "GPs remain a good choice … for small data sets".
- **Lakshminarayanan 2017.** M = 5 redes entrenadas con la verosimilitud negativa (NLL), con inicialización aleatoria y sin *bagging*. La evidencia fuera de distribución en regresión es solo un ejemplo 1-D de juguete.
- **Ross 2023.** Las redes "generalize poorly to new regimes", y un mejor R² fuera de línea llegó a asociarse con un peor desempeño en línea. Recomiendan "better physics baselines".
- **Protocolo sugerido:**
  - validación dejando fuera una corrida completa, con todo el ajuste (estandarización, número de EOF, núcleo) dentro de cada pliegue;
  - las 3 corridas a Tbot = 2.0 se tocan una sola vez;
  - pliegues adicionales que dejan fuera un z_c y que dejan fuera las corridas más cálidas del entrenamiento;
  - cobertura empírica de los intervalos al 68 % y al 95 %, y CRPS;
  - errores por corrida con intervalos *bootstrap*;
  - conocer el piso de ruido entre corridas.
