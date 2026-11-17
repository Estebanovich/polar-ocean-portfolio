# Revisión de literatura del Exp. 5: sensibilidad adjunta del transporte de calor sobre el talud

Fecha: 2026-11-14. Se hizo antes de diseñar el experimento y antes de correr nada.

**Pregunta del Exp. 5** (spec del portafolio, §Exp. 5): en un canal reentrante con talud, ¿a qué campos de viento y de batimetría es más sensible el transporte de calor hacia la plataforma? Se usaría el adjunto de MITgcm generado con Tapenade, como en `canyon_mx`.

**Fuentes.**
- **Consensus:** una búsqueda inicial de reconocimiento y diez búsquedas dirigidas (detalle en el registro, al final). Cada búsqueda devolvió entre 10 y 18 artículos, el número que pedí: la cuenta no limita la cobertura.
- **Biblioteca de la tesis** (`~/Documents/TESIS_InternalWaves-2/Literatura/`): tiene 3 de los artículos de aquí, marcados 📁: Heimbach et al. (2005), Gaikwad et al. (2024) y Tber et al. (2007). Varios artículos de esa biblioteca citan a Heimbach & Losch (2012), pero el artículo mismo no está.
- **Biblioteca del portafolio** (`literatura/`): no tiene nada de adjuntos ni del talud antártico.

Los marcados ⬜ no están en ninguna de las dos bibliotecas. Lo que dice este documento de ellos sale **solo de su resumen**, y hay que leerlos completos antes de citarlos en el README. Todos están en `readinglist.md`, con su DOI.

## 1. Lo esencial en seis puntos

1. **No encontré ningún estudio adjunto del transporte de calor a través del talud antártico.**
   - En las 11 búsquedas, todas las sensibilidades de ese transporte se obtienen con experimentos de perturbación: se escala el viento o se cambia la geometría y se compara (Dinniman et al., 2012; Stewart & Thompson, 2012; Daae et al., 2017; Bull et al., 2021, 2025; Morrison et al., 2023).
   - Los adjuntos polares que sí aparecen miden otra cosa:
     - el derretimiento bajo Pine Island (Heimbach & Losch, 2012);
     - el transporte por el paso de Drake respecto al viento (Mazloff, 2011);
     - la exportación de hielo por el Archipiélago Ártico Canadiense (Heimbach et al., 2010);
     - el contenido de calor del Pacífico sureste (Jones et al., 2019, 2020);
     - la estratificación del agua modal (Pimm et al., 2024, 2025);
     - la estimación de estado del mar de Amundsen (Nakayama et al., 2021).
   - La búsqueda por época desde 2021 tampoco devolvió ninguno. Ese es el hueco. Hay que decirlo con cuidado: es lo que devolvió Consensus, no una revisión exhaustiva.
2. **El transporte de calor hacia la plataforma lo llevan los remolinos, y eso choca con el adjunto.**
   - Sobre el talud, el flujo de calor neto hacia la costa es casi todo de remolinos que mezclan a lo largo de las isopicnas. La marea domina en el quiebre de la plataforma, pero la compensa el flujo medio (Stewart & Thompson, 2014; Stewart et al., 2018; Palóczy et al., 2018).
   - Resolver esos remolinos exige mallas de ~1 km (Stewart & Thompson, 2015, citado por Bull et al., 2025; no lo devolvió Consensus).
   - En un modelo con remolinos, el adjunto crece exponencialmente más allá del tiempo de predictibilidad (Lea et al., 2000; Köhl, 1999).
   - **El precedente directo es McLay (2006):** un adjunto de MITgcm en un **canal reentrante con remolinos**, el mismo marco del Exp. 5.
     - El tiempo no lineal del canal es de unos 200 días.
     - La estructura espacial de la sensibilidad sigue siendo útil a 278 días.
     - Las diferencias finitas dejan de servir como prueba de gradiente, porque también las afecta el caos.
     - Con una parametrización que reduce la energía de los remolinos, el adjunto coincide con las perturbaciones a 690 días.
3. **Hay una salida conocida: parametrizar los remolinos con un GM sensible al talud.**
   - Dettling et al. (2024, 2025; Nakayama es coautor del segundo) muestran que en un modelo grueso un GM/Redi reproduce el flujo de calor hacia la plataforma del modelo fino. La condición es bajar el coeficiente de GM un orden de magnitud sobre el talud.
   - Bull et al. (2025) encuentran que, en escalas de equilibrio, la respuesta de la temperatura de la plataforma al viento es "de primer orden lineal". Eso es justo lo que el adjunto necesita para ser válido.
   - Esto define la decisión central del diseño (§6).
4. **La herramienta existe y es abierta.**
   - Desde 2024, MITgcm tiene un adjunto con Tapenade en la rama principal, MITgcm-AD v2 (Gaikwad et al., 2024 📁).
   - Su Tabla 2 lista los paquetes para los que Tapenade da derivadas correctas. Entre ellos están `gmredi`, `kpp`, `ggl90`, `seaice`, `shelfice`, `exf` y `ptracers`. `obcs` no aparece en la tabla (lo leí del `.md` de la biblioteca, cuyo texto de tablas está desordenado; hay que confirmarlo en el PDF).
   - El criterio estándar de corrección es que el gradiente adjunto coincida con diferencias finitas y con el tangente lineal dentro del **1 %**, con una prueba previa de linealidad: perturbaciones ±Δu que den respuestas opuestas.
   - `tutorial_reentrant_channel` **no** tiene carpeta `code_tap` ni `code_ad` en el MITgcm local. Habría que construir el adjunto del canal desde cero, como se hizo en `canyon_mx`.
5. **El grupo de Mazloff tiene dos piezas que el Exp. 5 puede usar directamente.**
   - **Mazloff (2011):** con el adjunto de SOSE, el transporte por Drake responde al viento en escalas de ~100 días. En zonas de topografía compleja la sensibilidad se concentra en escalas de ~100 km, y ahí manda el rotacional del esfuerzo.
   - **Rodriguez, Mazloff & Gille (2016):** en SOSE, el transporte de calor hacia el mar de Amundsen lo sostiene, en el balance de vorticidad, el **rotacional del viento sobre el quiebre de la plataforma**. Es una hipótesis que un adjunto puede poner a prueba de forma directa.
   - **Lanham, Mazloff et al. (2025):** clasifican cuatro régimenes estacionales de intrusión de CDW. En el régimen 1 (Antártida oriental), los vientos del este controlan la intrusión vía Ekman.
   - **Gopalakrishnan, Cornuelle & Hoteit (2013):** un adjunto de MITgcm en el golfo de México, hecho en Scripps. Es un puente natural con tu experiencia mexicana.
6. **La topografía es un control tan importante como el viento, y el adjunto lo puede medir.**
   - Lösch & Heimbach (2007) calculan con el adjunto la sensibilidad a la batimetría del transporte por Drake y del transporte de calor a 32°S. Es tan grande como la sensibilidad a las condiciones de superficie, y se interpreta como esfuerzo de forma en el fondo.
   - En modelos de proceso, la geometría (cañones, ancho de la plataforma, embahiamientos) cambia el resultado tanto como el viento (Hogg & Munday, 2014; Daae et al., 2017; Ong et al., 2024; Bull et al., 2025).

## 2. Mecanismo físico: el transporte de CDW a través del talud

- **Revisión de partida:** Thompson et al. (2018), *The Antarctic Slope Current in a Changing Climate* (405 citas). Clasifica el frente del talud en plataforma fresca, densa y cálida. Muestra que los remolinos, la marea y los flujos densos son centrales en el intercambio, y que su evolución está mal resuelta.
- **Remolinos.**
  - En un modelo de proceso, el CDW cruza el quiebre solo por remolinos, mientras que el agua de superficie cruza por Ekman. El flujo de CDW es sensible a todo el forzamiento y a la geometría (Stewart & Thompson, 2014). Esa frase es casi una invitación al adjunto.
  - El agua densa que baja por el talud crea remolinos y chorros, y conecta las isopicnas de la plataforma con el CDW (Stewart & Thompson, 2016).
  - Observaciones con focas en la Antártida oriental: el transporte de CDW por remolinos se concentra en un tercio del talud (Foppert et al., 2019).
- **Viento.**
  - La celda profunda de vuelco es muy sensible a los vientos del este sobre el talud y casi insensible a los del oeste (Stewart & Thompson, 2012).
  - Más viento del este baja las isopicnas costeras y enfría la plataforma. El viento meridional, vía el hielo marino, calienta en profundidad (Morrison et al., 2023).
  - En la península, más viento lleva más CDW a la plataforma, pero no siempre más calor a las plataformas de hielo, porque también aumenta la mezcla (Dinniman et al., 2012).
- **Marea y flujo medio:** a 1/48°, la marea domina en el quiebre, pero la compensa el flujo medio. Neto: remolinos (Stewart et al., 2018, 2019).
- **Salinidad y deshielo:** el gradiente de salinidad a través del frente controla el flujo de calor. Más agua de deshielo da más flujo hacia la costa (Si et al., 2023). El deshielo también impulsa la contracorriente del talud (Si et al., 2024).
- **Topografía:**
  - Las intrusiones por cañones son episódicas incluso con forzamiento constante (Ong et al., 2024).
  - En una fosa, el torque de presión de fondo impulsa el flujo medio hacia la costa (Liu et al., 2022).
  - El efecto JEBAR también puede impulsar flujo hacia la costa en el mar de Ross (Wang et al., 2024).
  - Los régimenes de "flujo pasante" y de "giro" cambian cuánto responde la plataforma a un cambio de viento (Bull et al., 2025).
- **Modelos de canal reentrante:** son el marco estándar (Abernathey et al., 2011, la base de `tutorial_reentrant_channel` [verificar]; Morrison & Hogg, 2012; Munday & Zhai, 2015; Matsuta et al., 2026).

## 3. El adjunto de MITgcm

- **Origen:** Marotzke et al. (1999) construyeron el adjunto con TAMC. Con él calcularon la sensibilidad del transporte de calor a 29°N en el Atlántico a la hidrografía inicial, y separaron la sensibilidad "cinemática" de la "dinámica". El adjunto coincidió muy bien con perturbaciones directas.
- **Reglas y técnica:**
  - Giering & Kaminski (1998), *Recipes for adjoint code construction* (957 citas), es la base de TAMC/TAF.
  - Heimbach et al. (2005 📁) describen el *checkpointing*, las rutinas adjuntas escritas a mano para MPI y las pruebas del gradiente.
  - Vlasenko et al. (2015) comparan TAF, Tapenade, OpenAD y NAGWare. Las herramientas de transformación de código fuente son órdenes de magnitud más rápidas.
- **Tapenade:** Gaikwad et al. (2024 📁), MITgcm-AD v2; Tber et al. (2007 📁), el mismo método aplicado a OPA.
- **ECCO:**
  - Forget et al. (2015) describen ECCO v4.
  - Una práctica común es **subir la viscosidad para estabilizar el adjunto** y omitir en él algunas dependencias de las parametrizaciones (GGL90 y las pendientes isopicnas de GM).
  - Esa práctica es una advertencia para el Exp. 5: el adjunto que se usa en la práctica no es exactamente el del modelo directo.
- **Aplicaciones regionales:** golfo de México (Gopalakrishnan et al., 2013), mar Rojo (Zhan et al., 2018), ola de calor del Pacífico noreste (Wang et al., 2025), mar del Sur de China (Na, 2025).
- **Con ondas internas:** Shimizu (2024) muestra con MITgcm que el adjunto es estable bajo ondas internas no lineales. Advierte que la **definición de la función de costo** controla la direccionalidad y el contenido modal de la señal adjunta. Toca a la vez tu tesis y `canyon_mx`.

## 4. Adjuntos sobre la plataforma antártica y las plataformas de hielo

- **Heimbach & Losch (2012):** el derretimiento bajo Pine Island conecta con el interior en 30–60 días por advección. Las rutas de mayor sensibilidad pasan sobre un umbral, y hay "regiones de sombra" de sensibilidad muy baja. Es el antecedente más cercano en espíritu.
- **Nakayama et al. (2021):** primera estimación de estado con adjunto en el mar de Amundsen y Bellingshausen con cavidades. El costo baja 65 % en 20 iteraciones. Una advertencia del artículo: la resolución de LLC270 no resuelve los remolinos, que hay que parametrizar.
- **Sin adjunto pero con la misma pregunta:**
  - Bull et al. (2021): 28 perturbaciones de la corriente del talud; respuesta casi lineal a la salinidad y cuadrática a la temperatura.
  - Haigh et al. (2026): la variabilidad decadal de la contracorriente la controlan los flujos de agua dulce de las polinias.
  - Rieke et al. (2025): el contenido de calor en el Amundsen oriental lo controla el viento.

## 5. Límites del adjunto con remolinos

- **Lea, Allen & Haine (2000)**, en el sistema de Lorenz:
  - Sobre trayectorias largas, el adjunto diverge de la sensibilidad del clima.
  - Sobre trayectorias cortas, subestima la respuesta.
  - Existe una ventana intermedia en la que un **ensamble de adjuntos** da ~10 % de error.
- **Köhl (1999):** propone asimilar estadísticas en lugar de trayectorias, para ir más allá del límite de predictibilidad.
- **McLay (2006):** ver §1, punto 2. Es el único precedente en un canal reentrante de MITgcm; es una tesis y no tiene DOI.
- **Gebbie & Hsieh (2017):** para estimar estados, el criterio es la controlabilidad, no el caos. Esto vale para la estimación, no necesariamente para la sensibilidad.
- **Fuera de la oceanografía:**
  - Wang & Zaki (2026) cierran el adjunto promediado en ensamble con una viscosidad turbulenta.
  - Hay métodos de "sombra" para sensibilidad en sistemas caóticos (Thakur & Nadarajah, 2025, preprint).
  - Son ideas, no herramientas listas para MITgcm.
- **Ventanas usadas en la práctica:** de 2–7 días en análisis regionales con remolinos, a meses en cuencas (Köhl et al., 2026, citando a Gopalakrishnan et al., 2013, y a Verdy et al., 2023).

## 6. Lo que esto implica para el diseño (para hablarlo, no decidido)

1. **La decisión central es cómo tratar los remolinos.** Veo tres caminos; la elección es tuya:
   - **(a) Canal grueso con GM sensible al talud.** El adjunto es estable a escalas de meses o años (McLay, 2006; Bull et al., 2025). Cabe en la Mac (≤ 5 núcleos). La física de los remolinos queda parametrizada, y hay que decirlo así.
   - **(b) Canal que resuelve los remolinos con ventana corta**, más corta que el tiempo no lineal. La física es más fiel, pero la sensibilidad es de corto plazo. En la Mac solo es posible a una resolución intermedia y en un dominio pequeño.
   - **(c) Los dos, con una prueba cruzada.** El adjunto se calcula en el modelo grueso, y unas pocas perturbaciones directas en el fino comprueban si las sensibilidades sobreviven a los remolinos. Este es el estilo de Lea et al. (2000) y de Jones et al. (2019).
2. **Función de costo:** el transporte de calor a través de una isóbata (Palóczy et al., 2018, a 1000 m; Stewart et al., 2018), promediado en el tiempo. Shimizu (2024) y Marotzke et al. (1999) muestran que su definición importa.
3. **Controles:** el viento (τx, τy), la batimetría (Lösch & Heimbach, 2007), la hidrografía inicial y κᵥ. Coinciden con los mapas que ya sacó `canyon_mx`.
4. **Pruebas preregistradas:** linealidad con ±Δu; gradiente adjunto contra el tangente lineal y contra diferencias finitas dentro del 1 % (Gaikwad et al., 2024); perturbaciones finitas guiadas por el adjunto, con un tope de error fijado antes (Lösch & Heimbach, 2007, obtienen < 30 %).
5. **La hipótesis que conecta con Mazloff:** ¿el transporte de calor hacia la plataforma es más sensible al *rotacional* del viento sobre el quiebre que a su magnitud (Rodriguez et al., 2016; Mazloff, 2011)? Es una pregunta concreta, y el adjunto la contesta en una sola corrida.

## 7. Qué mejorar de `canyon_mx` según esta revisión

Sin haber vuelto a leer su código, solo contrastando su README con la literatura:
- **Prueba de gradiente.** Hay que confirmar si `grdchk` se hizo dentro del régimen lineal (prueba ±Δu) y con qué tolerancia. En un modelo no hidrostático con turbulencia, las diferencias finitas pueden no servir (McLay, 2006).
- **Ventana del costo.** Son 3 días (días 12.5–15.5) después de 11 días de corrida directa. Hay que mostrar que la sensibilidad no está dominada por crecimiento caótico, por ejemplo con un ensamble pequeño (Lea et al., 2000).
- **Definición de ω.** Shimizu (2024) advierte que, con ondas internas, la función de costo condiciona la señal adjunta. Hay que justificar ω en esos términos.
- **Perturbaciones guiadas por el adjunto.** El README dice "cross-check vs FD". Falta una perturbación finita grande, guiada por el mapa de sensibilidad, que muestre que el mapa predice la respuesta no lineal.

## 8. Quién domina el tema

| Grupo | Temas | Artículo representativo |
|---|---|---|
| Stewart (UCLA), con Thompson (Caltech) y Si | Modelos de proceso del talud, remolinos, marea y deshielo | Stewart & Thompson (2014); Si et al. (2023) |
| Heimbach (UT Austin), Lösch (AWI) | El adjunto de MITgcm y sus aplicaciones polares | Heimbach & Losch (2012); Lösch & Heimbach (2007) |
| Mazloff y Gille (Scripps) | SOSE, adjunto en el océano Austral, calor hacia la plataforma | Mazloff (2011); Rodriguez et al. (2016); Lanham et al. (2025) |
| Morrison y Hogg (ANU) | Modelos globales con remolinos, viento y agua densa | Morrison et al. (2020, 2023) |
| Köhl y Stammer (Hamburgo) | Adjunto en regímenes caóticos y reanálisis | Köhl (1999); Wang et al. (2025) |

## 9. Orden de lectura sugerido

1. **Thompson et al. (2018):** la revisión del frente del talud. Fíjate en la clasificación dinámica y en la sección de intercambio transversal.
2. **Stewart & Thompson (2014):** el mecanismo de remolinos y la frase "sensitive to all aspects of forcing and geometry". Mira cómo definen el transporte a través del quiebre.
3. **Marotzke et al. (1999):** qué es una sensibilidad adjunta de un transporte de calor y cómo se valida contra perturbaciones.
4. **Mazloff (2011):** el adjunto de Mazloff sobre un transporte y respecto al viento. Las escalas de 100 días y de 100 km son las que el Exp. 5 debería comparar.
5. **McLay (2006):** el canal reentrante con remolinos y el adjunto. Lee los capítulos de tiempo no lineal y de prueba de gradiente.
6. **Lösch & Heimbach (2007):** cómo se calcula y se interpreta la sensibilidad a la batimetría.
7. **Rodriguez, Mazloff & Gille (2016):** la hipótesis del rotacional del viento sobre el quiebre.
8. **Gaikwad et al. (2024) 📁:** qué paquetes soporta Tapenade y cómo se hacen las pruebas. Ya lo tienes.

## 10. Términos de búsqueda

- **Física:** Antarctic Slope Front/Current, cross-slope / shoreward / onshore heat transport, Circumpolar Deep Water intrusion, shelf break, eddy stirring, Ekman downwelling, dense shelf water export, reentrant channel.
- **Método:** adjoint sensitivity, algorithmic/automatic differentiation, TAF, TAMC, Tapenade, OpenAD, tangent linear, grdchk, 4D-Var, state estimation, ECCO, SOSE, B-SOSE.
- **Límites:** predictability limit, ensemble adjoint, shadowing, Lyapunov, climate sensitivity of chaotic systems.

```
("adjoint" OR "automatic differentiation" OR "Tapenade") AND ("Antarctic Slope Current" OR "continental shelf break" OR "shoreward heat")
("reentrant channel" OR "zonally periodic channel") AND ("adjoint" OR "sensitivity") AND ("eddy" OR "eddies")
("Circumpolar Deep Water" OR "CDW") AND ("wind stress curl" OR "easterlies") AND ("sensitivity" OR "perturbation experiment")
```

## Registro de búsquedas (Consensus, 2026-11-14)

| # | Búsqueda | Filtros | Artículos | Estado |
|---|---|---|---|---|
| 0 | adjoint sensitivity of onshore heat transport across the Antarctic continental slope | — | 18 | bien (ninguno con adjunto) |
| 1 | idealized reentrant channel model of wind-driven Antarctic Slope Current and cross-slope CDW transport | — | 12 | bien |
| 2 | adjoint sensitivity analysis with the MITgcm adjoint model ocean circulation | — | 12 | bien |
| 3 | adjoint sensitivity of ice shelf melt or Antarctic continental shelf heat content to wind and ocean state | — | 11 | bien |
| 4 | adjoint sensitivity in chaotic eddy-resolving ocean models: exponential gradient growth and assimilation window limits | — | 12 | bien |
| 5 | Mazloff Southern Ocean State Estimate adjoint sensitivity heat transport | — | 12 | bien |
| 6 | review of cross-shelf exchange and heat delivery to the Antarctic continental shelf | — | 10 | bien |
| 7 | review of adjoint methods and algorithmic differentiation for ocean sensitivity analysis and state estimation | — | 10 | bien |
| 8 | adjoint sensitivity Southern Ocean transport to wind stress and topography | year_max 2015 | 10 | bien |
| 9 | adjoint sensitivity Southern Ocean Antarctic shelf heat transport wind stress | year_min 2021 | 10 | bien (ninguno con adjunto sobre la plataforma) |
| 10 | sensitivity of eddy-driven shoreward heat flux across the Antarctic shelf break to easterly wind strength and continental slope steepness | year_min 2015 | 10 | bien |

- **Búsquedas:** 11 hechas, 11 con resultados y ninguna fallida.
- **Artículos:** 127 devueltos en total y unos 95 distintos. La bibliografía tiene 59 entradas: 58 de Consensus y una (Tber et al., 2007) de la biblioteca de la tesis.
- **Tomado de citas y no devuelto por Consensus:** Stewart & Thompson (2015), citado por Bull et al. (2025); Verdy et al. (2023), citado por Köhl et al. (2026). Hay que verificarlos antes de citarlos.
- **Cambio de vocabulario:** antes de 2015 se habla de "ACC transport sensitivity" y "eddy saturation" en canales. Desde 2018 domina "Antarctic Slope Current/Front" y "shoreward heat transport". En los dos periodos la palabra "adjoint" casi no aparece junto a la plataforma.

## Bibliografía

Cada entrada lleva su enlace de Consensus. Los DOI están en `readinglist.md`.

- Abernathey, Marshall & Ferreira (2011). The dependence of Southern Ocean meridional overturning on wind stress. *JPO*. [Consensus](https://consensus.app/papers/details/w:4bf12eb9-c5cb-48b8-8899-6264d73545a6/)
- Bull et al. (2021). Remote control of Filchner-Ronne Ice Shelf melt rates by the Antarctic Slope Current. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:a837e130-78ae-4a6a-9e10-281a2fd54250/)
- Bull, Munday & Jenkins (2025). Influence of topography and winds on the distribution of water masses on the Antarctic continental shelf. *JPO*. [Consensus](https://consensus.app/papers/details/w:ad9f93a3-c781-4e40-a99d-acccdcde431e/)
- Daae et al. (2017). On the effect of topography and wind on warm water inflow: an idealized study of the southern Weddell Sea. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:bc68af9c-9cb3-4bea-a741-b3de5c277f42/)
- Dettling et al. (2024). Toward parameterizing eddy-mediated transport of warm deep water across the Weddell Sea continental slope. *JPO*. [Consensus](https://consensus.app/papers/details/w:3869c682-5cf1-4450-8261-fa83b79f41c6/)
- Dettling, Nakayama et al. (2025). A topography-aware eddy parameterization improves warm water transport across the Cape Darnley continental slope. *JAMES*. [Consensus](https://consensus.app/papers/details/w:264f1a1f-4aea-4ed6-be1f-fa120ed6b191/)
- Dinniman, Klinck & Hofmann (2012). Sensitivity of CDW transport and ice shelf basal melt along the West Antarctic Peninsula to changes in the winds. *J. Climate*. [Consensus](https://consensus.app/papers/details/w:01e385b0-1745-40f4-8f54-914b7c92efc7/)
- Foppert, Rintoul & England (2019). Along-slope variability of cross-slope eddy transport in East Antarctica. *GRL*. [Consensus](https://consensus.app/papers/details/w:1210eae2-d049-4123-84ff-93adf84220a3/)
- Forget et al. (2015). ECCO version 4. *GMD*. [Consensus](https://consensus.app/papers/details/w:e68c0f4c-0b98-4184-8943-45c75e7ffde6/)
- Gaikwad et al. (2024). MITgcm-AD v2. *Future Gener. Comput. Syst.* 📁 [Consensus](https://consensus.app/papers/details/w:40c74e63-e420-48e2-ba9b-08427c874b45/)
- Gebbie & Hsieh (2017). Controllability, not chaos, key criterion for ocean state estimation. *NPG*. [Consensus](https://consensus.app/papers/details/w:bef055f1-5f78-4c67-9b48-a8169e2e8f10/)
- Giering & Kaminski (1998). Recipes for adjoint code construction. *ACM TOMS*. [Consensus](https://consensus.app/papers/details/w:6cb50592-74f0-45fa-b131-e19db1075552/)
- Gopalakrishnan, Cornuelle & Hoteit (2013). Adjoint sensitivity studies of loop current and eddy shedding in the Gulf of Mexico. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:d65cc39a-8df7-494a-93ab-68cf6f707fe1/)
- Haigh et al. (2026). Wind-driven coastal polynya variability drives decadal ice-shelf melt variability in the Amundsen Sea. *GRL*. [Consensus](https://consensus.app/papers/details/w:f918e275-771b-4649-ba06-76dab7ff2d32/)
- Heimbach, Hill & Giering (2005). An efficient exact adjoint of the parallel MIT General Circulation Model. *Future Gener. Comput. Syst.* 📁 [Consensus](https://consensus.app/papers/details/w:23301272-35f0-4958-a7c1-3706c8faf6c3/)
- Heimbach et al. (2010). On the formulation of sea-ice models. Part 2: adjoint sea-ice export sensitivities through the CAA. *Ocean Modelling*. [Consensus](https://consensus.app/papers/details/w:020572f2-1b4a-46c4-996d-3c95cf8408ec/)
- Heimbach & Losch (2012). Adjoint sensitivities of sub-ice-shelf melt rates to ocean circulation under the Pine Island Ice Shelf. *Ann. Glaciol.* [Consensus](https://consensus.app/papers/details/w:6341e7ce-98ef-4fdb-9799-425e3959f0f8/)
- Hogg & Munday (2014). Does the sensitivity of Southern Ocean circulation depend upon bathymetric details? *Phil. Trans. R. Soc. A*. [Consensus](https://consensus.app/papers/details/w:5cb69614-4a4e-4a80-96c6-aa19f00c00c2/)
- Jones et al. (2019). Heat distribution in the Southeast Pacific is only weakly sensitive to high-latitude heat flux and wind stress. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:f0b5fdf1-438a-44e9-91cc-0152d16d2317/)
- Jones et al. (2020). The sensitivity of Southeast Pacific heat distribution to local and remote changes in ocean properties. *JPO*. [Consensus](https://consensus.app/papers/details/w:1ec698af-3422-429d-8c3d-6839e708fe1e/)
- Köhl (1999). An adjoint method for the assimilation of statistical characteristics into eddy-resolving ocean models. Tesis, Kiel. [Consensus](https://consensus.app/papers/details/w:be3a19bb-5013-4c43-acab-79395967b06f/)
- Köhl et al. (2026). A pilot variational coupled reanalysis based on the CESAM climate model. *QJRMS*. [Consensus](https://consensus.app/papers/details/w:ab70834f-93e6-4852-80a2-e589eec01662/)
- Lanham, Mazloff et al. (2025). Seasonal regimes of warm CDW intrusion toward Antarctic ice shelves. *Commun. Earth Environ.* [Consensus](https://consensus.app/papers/details/w:0639a17a-2b8c-4ede-a10c-3e5eb96b0932/)
- Lea, Allen & Haine (2000). Sensitivity analysis of the climate of a chaotic system. *Tellus A*. [Consensus](https://consensus.app/papers/details/w:21df84d7-34b1-4d05-83ab-3503dd06f4e6/)
- Liu et al. (2022). Topography-mediated transport of warm deep water across the continental shelf slope, East Antarctica. *JPO*. [Consensus](https://consensus.app/papers/details/w:26cf84da-374b-4987-99ea-3495c1ac351b/)
- Lösch & Heimbach (2007). Adjoint sensitivity of an ocean general circulation model to bottom topography. *JPO*. [Consensus](https://consensus.app/papers/details/w:3badd2e5-1fc8-43f9-bced-a513b2c1236e/)
- Marotzke et al. (1999). Construction of the adjoint MIT ocean GCM and application to Atlantic heat transport sensitivity. *JGR*. [Consensus](https://consensus.app/papers/details/w:969b0720-5dfa-4965-b9ad-ff5fc4810827/)
- Matsuta, Tanaka & Kubokawa (2026). Drag-controlled regime transitions in the eddy saturation mechanism of the ACC. *GRL*. [Consensus](https://consensus.app/papers/details/w:b9337cd8-afd6-4a7e-adce-7c352403d821/)
- Mazloff, Heimbach & Wunsch (2010). An eddy-permitting Southern Ocean State Estimate. *JPO*. [Consensus](https://consensus.app/papers/details/w:baf8be8c-6fe1-4b43-a1bc-6461ca2852be/)
- Mazloff (2011). On the sensitivity of the Drake Passage transport to air–sea momentum flux. *J. Climate*. [Consensus](https://consensus.app/papers/details/w:09a3b2a3-d980-46fc-b0ae-d73e8dadc683/)
- McLay (2006). The use of adjoint models for determining the sensitivity of integral quantities in an eddy resolving ocean GCM. Tesis, Southampton. [Consensus](https://consensus.app/papers/details/w:524ad0a0-52d2-437d-845b-dc342b6a530c/)
- Morrison & Hogg (2012). On the relationship between Southern Ocean overturning and ACC transport. *JPO*. [Consensus](https://consensus.app/papers/details/w:0beaf48a-43ca-478c-88d7-aa241bf81309/)
- Morrison et al. (2020). Warm CDW transport toward Antarctica driven by local dense water export in canyons. *Sci. Adv.* [Consensus](https://consensus.app/papers/details/w:7893a533-b26b-4cd6-b80f-76961bfa3b04/)
- Morrison et al. (2023). Sensitivity of Antarctic shelf waters and abyssal overturning to local winds. *J. Climate*. [Consensus](https://consensus.app/papers/details/w:e3b1a1f7-d87e-4dd2-9c09-30652fb28c41/)
- Munday & Zhai (2015). Sensitivity of Southern Ocean circulation to wind stress changes: role of relative wind stress. *Ocean Modelling*. [Consensus](https://consensus.app/papers/details/w:985c16ec-c126-4548-aa6a-cfdc07d950a2/)
- Na (2025). State estimation and adjoint sensitivity analysis in the South China Sea. Tesis, HKUST. [Consensus](https://consensus.app/papers/details/w:7afbc492-501a-4aad-a444-d33a15c0f838/)
- Nakayama et al. (2021). Development of adjoint-based ocean state estimation for the Amundsen and Bellingshausen seas and ice shelf cavities using MITgcm–ECCO (66j). *GMD*. [Consensus](https://consensus.app/papers/details/w:d3d40f67-7df5-4e9a-bba5-75b606e1e65e/)
- Ong et al. (2024). Intrinsically episodic Antarctic shelf intrusions of CDW via canyons. *JPO*. [Consensus](https://consensus.app/papers/details/w:8b52513d-3586-4af7-a8f7-c99655a9e55f/)
- Palóczy, Gille & McClean (2018). Oceanic heat delivery to the Antarctic continental shelf: large-scale, low-frequency variability. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:69120f13-e28e-403f-b0d7-cfab2a77c262/)
- Pimm et al. (2024). Surface heat fluxes drive a two-phase response in Southern Ocean mode water stratification. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:d86a5612-62a8-4647-8e64-d6e85deb1c55/)
- Pimm et al. (2025). Local versus far-field control on South Pacific Subantarctic mode water variability. *Ocean Sci.* [Consensus](https://consensus.app/papers/details/w:ae74ff2b-1437-457a-b4d3-a1376733aa66/)
- Rieke et al. (2025). Interannual variability of ocean heat content in the eastern Amundsen Sea controlled by wind. *GRL*. [Consensus](https://consensus.app/papers/details/w:a2cabee4-df97-4119-81ea-f11add31b318/)
- Rodriguez, Mazloff & Gille (2016). An oceanic heat transport pathway to the Amundsen Sea Embayment. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:eeabb40f-d4c1-4208-8a4b-c413e262f11e/)
- Shimizu (2024). Fundamental properties of adjoint model and adjoint sensitivity under fully nonlinear hydrostatic internal gravity waves. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:c82dcd3d-a14f-4837-9d92-184ce0edb9a0/)
- Si, Stewart & Eisenman (2023). Heat transport across the Antarctic Slope Front controlled by cross-slope salinity gradients. *Sci. Adv.* [Consensus](https://consensus.app/papers/details/w:f8de910b-e66e-45b3-ba7b-10da13f73572/)
- Si et al. (2024). Antarctic Slope Undercurrent and onshore heat transport driven by ice shelf melting. *Sci. Adv.* [Consensus](https://consensus.app/papers/details/w:a4f4b9c4-5a9c-4850-9146-fecce8b1cb90/)
- Stewart & Thompson (2012). Sensitivity of the ocean's deep overturning circulation to easterly Antarctic winds. *GRL*. [Consensus](https://consensus.app/papers/details/w:f5ef094e-ca75-4f43-92ae-a3dfc684cb27/)
- Stewart & Thompson (2014). Eddy-mediated transport of warm CDW across the Antarctic shelf break. *GRL*. [Consensus](https://consensus.app/papers/details/w:3b91cb90-f980-49d5-82ff-d393e3f31b81/)
- Stewart & Thompson (2016). Eddy generation and jet formation via dense water outflows across the Antarctic continental slope. *JPO*. [Consensus](https://consensus.app/papers/details/w:b4567220-8935-4c52-97e3-3df8a49d7401/)
- Stewart, Klocker & Menemenlis (2018). Circum-Antarctic shoreward heat transport derived from an eddy- and tide-resolving simulation. *GRL*. [Consensus](https://consensus.app/papers/details/w:b1b3fd4d-9476-4b8b-87c2-00320ae5c0f8/)
- Stewart, Klocker & Menemenlis (2019). Acceleration and overturning of the Antarctic Slope Current by winds, eddies, and tides. *JPO*. [Consensus](https://consensus.app/papers/details/w:99a38c06-60e7-4477-a891-1f41e68aa707/)
- Tber, Hascoët, Vidard & Dauvergne (2007). Building the tangent and adjoint codes of OPA with TAPENADE. 📁 (no lo devolvió Consensus)
- Thakur & Nadarajah (2025). A stabilized march approach to adjoint-based sensitivity analysis of chaotic flows. arXiv (preprint, sin revisión por pares). [Consensus](https://consensus.app/papers/details/w:7ebab769-b1e6-4857-b5f5-a8cf3e98bacf/)
- Thompson et al. (2018). The Antarctic Slope Current in a changing climate. *Rev. Geophys.* [Consensus](https://consensus.app/papers/details/w:f9e09b78-1830-485b-9856-6a82f3025ab8/)
- Vlasenko, Köhl & Stammer (2015). The efficiency of geophysical adjoint codes generated by automatic differentiation tools. *Comput. Phys. Commun.* [Consensus](https://consensus.app/papers/details/w:7c0875d8-479b-443c-9729-b910cc49046e/)
- Wang, Köhl & Stammer (2025). Northeast Pacific marine heatwave mechanism inferred from adjoint sensitivities. *J. Climate*. [Consensus](https://consensus.app/papers/details/w:305435cb-206d-4479-ac2b-715eca4525da/)
- Wang et al. (2024). A model study of buoyancy driven cross-isobath transport over the Ross Sea continental shelf break. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:1d99d4a0-6d65-449b-8adf-38347c2031a9/)
- Wang & Zaki (2026). Mitigating adjoint chaos in wall turbulence. *JFM*. [Consensus](https://consensus.app/papers/details/w:487a27c0-f076-4ad5-a90c-634e9b8c242c/)
- Zhan et al. (2018). Sensitivity studies of the Red Sea eddies using adjoint method. *JGR Oceans*. [Consensus](https://consensus.app/papers/details/w:4702a708-62f7-485b-8be9-17bc631bb863/)
