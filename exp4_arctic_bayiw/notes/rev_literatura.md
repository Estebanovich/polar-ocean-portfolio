# Revisión de literatura del Exp. 4: BayIW ártico

Fecha: 2026-11-11. Se hizo antes de diseñar el experimento y antes de correr nada.

**Pregunta del Exp. 4** (spec del portafolio, §Exp. 4): cuando *f* es grande (~75–80°N), hay una haloclina ártica y la marea es débil, ¿cambia la partición entre la energía cercana a la inercial y la energía de seiche que la bahía irradia? ¿Qué pasa si una capa de hielo frena la superficie?

**Fuente.**
- **Consensus:** diez búsquedas, más la búsqueda inicial de reconocimiento (detalle en el registro de búsquedas, al final).
- **Biblioteca del portafolio** (`literatura/`): no tiene nada del tema.
- **Biblioteca de la tesis** (`~/Documents/TESIS_InternalWaves-2/Literatura/`, 268 artículos en `.md`): se revisó con `litsearch.py` después de las búsquedas, cuando Esteban la mencionó. Tiene 12 de los artículos citados aquí (marcados 📁). También tiene varios pertinentes que Consensus no devolvió (lista al final de §2). No tiene nada específico del Ártico ni del hielo marino.

Los marcados ⬜ no están en ninguna de las dos bibliotecas. Lo que dice este documento de ellos sale **solo de su resumen**, y hay que leerlos completos antes de citarlos en el README. Todos están en `readinglist.md`, con su DOI.

**Un dato de la tesis que condiciona todo.**
- **Tesis:** f₀ = 6.97×10⁻⁵ s⁻¹ (unos 28.6°N), con un periodo inercial de unas 25 h. La seiche barotrópica principal dura 3.4–4.3 h, y el modo baroclínico de seiche 4.5–6.5 h (postproceso de `Output_IWBay`).
- **A 75°N y 80°N:** el periodo inercial es de 12.4 h y 12.2 h.
- **Consecuencia:** la seiche sigue por encima de f y puede seguir irradiando como onda interna libre, pero las dos bandas se acercan, de un factor ~6 a uno ~3 en periodo.
- **Corrección:** en el checkpoint de esta revisión dije que la seiche podría quedar por debajo de f a esa latitud. No es así para esta bahía.
- **La latitud crítica de la marea M2 (74.5°N)** importa para la marea, no para la seiche.

## 1. Lo esencial en cinco puntos

1. **El hielo explica solo una parte de por qué el Ártico tiene poca energía de ondas internas.**
   - La explicación clásica es que el hielo impide que el viento meta energía (Rainville & Woodgate, 2009) y que amortigua las ondas por fricción en su cara inferior (Rippeth & Fine, 2022).
   - Con menos hielo, la energía cercana a la inercial aumenta (Dosser & Rainville, 2016).
   - Pero "Not just sea ice": también importan otros factores (Guthrie & Morison, 2021). Cole et al. (2018) separan la concentración del hielo, su rugosidad y la estratificación de la capa superficial.
   - La mezcla en el interior sigue siendo baja, aunque cambia según la región (Guthrie et al., 2013; Fer, 2014; Lincoln et al., 2016; Dosser et al., 2021; Fine & Cole, 2022).
2. **El hielo que se mueve no siempre frena.**
   - La oscilación inercial del propio hielo puede *generar* ondas cercanas a la inercial, sobre todo cuando el hielo se forma y entra en resonancia con el agua (Kawaguchi et al., 2019).
   - El hielo en deriva libre también genera ondas internas (De Abreu & Timmermans, 2026; McPhee, 1978).
   - El hielo fijo a la costa (*landfast*) sí controla y reduce la turbulencia superficial (Inall et al., 2021, en la Antártida).
   - Para el Exp. 4 esto significa que "con hielo" no es un solo caso: hay que separar el hielo **fijo** del hielo **libre**.
3. **La generación costera de ondas cercanas a la inercial tiene teoría y observaciones recientes, pero solo en latitudes medias.**
   - La pared de la costa obliga a que la oscilación inercial produzca velocidad vertical, y esa velocidad bombea ondas internas y superficiales que salen mar adentro (Kelly, 2019, con un estimado de ~10 % del trabajo del viento global en ondas costeras barotrópicas).
   - En el lago Superior, las velocidades cercanas a la inercial se inhiben dentro de dos radios de Rossby internos de la costa, y la onda que sale lleva unos 10 W por metro de costa (Kelly et al., 2024). Hay observaciones en una cuenca cerrada (Green et al., 2025).
   - Los clásicos son Kundu, Chao & McCreary (1983) y Kundu (1984). Los modelos con frentes son de Davies & Xing (2002, 2005) y Xing & Davies (2004).
   - El trabajo de tu tesis cae justo aquí, y el grupo de Kelly (Duluth, con Scripps en el artículo de 2024) es el más cercano.
4. **En el Ártico, la generación costera casi no se ha estudiado.**
   - Lo más cercano es Carmack & Kulikov (1998): surgencia forzada por viento y generación de **ondas de Kelvin internas** en el cañón Mackenzie, mar de Beaufort, en aguas libres de hielo.
   - Hughes & Klymak (2019) muestran que, en canales del Archipiélago Ártico Canadiense, al norte de la latitud crítica, las paredes permiten que la marea interna irradie como ondas de Kelvin.
   - Ninguna de las búsquedas devolvió un estudio de una bahía o costa ártica que irradie ondas cercanas a la inercial con y sin hielo. Ese es el hueco que puede llenar el Exp. 4.
5. **La teoría de cuencas en rotación da el marco para interpretar el cambio de f.**
   - En lagos grandes, el número de Burger (radio de Rossby interno entre el tamaño de la cuenca) decide si dominan las ondas de Kelvin internas o las de Poincaré, y cómo se reparte la energía (Antenucci & Imberger, 2001; Antenucci et al., 2000; Shimizu & Imberger, 2009; Ulloa et al., 2015).
   - A 75–80°N, f se duplica respecto a la tesis, así que el radio de Rossby se reduce más o menos a la mitad para la misma estratificación. La haloclina, además, cambia N(z).
   - La tesis ya identificó una onda de Kelvin interna (§3.2.4). En el Ártico debería estar más atrapada a la costa.

## 2. Orden de lectura sugerido

1. **Alford, MacKinnon et al. (2016)**, *Near-inertial internal gravity waves in the ocean* (*Annu. Rev. Mar. Sci.*, 📁). Es la mejor visión de conjunto y es de la investigadora objetivo. Leer la parte de generación por viento y la de modos verticales.
2. **Rippeth & Fine (2022)**, *Turbulent mixing in a changing Arctic Ocean* (*Oceanography*, ⬜, acceso abierto). Lo esencial del Ártico: por qué hay poca energía y cómo cambia con el hielo. Ver la diferencia entre plataformas y cuencas.
3. **Kelly (2019)**, *Coastally generated near-inertial waves* (*JPO*, 📁). Es la teoría más cercana a la bahía. Ver la estimación del flujo de energía en una costa vertical y cómo se extiende a una topografía cualquiera.
4. **Kelly et al. (2024)**, *Direct observations of coastally generated near-inertial waves* (*JGR*, 📁). La separación entre "flujo forzado" y "flujo de onda" es una forma de diagnóstico que el Exp. 4 puede copiar.
5. **Dosser & Rainville (2016)** y **Martini et al. (2014)** (*JPO*, ⬜). Cómo cambia el campo de ondas cercanas a la inercial entre verano sin hielo e invierno con hielo, en el mar de Beaufort.
6. **Guthrie & Morison (2021)** (*GRL*, ⬜). El contrapunto: no todo es hielo.
7. **Antenucci & Imberger (2001)** (*L&O*, ⬜; Antenucci et al. 2000 sí está 📁). Dispersión y energía de ondas internas de cuenca con rotación. Sirve para predecir los modos del caso ártico antes de correr.
8. **Carmack & Kulikov (1998)** y **Hughes & Klymak (2019)** (⬜). Los dos únicos precedentes árticos costeros que salieron.

**Ya en la biblioteca de la tesis y pertinentes, aunque Consensus no los devolvió** (📁):
- **Davidson (2001):** la asimetría en la respuesta de una bahía estratificada al viento. Es el más cercano a la configuración de la bahía.
- **Alford (2020a):** el trabajo del viento cercano a la inercial con modelos de losa y esfuerzo relativo. El esfuerzo relativo es la forma natural de meter el hielo libre.
- **Pollard (1970) y Pollard & Millard (1970):** el modelo de losa.
- **Gill (1984):** ondas internas en la estela de tormentas, con la proyección modal.
- **Zhai et al. (2009):** la pérdida de energía cercana a la inercial en mezcla.
- **Ondas internas de cuenca en lagos:** Beletsky et al. (1997), Hodges et al. (2000) y Lemmin et al. (2005), con ondas de Kelvin y Poincaré.
- **Seiches:** Arneborg & Liljebladh (2001), en el fiordo Gullmar, y Cushman-Roisin et al. (2005), sobre seiches con estratificación.

## 3. Cómo llegó el campo aquí

| Año | Hito | Por qué importa |
|---|---|---|
| 1978 | McPhee: oscilación inercial en hielo a la deriva | El hielo libre oscila a f y transmite esa oscilación al agua |
| 1983–1984 | Kundu et al.; Kundu: corrientes costeras transitorias y oscilaciones inerciales costeras | Teoría clásica de la respuesta costera al viento con rotación |
| 1998 | Carmack & Kulikov: ondas de Kelvin internas en el cañón Mackenzie | Primer caso ártico costero forzado por viento |
| 2000–2001 | Antenucci & Imberger: ondas internas de cuenca en lagos con rotación | Marco Kelvin/Poincaré y número de Burger |
| 2005 | Pinkel: propagación de ondas cercanas a la inercial en el Ártico occidental; MacKinnon & Gregg: ondas cercanas a la inercial en la plataforma de Nueva Inglaterra | El Ártico como océano de baja energía |
| 2009 | Rainville & Woodgate: generación de ondas internas en el Ártico estacionalmente sin hielo | Sin hielo, el viento sí mete energía |
| 2014–2016 | Martini et al.; Fer; Dosser & Rainville; Lincoln et al. | Primeras series largas; la respuesta depende de la estación y del hielo |
| 2016 | Alford, MacKinnon et al., revisión en *Annual Review* | Síntesis de la física de las ondas cercanas a la inercial |
| 2019 | Kelly: ondas cercanas a la inercial generadas en la costa (teoría) | La costa como fuente de ondas cercanas a la inercial |
| 2021–2022 | Guthrie & Morison ("not just sea ice"); Dosser et al.; Fine & Cole; Rippeth & Fine (revisión) | Se matiza el papel del hielo |
| 2024–2026 | Kelly et al. y Green et al. (lago Superior); Muilwijk et al. (hielo más débil, más esfuerzo); De Abreu & Timmermans (deriva libre) | Observación directa de generación costera; el hielo móvil como transmisor |

**Términos.**
- **Ondas:** antes de 2010 se habla de *inertial oscillations* y *inertio-gravity waves*. Después, de *near-inertial waves (NIWs)*, *wind work* y *energy flux*.
- **Esfuerzo bajo el hielo:** se busca como *ice-ocean stress*, *ocean surface stress* y *momentum transfer*.
- **Bahías y lagos:** la literatura de cuencas usa *basin-scale internal waves*, *internal seiche*, *Kelvin/Poincaré waves* y *Burger number*.

## 4. Sub-áreas

### 4.1 Ondas cercanas a la inercial en el Ártico y el papel del hielo

**Qué muestra.**
- **Pocas ondas por el hielo, pero no solo por él:** la energía de ondas internas del Ártico es un orden de magnitud menor que en latitudes medias, y el hielo es la causa principal pero no la única (Rainville & Woodgate, 2009; Guthrie & Morison, 2021).
- **Con hielo libre:** si el hielo se mueve libremente, su oscilación inercial genera ondas cercanas a la inercial (Kawaguchi et al., 2019; De Abreu & Timmermans, 2026).
- **Mezcla:** sigue siendo baja en las cuencas, pero aumenta a profundidades intermedias cuando no hay hielo (Lincoln et al., 2016). En registros largos se observan cambios regionales (Dosser et al., 2021; Fine & Cole, 2022).

**Artículos clave:**
- Rainville & Woodgate (2009). *GRL*. 161 citas. https://doi.org/10.1029/2009GL041291
- Dosser & Rainville (2016). *JPO*. 49 citas. https://doi.org/10.1175/JPO-D-15-0056.1
- Martini et al. (2014). *JPO*. 40 citas; del verano sin hielo al invierno con hielo. https://doi.org/10.1175/JPO-D-13-0160.1
- Guthrie & Morison (2021). *GRL*. https://doi.org/10.1029/2020GL090508
- Cole et al. (2018). *JGR Oceans*; concentración, rugosidad y estratificación. https://doi.org/10.1029/2018JC014096
- Kawaguchi et al. (2019). *Polar Science*; resonancia hielo–agua. https://doi.org/10.1016/j.polar.2019.01.004
- Fer (2014). *JPO*; mezcla cercana a la inercial en el Ártico central. https://doi.org/10.1175/JPO-D-13-0133.1
- Lincoln et al. (2016). *GRL*; mezcla por viento en un Ártico sin hielo. https://doi.org/10.1002/2016GL070454

**Términos de búsqueda:** near-inertial waves, NIW, Arctic Ocean, sea ice, ice-tethered profiler, Beaufort Sea, Canada Basin, wind work, internal wave energy, inertial oscillation.

```
("near-inertial" OR "inertial oscillation*") AND ("Arctic" OR "Beaufort" OR "Canada Basin") AND ("sea ice" OR "ice-free")
("internal wave*" AND "Arctic") AND ("mixing" OR "dissipation") AND ("sea ice decline" OR "ice loss")
```

### 4.2 Generación costera de ondas internas por viento

**Qué muestra.**
- **Mecanismo:** la condición de no flujo normal en la costa convierte la oscilación inercial de la capa de mezcla en bombeo vertical, que genera ondas internas y superficiales que salen mar adentro (Kelly, 2019).
- **Observaciones:** en el lago Superior, la oscilación se inhibe dentro de ~2 radios de Rossby internos y una fracción pequeña sale como onda (Kelly et al., 2024). La mayor parte del trabajo del viento se pierde directamente en turbulencia.
- **Modelos:** los modelos 2D y 3D con fronteras y frentes muestran ondas superinerciales que salen de la costa y ondas subinerciales atrapadas (Davies & Xing, 2002, 2005; Xing & Davies, 2004).
- **Seiches con inercia:** Chen, Chen & Xing (2017) simulan, en una cuenca 2D, seiches, oscilaciones inerciales y ondas cercanas a la inercial a la vez. Es un análogo directo de la partición de tu tesis.

**Artículos clave:**
- Kelly (2019). *JPO*. 📁 https://doi.org/10.1175/JPO-D-18-0148.1
- Kelly et al. (2024). *JGR Oceans*. 📁 https://doi.org/10.1029/2024JC020932
- Green et al. (2025). *JGR Oceans*; cuenca cerrada. 📁 https://doi.org/10.1029/2024JC021828
- Chen, Chen & Xing (2017). *Ocean Science*; cuenca 2D con seiches e inercia. 📁 https://doi.org/10.5194/os-13-829-2017
- Kundu, Chao & McCreary (1983). *Deep-Sea Res.* 📁 https://doi.org/10.1016/0198-0149(83)90061-4
- Davies & Xing (2005). *JPO*; frente de fondo. https://doi.org/10.1175/JPO2732.1
- MacKinnon & Gregg (2005). *JPO*; ondas cercanas a la inercial en una plataforma con estratificación que cambia. https://doi.org/10.1175/JPO2822.1
- Carmack & Kulikov (1998). *JGR*; caso ártico. https://doi.org/10.1029/98JC00113

**Términos de búsqueda:** coastally generated near-inertial waves, coastal inhibition, inertial pumping, offshore energy flux, Rossby radius, coastal boundary, wind impulse, enclosed basin, embayment.

```
("near-inertial" OR "inertial") AND ("coast*" OR "bay" OR "embayment" OR "enclosed basin") AND ("generation" OR "radiation")
("internal Kelvin wave*" OR "coastal-trapped wave*") AND ("wind-forced" OR "wind event") AND ("Arctic" OR "high latitude")
```

### 4.3 Cómo transmite el hielo el viento al océano

**Qué muestra.**
- **Dos regímenes de hielo:** el esfuerzo que llega al océano depende de la concentración, del grosor y de la rugosidad del hielo. Un hielo más delgado y débil transmite más esfuerzo (Martin et al., 2014, 2016; Muilwijk et al., 2024). Con concentración < 85 % el hielo está en deriva libre (Heorton et al., 2019).
- **Hielo fijo:** el hielo fijo a la costa reduce la turbulencia de la capa superficial (Inall et al., 2021).
- **Coeficientes:** el coeficiente de arrastre hielo–océano y el ángulo de giro se miden en campañas como MOSAiC (Kawaguchi et al., 2024; Fer et al., 2022) y en mediciones de rugosidad (Cole et al., 2017).
- **Para un modelo idealizado, hay dos extremos simples:**
  - **Hielo fijo:** sin esfuerzo del viento y con arrastre hielo–agua en la superficie.
  - **Hielo libre:** una losa con inercia y arrastre, al estilo de McPhee (1978).

**Artículos clave:**
- McPhee (1978). *Dyn. Atmos. Oceans*; oscilación inercial en hielo a la deriva. https://doi.org/10.1016/0377-0265(78)90005-2
- Martin, Steele & Zhang (2014). *JGR*. 153 citas. https://doi.org/10.1002/2013JC009425
- Martin et al. (2016). *JGR*; rugosidad variable. https://doi.org/10.1002/2015JC011186
- Muilwijk et al. (2024). *Nature Communications*. https://doi.org/10.1038/s41467-024-50874-0
- Cole et al. (2017). *Elementa*; rugosidad y transferencia de momento en la zona marginal. https://doi.org/10.1525/elementa.241
- Heorton et al. (2019). *JGR Oceans*; coeficientes de arrastre. https://doi.org/10.1029/2018JC014881
- Inall et al. (2021). *JGR Oceans*; hielo fijo y turbulencia. https://doi.org/10.1029/2021JC017963
- De Abreu & Timmermans (2026). *JPO*; deriva libre. https://doi.org/10.1175/JPO-D-25-0165.1

**Términos de búsqueda:** ice–ocean stress, ocean surface stress, momentum transfer, free drift, landfast ice, internal ice stress, ice–ocean drag coefficient, turning angle.

```
("landfast ice" OR "fast ice" OR "free drift") AND ("ocean surface stress" OR "momentum transfer") AND ("inertial" OR "internal wave*")
("ice-ocean drag" OR "ice–ocean stress") AND ("model*" OR "parameterization") AND "Arctic"
```

### 4.4 Teoría a latitud alta: cuencas en rotación y latitud crítica

**Qué muestra.**
- **Cuencas en rotación:** las ondas internas de cuenca se separan en ondas de Kelvin, atrapadas a la costa y subinerciales, y ondas de Poincaré, superinerciales. El número de Burger fija su frecuencia y su reparto de energía (Antenucci & Imberger, 2001). La amortiguación y la degeneración no lineal de las de Kelvin están estudiadas en laboratorio y en simulación numérica directa (DNS) (Shimizu & Imberger, 2009; Ulloa et al., 2015).
- **Latitud crítica:** al norte de 74.5°N la marea M2 no puede irradiar como onda interna libre (Vlasenko et al., 2003). Aun así hay conversión y mezcla, por ondas de sotavento inestables (Rippeth et al., 2017) o por ondas de Kelvin en canales (Hughes & Klymak, 2019).
- **Para el Exp. 4:** la marea no entra (no hay marea en la configuración). Lo que importa es que, a 75–80°N, f queda muy cerca de la frecuencia de la seiche baroclínica de modo 1 de la tesis (4.5–6.5 h, contra 12.2–12.4 h de periodo inercial). Hay que verificar con N(z) de la haloclina qué modos quedan por arriba y por abajo de f.

**Artículos clave:**
- Antenucci & Imberger (2001). *Limnol. Oceanogr.*; energía de ondas internas largas en lagos grandes. https://doi.org/10.4319/lo.2001.46.7.1760
- Antenucci, Imberger & Saggio (2000). *Limnol. Oceanogr.*; lago Kinneret. 📁 https://doi.org/10.4319/lo.2000.45.7.1621
- Shimizu & Imberger (2009). *JFM*; amortiguación en cuencas con rotación. https://doi.org/10.1017/S0022112009008039
- Ulloa et al. (2015). *JFM*; degeneración de ondas de Kelvin internas. https://doi.org/10.1017/jfm.2015.311
- Hughes & Klymak (2019). *JPO*; canales al norte de la latitud crítica. https://doi.org/10.1175/JPO-D-18-0132.1
- Rippeth et al. (2017). *GRL*; conversión de marea al norte de la latitud crítica. https://doi.org/10.1002/2017GL075310
- Vlasenko et al. (2003). *Deep-Sea Res. I*; marea interna cerca de la latitud crítica. https://doi.org/10.1016/S0967-0637(03)00018-9
- Zemskova & Lerczak (2025). *JPO*; marea interna en cuencas semicerradas. https://doi.org/10.1175/JPO-D-24-0212.1

**Términos de búsqueda:** Burger number, internal Kelvin wave, Poincaré wave, basin-scale internal wave, internal seiche, rotating basin, critical latitude, subinertial, superinertial.

```
("internal seiche*" OR "basin-scale internal wave*") AND ("rotation" OR "Kelvin" OR "Poincaré") AND ("energy" OR "dispersion")
("critical latitude" OR "subinertial") AND ("internal wave*" OR "internal tide*") AND ("Arctic" OR "Barents" OR "Kara")
```

### 4.5 Modelado de ondas cercanas a la inercial (cobertura débil)

**Qué muestra.**
- **Cobertura:** la búsqueda de modelado bajo hielo se desvió hacia el oleaje superficial en hielo. No salió ningún estudio idealizado de ondas cercanas a la inercial bajo hielo con MITgcm.
- **Lo útil:**
  - Un análisis modal global de las ondas cercanas a la inercial en HYCOM: solo el ~30 % del aporte del viento proyecta en los 5 primeros modos, y el 79 % se disipa localmente (Raja et al., 2022).
  - Un cálculo simple de la fracción local q por proyección modal de una capa de mezcla tipo losa (Alford, 2020).
  - Modelos idealizados de la generación por una tormenta (Voelker et al., 2019) y un modelo de losa con radiación de ondas (Olbers et al., 2020).
  - En el Ártico profundo, Bracamontes-Ramírez, Walter & Losch (2024) estudian la propagación hacia el fondo de la cuenca de Canadá (Losch es desarrollador de MITgcm).

**Artículos clave:**
- Raja et al. (2022). *JPO*. 📁 https://doi.org/10.1175/JPO-D-21-0130.1
- Alford (2020). *JPO*. 📁 https://doi.org/10.1175/JPO-D-20-0106.1
- Voelker et al. (2019). *Dyn. Atmos. Oceans*. https://doi.org/10.1016/j.dynatmoce.2019.03.005
- Olbers, Jurgenowski & Eden (2020). *Ocean Dynamics*. https://doi.org/10.1007/s10236-020-01376-2
- Bracamontes-Ramírez, Walter & Losch (2024). *JGR Oceans*. https://doi.org/10.1029/2023JC020228

**Términos de búsqueda:** modal decomposition, near-inertial energy flux, slab model, wind power input, idealized simulation, MITgcm, seaice package.

```
("near-inertial" AND ("modal decomposition" OR "vertical modes")) AND ("energy flux" OR "wind power")
("MITgcm" OR "idealized model") AND ("near-inertial" OR "inertial oscillation*") AND ("sea ice" OR "ice cover")
```

## 5. Grupos que dominan el tema

- **Alford, MacKinnon y colaboradores (Scripps).** Revisión de las ondas cercanas a la inercial (2016), el Climate Process Team de mezcla (MacKinnon et al., 2017), la fracción local q (Alford, 2020) y la plataforma de Nueva Inglaterra (MacKinnon & Gregg, 2005). Alford apareció en 6 búsquedas y MacKinnon en 7. Es la investigadora objetivo.
- **Kelly (Univ. de Minnesota Duluth).** Generación costera: teoría (2019), lago Superior (2024, 2025) y generación de modos bajos (Stokes, Kelly & Lucas, 2025). Apareció en 4 búsquedas.
- **Rippeth, Fine, Lincoln (Bangor y WHOI).** Mezcla en un Ártico cambiante y latitud crítica. Rippeth apareció en 5 búsquedas.
- **Kawaguchi (Univ. de Tokio).** Hielo, ondas cercanas a la inercial y MOSAiC. Apareció en 5 búsquedas.
- **Cole, Toole, Timmermans (WHOI y Yale).** Hielo, rugosidad, perfiladores atados al hielo (ITP). Cole apareció en 3 búsquedas y Timmermans en 5.

## 6. Huecos y lo que significan para el Exp. 4

- **La costa ártica como fuente de ondas cercanas a la inercial.** La teoría de Kelly (2019) y las observaciones del lago Superior están en latitudes medias, y en el Ártico solo salió Carmack & Kulikov (1998), de ondas de Kelvin internas. Si nadie ha cuantificado cuánta energía irradia una costa ártica con y sin hielo, no sabemos si la pérdida de hielo cerca de la costa aumenta la energía de las ondas internas igual que en mar abierto. El Exp. 4 es una prueba idealizada de eso.
- **El hielo fijo frente al hielo libre.** Las observaciones mezclan los dos regímenes. Si el hielo libre transmite oscilación inercial mientras el fijo la bloquea, la respuesta costera depende del tipo de hielo y no solo de su presencia. Un modelo idealizado puede separar los dos casos de forma limpia. Las observaciones no pueden.
- **La partición seiche / cercana a la inercial cuando las dos bandas se acercan.** La teoría de cuencas en rotación viene de lagos sin bahía abierta al mar. No sabemos cómo cambia, con f grande y una haloclina somera, el reparto que encontró la tesis (la respuesta cercana a la inercial lleva la mayor parte de la energía baroclínica). Esa es la pregunta central del Exp. 4.
- **Modelado bajo hielo.** No salió ningún estudio idealizado con MITgcm de ondas cercanas a la inercial bajo hielo. Representar el hielo como un esfuerzo superficial modificado es una simplificación que hay que declarar. Las bases para elegir la forma están en la sub-área 4.3.

**Implicaciones para el diseño (propuestas para decidir en el brainstorming, no decididas):**
1. **Tres casos de superficie:** agua abierta, hielo fijo (sin esfuerzo del viento y con arrastre hielo–agua) y hielo libre (losa con inercia y arrastre, al estilo de McPhee). El spec solo tenía "con capa de hielo".
2. **Teoría antes de correr:** calcular los modos verticales con N(z) de la haloclina y f a 75°N para predecir qué bandas quedan por arriba de f, cuál es el radio de Rossby interno y el número de Burger de la bahía. Eso da criterios preregistrados.
3. **Diagnóstico tipo Kelly et al. (2024):** separar el flujo forzado del flujo de onda y medir el flujo de energía por metro de costa, comparable con los 10 W m⁻¹ del lago Superior.
4. **Fracción local q por proyección modal (Alford, 2020)**, comparable entre la tesis y el caso ártico.

## 7. Registro de búsquedas

El plan de Consensus permite hasta 300 artículos por búsqueda. Se pidieron 20–30 por búsqueda.

| # | Búsqueda | Filtros | Artículos | Estado |
|---|---|---|---|---|
| 0 | wind-driven near-inertial internal waves Arctic Ocean sea ice halocline energy flux | — | 30 | OK (reconocimiento) |
| 1 | near-inertial internal wave energy Arctic Ocean changes with declining sea ice cover observations moorings | — | 30 | OK |
| 2 | coastal generation of near-inertial internal waves by wind near a coastline or bay, offshore radiation, coastal inhibition of inertial oscillations | — | 30 | OK |
| 3 | momentum transfer from wind to ocean through sea ice, ice internal stress, landfast ice versus free drift, effect on inertial oscillations and ocean surface stress | — | 30 | OK |
| 4 | critical latitude subinertial internal waves high latitudes, M2 internal tide trapping north of 74.5 degrees, superinertial versus subinertial basin seiche radiation | — | 30 | OK |
| 5 | numerical model simulation of wind-generated near-inertial waves under sea ice in the Arctic Ocean, ice-ocean model MITgcm idealized | — | 30 | OK (fuera de tema: casi todo oleaje en hielo) |
| 6 | review of internal waves and turbulent mixing in the Arctic Ocean with sea ice decline | revisiones | 20 | OK |
| 7 | review basin-scale internal waves in rotating basins: internal Kelvin and Poincaré waves, effect of Earth's rotation on internal seiches, Burger number | — | 20 | OK |
| 8 | wind-forced coastal response, inertial oscillations near a coast, generation of internal waves at a coastal boundary, stratified shelf | year_max 2010 | 20 | OK |
| 9 | wind-forced coastal generation of near-inertial and internal Kelvin waves on Arctic shelves and coasts, landfast ice, ice-free season | year_min 2021 | 20 | OK (casi nada ártico costero: es el hueco) |
| 10 | near-inertial internal gravity waves wind generation, energy flux into the ocean interior, mixed layer, vertical modes, dissipation (seguimiento de Alford et al., 2016) | year_min 2017 | 20 | OK |

- **Búsquedas:** 11 ejecutadas, la inicial más las 10 del plan estándar:
  - 5 por sub-área;
  - 2 de revisiones;
  - 2 por época, sobre la generación costera;
  - 1 de seguimiento del más citado.

  Las 11 dieron resultados y ninguna falló.
- **Artículos únicos recibidos:** 218 (deduplicados por título).
- **Citados aquí:** 50 artículos de Consensus: 47 nuevos en `readinglist.md` y 3 que ya estaban (Rainville & Woodgate 2009, Dosser & Rainville 2016, Kawaguchi et al. 2019). Además, 11 de la biblioteca de la tesis que Consensus no devolvió (§2).
- **Cobertura débil:**
  - La sub-área 4.5 (modelado bajo hielo) devolvió casi solo oleaje en hielo.
  - La búsqueda 9 no encontró generación costera ártica. Esto puede ser un hueco real o un problema de terminología; conviene una búsqueda a mano con "Kelvin wave" + "Laptev" / "Kara" / "Chukchi".
- **Bibliotecas locales:**
  - En `literatura/` del portafolio, `litsearch.py` no encontró nada para "near-inertial" ni para "internal waves Arctic".
  - En la biblioteca de la tesis se buscó: near-inertial, inertial oscillation, sea ice, Arctic, Kelvin wave, Poincare, Burger number, coastal, seiche, critical latitude, landfast y halocline. Resultado: 12 de los citados ya están ahí, más los pertinentes listados en §2. Nada de hielo ni del Ártico. "landfast" no dio ningún resultado.

## Para descargar (prioridad)

Rippeth & Fine 2022 · Dosser & Rainville 2016 · Martini et al. 2014 · Guthrie & Morison 2021 · Antenucci & Imberger 2001 · Carmack & Kulikov 1998 · Hughes & Klymak 2019 · McPhee 1978 · Rainville & Woodgate 2009 · Cole et al. 2018. (Alford 2016, Kelly 2019, Kelly 2024 y Chen 2017 ya están en la biblioteca de la tesis.) Los DOI están en `readinglist.md`, sección Exp. 4.
