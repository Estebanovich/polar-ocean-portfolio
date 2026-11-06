# Posicionamiento del Exp. 1 frente a la literatura (6 nov 2026)

El informe completo, con citas textuales y 13 afirmaciones a evitar, es un documento de trabajo local, no incluido en el repositorio. Esta nota añade la lectura propia de Kimura et al. (2014) y Spall et al. (2017). Todo lo de aquí sale de los textos en `literatura/`.

## (b) Dos plumas contra la separación

- **Cenedese & Gatto (2016a):** experimento de laboratorio con fluido de dos capas y 4 separaciones. Fuentes muy juntas derriten ≈ 0.5 de lo que derriten dos plumas independientes.
- **Kimura et al. (2014, §3f, Figs. 11–13):** MITgcm no, modelo no hidrostático de malla no estructurada.
  - Dos canales de 20 m³/s cada uno, con dx = 0, 10, 50 y 200 m.
  - **Ambiente sin estratificación** (l. 18) y corridas de 2000 s.
  - Resultado: con plumas muy juntas, la tasa de derretimiento máxima es mayor que con una sola pluma, pero el promedio sobre la cara es menor. El promedio aumenta con la separación.
  - En 2014 escriben: "There have been no attempts to study the effects of interacting meltwater plumes and their effect on the melt rate".
- **Slater et al. (2015):** de 1 a 50 canales con descarga total fija. Es distribución de la descarga, no la separación de dos plumas.
- **Carroll et al. (2015):** plumas separadas 100–400 m; se fusionan a 100 m. No reporta derretimiento.
- **Lo que aporta (b):**
  - la interacción de dos plumas **con estratificación**, de dos capas y continua, en 3D numérico, que ninguno de los anteriores tiene;
  - un barrido más denso de x0, con Q0 fijo por fuente y una referencia de pluma única;
  - la altura de despegue frente a la teoría semicónica.
- **Frase sugerida:** "extends the laboratory result of Cenedese & Gatto (2016a) and the unstratified simulations of Kimura et al. (2014) to two-layer and continuous stratification in 3-D".
- **No decir:** "first numerical study of interacting plumes" (Kimura 2014 ya lo hizo).

## (c) Pluma, viento y plataforma en un fiordo idealizado

- **Carroll et al. (2017):** casi nuestra configuración (MITgcm hidrostático, 200 m, f igual, pluma semicónica de Cowton). Muestra que el viento amplifica el intercambio de la pluma 7.5–8.6× (l. 825). No tiene oscilación de la picnoclina de plataforma.
- **Sanchez et al. (2024), coautora Straneo:** Sermilik realista, corridas con y sin pluma (NG/WG, l. 297–301). La pluma duplica el intercambio en verano. El viento y las ondas costeras quedan mezclados en NG. No reporta derretimiento.
- **Spall et al. (2017), coautores Jackson y Straneo:** MITgcm a 250 m, plano f, viento de 3 días, **sin pluma ni derretimiento**. Un evento descarga ~10 % de la capa superior.
- **Jackson et al. (2018):** ROMS sin pluma, forzamiento de plataforma (ondas de Kelvin y resonancia). Critican el forzamiento impuesto por la frontera, que es justo lo que hace (c).
- **Lo que aporta (c):**
  - corridas controladas con cada forzamiento por separado y en combinación;
  - la **partición del derretimiento** en el frente;
  - el intercambio en banda sinóptica.
- **Lo que duplica:** "el viento amplifica el intercambio de la pluma" (Carroll 2017) y "la pluma duplica el intercambio" (Sanchez 2024). Se citan como antecedentes y no se presentan como hallazgos.
- **No decir:**
  - "coastal-trapped wave" para nuestra oscilación impuesta en la OB;
  - "realistic Sermilik";
  - que el derretimiento a 200 m es robusto, mientras no haya una prueba de sensibilidad.

## Pendientes de diseño que surgieron de la lectura

1. Agregar un caso de control sin forzamiento para (c), que mide el intercambio espurio que genera la OB. Requiere aprobación de Esteban.
2. El pulso de viento de 12 h está cerca del periodo inercial (12.7 h). Se declara en el README.
3. En (b), x0c ≈ 57 m equivale a solo ~6 celdas. Las separaciones chicas quedan mal resueltas; se declara como limitación.
