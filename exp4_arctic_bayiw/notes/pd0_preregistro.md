# PD0 — Preregistro del Exp. 4 (APROBADO por Esteban el 2026-11-11 (PD0))

Fecha: 2026-11-11. No hay ninguna salida de MITgcm del Exp. 4 al escribir esto, ni del fiordo ni de la bahía. Esteban aprobó el borrador con las cuatro decisiones del §0, que ya están incorporadas. **Desde ahora los umbrales y las tolerancias están congelados.** Cambiar cualquiera de ellos exige una enmienda fechada al final de esta nota, nunca una edición silenciosa.

Los números salen de `analysis/theory_numbers.py` (`figures/theory.json`, `figures/profiles_vmf.json`) y la figura de control es `figures/profiles_vmf.{png,pdf}`. Lo que viene de la memoria y no se pudo comprobar lleva [verificar].

**Fuentes leídas para este preregistro**
- UNIS HD (Skogseth, Ellingsen et al., doi:10.21334/unis-hydrography; CC BY 4.0). El resumen de cada archivo advierte que la calibración no está revisada en todos los datos, "salinity values in particular".
- Støylen & Fer 2014 (S&F14), *NPG* 21:87–100, versión en línea: §4.1 (p. 92, configuración de MITgcm) y §4.2.1 (p. 93–94, velocidad de la onda de Kelvin interna).
- Large, McWilliams & Doney 1994 (L94), *Rev. Geophys.* 32:363–403, en la biblioteca de la tesis: ec. 28–29 (p. 373) y la elección Ric = 0.3 (§5, p. 377).
- Código y documentación de MITgcm (`/Volumes/ESTEBAN_MIT/MITgcm`, solo lectura), citados como archivo:línea.
- Frank, Jonassen, Skogseth & Vihma 2023 (F23), *JGR Atmos.*, doi:10.1029/2022JD038011, copia de acceso abierto del repositorio de la U. de Bergen: §3.1 (p. 4–5) y §3.3.2 (p. 9).
- Fer et al. 2022 (F22), *JGR Oceans*, doi:10.1029/2022JC018751, copia de acceso abierto (BORA): Key Points (p. 1), §4.3 y §5 (p. 11–17), §6 (p. 18).
- Kawaguchi et al. 2024 (K24), *Sci. Rep.*, doi:10.1038/s41598-024-66124-8, texto completo en Europe PMC (PMC11224423): introducción y discusión.
- TEOS-10 mediante gsw 3.6.21 (entorno `pyocean`).

Ningún PDF se guardó en el repositorio ni en `literatura/`.

---

## 0. Decisiones de Esteban (2026-11-11)

1. **Resolución en invierno.** Se decide en el PD1, con el costo medido. Si F2 y F3 caben a 100 m, corren a 100 m; si no, la limitación se declara (§2.3, punto 1). Ahora no se cambia la malla.
2. **Duración.** ~280 h si el cronometraje del PD1 lo permite; si no, 240 h, y se declara que con H medio solo caben 1.7 vueltas de Kelvin (§2.3, punto 2). Todos los umbrales se integran sobre la corrida completa, 0–T_fin, y T_fin se fija en el PD1. R1k frente a R1 sigue en 0–120 h. Las bandas y Δf se comprobaron para 240 h y para 280 h (§3).
3. **Ecuación de estado.** Todas las corridas del fiordo y las del puente R2/R3 usan una EOS lineal ajustada a Van Mijenfjorden (§1.3; coeficientes y procedencia en `cases.yaml`, `configs.fjord_vmf.eos`). R1 y R1k conservan la EOS de la tesis (tα = 2×10⁻⁴ K⁻¹, sβ = 7.4×10⁻⁴). Todos los números de la teoría se recalcularon con la EOS ajustada. Los valores con la EOS de la tesis se reportan como sensibilidad.
4. **El resto del borrador se aprueba tal como está.**

**Puntos que la aprobación vuelve explícitos:**
- **V2a se mide en la estación de la plataforma, no dentro del fiordo** (§4). Dentro del fiordo, la teoría predice picos lejos de f (§2.3, punto 3).
- **Para la Task 3:** la profundidad efectiva fuera de la bahía de la tesis es de 600 m, no los 1000 m de su batimetría (§4, V1).

---

## 1. Perfiles representativos (spec §3.3)

### 1.1 Criterio (en `cases.yaml`, `configs.fjord_vmf.profiles`)

- **Región:** Van Mijenfjorden interior, al este de Akseløya. El polígono se definió en UTM 33N (`polygon_utm`) y se pasó a lon/lat (`polygon_lonlat`). Su borde oeste pasa por las secciones de Akselsundet y Mariasundet de la Task 1b; su borde sur (y = 8618 km) queda al norte de Van Keulenfjorden.
  - **Control:** el polígono contiene 13 167 de las 13 169 celdas del fiordo interior de la malla (`bathy.inner_fjord_mask`) y ninguna otra celda de mar fuera de los discos de cierre de los estrechos (`theory.json`, `fjord_geometry.polygon_control`).
  - Incluye los discos de cierre de los dos estrechos, así que entran los perfiles tomados en los estrechos.
- **Meses:** verano = agosto–septiembre; invierno tardío = marzo–abril (spec §3.3).
- **Años:** 1990–2019, en las dos estaciones. Es la era de CTD electrónico (excluye las estaciones de botella de 1912 y 1926 y los datos de 1981), dentro de las tres décadas del calentamiento descrito por Skogseth et al. 2020 [verificar]. La descarga de UNIS HD termina en 2019. El efecto de la ventana es pequeño: con todos los años habría 878 y 53 perfiles en lugar de 871 y 46.
- **Control de calidad** (`gendata/profiles.qc_check`, que devuelve también el motivo de rechazo), sobre el perfil completo, en este orden:
  1. se rechazan los 68 archivos sin `PSAL` (solo temperatura) y se quitan los niveles con valores faltantes;
  2. se rechaza si quedan menos de 5 niveles;
  3. se rechaza si algún S está fuera de [0, 35.5] o algún T fuera de [−2.1, 10] °C;
  4. se rechaza si más del 10 % de los pasos de 1 m son inversiones de densidad, con σ0 de TEOS-10 (gsw) y una media móvil de 3 m. Una inversión es un paso con Δσ0 < −0.004 kg m⁻³, que es el efecto en σ0 de un error de salinidad de 0.005, del orden de la exactitud de un CTD.

  **Declaración (*post hoc*):** la tolerancia de 0.004 kg m⁻³ se eligió **después de ver los datos**. Con 0.001 kg m⁻³ se rechazaban 31 de 46 perfiles de invierno (67 %), porque la columna de invierno es casi homogénea (Δσ0 de superficie a fondo: mediana 0.057 kg m⁻³) y el ruido del sensor cuenta como inversión. Es una elección de control de calidad de los datos, no un criterio sobre el modelo, y está etiquetada *post hoc* en `cases.yaml` y en `qc_check`. La sensibilidad completa está en `theory.json` (`qc_pass_vs_inversion_tol`):

  | Tolerancia (kg m⁻³) | 0.001 | 0.002 | **0.004** | 0.01 |
  |---|---|---|---|---|
  | Perfiles de verano aceptados (de 871) | 821 | 831 | **842** | 854 |
  | Perfiles de invierno aceptados (de 46) | 15 | 22 | **36** | 44 |

- **Mínimos:** al menos 20 perfiles por estación tras el control de calidad. Un nivel entra en la mediana si lo alcanzan al menos 5 perfiles.
- **Perfil representativo** (`profiles.representative`):
  - mediana por nivel de la malla vertical del fiordo (`delR` de `cases.yaml`, 40 niveles, ruling R2);
  - cada perfil se interpola linealmente dentro de su rango de datos, y se extiende hasta la superficie si su primer dato está a menos de 5 m;
  - bajo el último nivel con al menos 5 perfiles se repite el último valor; por encima del primero, el primer valor;
  - después se impone la estabilidad con la **EOS lineal ajustada** (§1.3), con el ajuste mínimo: ajuste por mínimos cuadrados no decreciente de sβS − tαT, ponderado por el espesor de la celda (pool-adjacent-violators), aplicado solo a S.

### 1.2 Resultado

| | Verano (ago–sep) | Invierno (mar–abr) |
|---|---|---|
| Seleccionados (polígono, meses, 1990–2019) | 871 | 46 |
| Rechazados por control de calidad | 29 (5 sin salinidad, 23 por inversiones, 1 con pocos niveles) | 10 (inversiones) |
| **Aceptados (n)** | **842** (agosto 521, septiembre 321) | **36** (marzo 9, abril 27) |
| Años | 18 años, 1993–2017; 2010 aporta 344 perfiles (41 %) | 2000 (12), 2007 (1), 2008 (2), 2016 (5), 2017 (16) |
| Longitud media de los perfiles | 15.48°E | 15.14°E |
| Instrumentos | CTD (C) y algunas botellas (B) | CTD |
| Primer nivel con ≥ 5 perfiles (extensión hacia arriba) | −0.75 m (0 niveles extendidos) | −0.75 m (0 niveles extendidos) |
| Último nivel con ≥ 5 perfiles (extensión hacia abajo) | −107.9 m (8 niveles extendidos) | −98.0 m (9 niveles extendidos) |
| Niveles ajustados por estabilidad (EOS ajustada) | 2 (\|ΔS\| máx. 0.0015) | 19 (\|ΔS\| máx. 0.0046) |

- **Figura de control:** `figures/profiles_vmf.png`: T, S, N (EOS ajustada y TEOS-10) y los modos w 1–3 a H = 119 m. En el PDF, las líneas finas de los perfiles están rasterizadas.
- **Comparación (descriptiva):** S&F14 usan c_i = 0.42 m s⁻¹ para agosto de 2010 en el fiordo exterior (§4.2.1, p. 94). Nuestro c₁ de verano está entre 0.35 (H medio) y 0.52 m s⁻¹ (H máximo); ver §2.
- **Limitación:** el invierno se apoya en 36 perfiles de 5 años, y dos años (2000 y 2017) aportan el 78 %. El c₁ de invierno depende mucho de unos pocos perfiles profundos (§2, sensibilidad "un voto por año").

### 1.3 Ecuación de estado (decisión 3)

**Ajuste** (`profiles.fit_linear_eos`, `theory_numbers.py --fit-eos`):
- Mínimos cuadrados de σ0 de TEOS-10 (gsw, p = 0) contra T in situ y S práctica, en los 66 niveles por encima de 120 m de los perfiles representativos de verano e invierno juntos (T de −1.60 a 4.74 °C; S de 30.81 a 34.76).
- La forma es la EOS lineal de MITgcm, ρ = ρNil(sβ(S − sRef) − tα(T − tRef)) + ρNil − ρConst (`model/src/find_rho.F:104-107`).
- Para el ajuste, los perfiles se estabilizaron con la EOS de la tesis (cambio de S < 0.004). Al reajustar sobre los perfiles finales los coeficientes cambian en 3×10⁻⁵ relativo (`theory.json`, `eos.check_on_final_profiles`): es un punto fijo.

**Coeficientes** (en `cases.yaml`, `configs.fjord_vmf.eos`, con su procedencia):

| tα (K⁻¹) | sβ (PSU⁻¹) | tRef (°C) | sRef | ρNil (kg m⁻³) |
|---|---|---|---|---|
| 3.0456×10⁻⁵ | 8.4068×10⁻⁴ | 0.4676 | 33.6395 | 1026.96 |

Como referencia, α de TEOS-10 va de 3.3×10⁻⁵ (−1.5 °C) a 9.7×10⁻⁵ K⁻¹ (4.5 °C, S = 31). El tα de la tesis (2×10⁻⁴) es de 2 a 6 veces mayor.

**Tolerancia (fijada aquí y probada):** |σ_lineal − σ0_TEOS-10| ≤ 0.03 kg m⁻³ en el rango de los perfiles, el 1 % del contraste de densidad de la columna de verano (~3 kg m⁻³).
- Error obtenido: máx. 0.017 y rms 0.006 kg m⁻³ (`test_fit_linear_eos_recovers_teos10` y `test_real_fitted_eos_in_cases_yaml`).
- Fuera de ese rango, en combinaciones T–S que no se dan en el fiordo, el error crece (hasta 0.28 kg m⁻³ en las esquinas del rectángulo T×S), así que la EOS no se debe usar fuera de él.

**Efecto en c₁** (H = 119 m):

| | EOS ajustada | TEOS-10 | EOS de la tesis (sensibilidad) |
|---|---|---|---|
| Verano | 0.517 | 0.520 | 0.574 (+11 %) |
| Invierno | 0.189 | 0.183 | 0.185 |

La EOS ajustada reproduce TEOS-10 en c₁ dentro del 3 % en las dos estaciones.

**Puente:** R2/R3 usan la EOS ajustada con el perfil de verano. El paso R1→R2 cambia, entonces, f, el perfil **y** la EOS. Las tres cosas se declaran juntas como "el océano de Van Mijenfjorden".

---

## 2. Números de la teoría (`figures/theory.json`, EOS ajustada)

**Ecuaciones** (`theory/modes.py`):
- modos w hidrostáticos con tapa rígida, W'' + (N²/c²)W = 0 con W(0) = W(−H) = 0 (Gill 1982 §6.11 [verificar]);
- R_n = c_n/f (m);
- Bu = R₁/L y R₁/W (adimensional);
- T_K = P/c₁ (s), con P el perímetro;
- Merian: T = 4L/((2n−1)√(gH)) para una cuenca cerrada en un extremo y abierta en el otro, y T = 2L/(n√(gH)) para una cuenca cerrada en ambos (Gill 1982 cap. 5 [verificar]);
- Poincaré transversal, ω = (f² + c_n²(π/W)²)^½ (Gill 1982 §10.4 [verificar]);
- f = 2Ω sin φ.

N² se calcula con la EOS lineal en las interfaces de la malla del modelo (`modes_column`). En la sensibilidad TEOS-10 (`modes_column_teos10`) se usa la longitud media de los perfiles de cada estación (15.48 y 15.14°E). Los N² negativos de TEOS-10 se recortan a 0 y se cuentan: 0 en verano y 4 interfaces en invierno.

### 2.1 Geometría

**Fiordo interior** (malla de la Task 1b, 200 m; `bathy.inner_fjord_mask`):
- área 526.8 km²;
- profundidad máxima 119.0 m y media 62.6 m;
- longitud a lo largo del eje L = 62.8 km. Es el camino más largo por agua (8-conexo, Dijkstra) desde las celdas junto a los estrechos;
- ancho medio W = área/L = 8.39 km;
- eje principal a 12.4° al norte del este (ENE–OSO);
- **perímetro P = 176.7 km**: lazo exterior del contorno 0.5 de la máscara (marching squares, `contourpy`), incluidas las bocas. La parte de costa es 170.0 km. El perímetro en escalera (aristas de celda) da 206.2 km (+17 %) y se descartó porque exagera las costas diagonales.

**Bahía de la tesis** (`bahia_01_expand_bat.bin` y `_dx/_dy.bin`, solo lectura): rectángulo de L = 121.0 km por W = 60.6 km y H = 164 m. Costa de tres lados 2L + W = 302.6 km; lazo cerrado 363.2 km.

**Por qué P se define distinto en el fiordo y en la bahía:**
- En el fiordo, las bocas son los dos estrechos (secciones de 1.1 y 1.8 km), del orden de R₁ = 2.5–3.6 km o menores. La onda de Kelvin puede rodear la cuenca cruzándolos, y las bocas son solo el 4 % del lazo. Por eso P es el lazo cerrado.
- En la bahía, la boca mide 60.6 km, unas 9–14 R₁. Una onda de Kelvin no cruza una boca así: sale a la plataforma. Por eso P es la costa de tres lados.
- El lazo cerrado de la bahía se reporta como sensibilidad (`kelvin_circuit_period_closed_loop_h`).

### 2.2 Tabla (EOS ajustada; R1 con la EOS de la tesis)

| Caso | f (s⁻¹) / T_i (h) | H modos (m) | c₁, c₂, c₃ (m s⁻¹) | R₁ (km) | Bu = R₁/L; R₁/W | T_K = P/c₁ (h) | ω_Poincaré,1/f (modos 1, 2, 3) | Merian (h): cerrada–abierta; cerrada–cerrada |
|---|---|---|---|---|---|---|---|---|
| VMF verano | 1.4266e-4 / 12.23 | 119.0 | 0.517, 0.307, 0.189 | 3.62 | 0.058; 0.43 | 95.0 | 1.69, 1.28, 1.12 | 2.82; 1.41 (H̄ = 62.6 m, L = 62.8 km) |
| VMF verano, H medio | ″ | 62.6 | 0.352, 0.178, 0.118 | 2.47 | 0.039; 0.29 | 139.6 | 1.36, 1.10, 1.05 | ″ |
| VMF invierno | ″ | 119.0 | 0.189, 0.058, 0.039 | 1.32 | 0.021; 0.16 | 259.8 | 1.12, 1.01, 1.01 | ″ |
| VMF invierno, H medio | ″ | 62.6 | 0.049, 0.029, 0.022 | 0.34 | 0.005; 0.041 | 998 | 1.01, 1.00, 1.00 | ″ |
| Bahía de la tesis, f de la tesis, perfil "realistic", EOS de la tesis (R1) | 6.97e-5 / 25.04 | 164 | 0.480, 0.215, 0.141 | 6.89 | 0.057; 0.11 | 174.9 (costa) / 210.0 (lazo) | 1.06, 1.01, 1.01 | 3.35; 1.68 |
| Bahía de la tesis a 78°N, perfil de verano VMF (R2) | 1.4266e-4 / 12.23 | 164 | 0.617, 0.327, 0.200 | 4.32 | 0.036; 0.071 | 136.3 / 163.6 | 1.02, 1.01, 1.00 | 3.35; 1.68 |

**Sensibilidad de c₁ a la EOS** (los números del borrador, con la EOS de la tesis):
- verano 0.574 (H máx.) y 0.374 (H medio);
- invierno 0.185 y 0.063;
- R2: 0.698.

Con la EOS de la tesis, R₁ de verano sería 4.03 km y T_K 85.5 h.

**Sensibilidad "un voto por año"** (mediana de las medianas anuales, H = 119 m): verano 0.523; invierno **0.086**. El c₁ de invierno es incierto por un factor de ~2 (0.09–0.19 m s⁻¹): depende del escalón de salinidad cerca del fondo, que dan pocos perfiles.

### 2.3 Lo que estos números predicen (antes de correr)

1. **Resolución** (R₁/Δx; `theory.json`, `resolution`):

   | Δx | Verano, H máx. | Verano, H medio | Invierno, H máx. | Invierno, H medio |
   |---|---|---|---|---|
   | 100 m | 36 | 25 | 13 | 3.4 |
   | 200 m | 18 | 12 | 6.6 | 1.7 |
   | 400 m | 9.1 | 6.2 | 3.3 | 0.9 |

   **En invierno el modo 1 está mal resuelto a 200 m.** V4 solo prueba el verano. Regla de decisión 1: si el PD1 muestra que caben, F2 y F3 corren a 100 m (R₁/Δx = 3.4–13); si no, se declara la limitación.
2. **Kelvin** (vueltas en el registro):

   | | 240 h | 280 h |
   |---|---|---|
   | Verano, H máx. | 2.5 | 3.0 |
   | Verano, H medio | **1.7** | 2.0 |
   | Invierno, H máx. | 0.9 | 1.1 |
   | Invierno, H medio | 0.2 | 0.3 |

   - La spec pide al menos dos vueltas: con 280 h se cumple justo en el peor caso de verano; con 240 h se declara (decisión 2).
   - En invierno la onda de Kelvin interna no es resoluble espectralmente con ninguna de las dos duraciones.
3. **Poincaré y la banda inercial:** el fiordo es estrecho frente a πR₁ (W = 8.4 km). El primer modo de Poincaré transversal del modo vertical 1 está en 1.36–1.69 f (banda superinercial). Solo los modos verticales ≥ 3 (1.05–1.12 f) y la respuesta inercial forzada caen en la banda cercana a la inercial. En la bahía de la tesis (W = 60.6 km) el primer Poincaré está en 1.02–1.06 f, dentro de la banda cercana a la inercial.

   **Predicción del puente:** el paso de la bahía al fiordo real saca energía del modo 1 de la banda cercana a la inercial. Esto favorece la alternativa de H4a (Kelvin) por la geometría, no por f.
4. **Seiche de Merian** del fiordo interior: 1.4–2.8 h, muy por encima de f; no se mezcla con las bandas internas.

---

## 3. Bandas espectrales (`spectral_bands` en `theory_numbers.py`)

Las bandas se definen relativas a f y a la seiche de Merian más lenta del fiordo interior (cuenca cerrada en un extremo y abierta en el otro, T_M = 2.82 h). Son contiguas y no se solapan. La duración del registro, T_fin, se fija en el PD1 (decisión 2); se comprobaron las dos opciones:

| Banda | Límites | Fiordo a 78°N (cph; periodos) con 240 h | Bins con 240 h (Δf = 1/240 cph) | Bins con 280 h (Δf = 1/280 cph) | Bins con 180 h |
|---|---|---|---|---|---|
| baja (no resuelta) | [0, 2Δf) | 0–0.0083; > 120 h (con 280 h: 0–0.0071; > 140 h) | 2 | 2 | 2 |
| **Kelvin / subinercial** | [2Δf, 0.85 f) | 0.0083–0.0695; 14.4–120 h | 14.7 | 17.5 | 10.5 |
| **cercana a la inercial** | [0.85 f, 1.20 f) | 0.0695–0.0981; 10.2–14.4 h | 6.9 | 8.0 | 5.1 |
| superinercial | [1.20 f, 0.5/T_M) | 0.0981–0.1775; 5.6–10.2 h | 19.1 | 22.2 | 14.3 |
| **seiche** | [0.5/T_M, Nyquist] | 0.1775–1.5; 0.67–5.6 h | 317 | 370 | 238 |

Los límites relativos a f y a T_M no dependen de la EOS. Lo único que la EOS mueve es dónde caen los periodos de Kelvin.

**Justificación de cada límite:**
- **2Δf:** con una ventana Hann, la frecuencia más baja separable del cero.
- **[0.85, 1.20] f:** deja ≥ 5 bins en la banda cercana a la inercial incluso con 180 h, y deja un margen de ≥ 2Δf a cada lado de f para la fuga espectral. Es asimétrica hacia arriba porque las ondas cercanas a la inercial se desplazan al azul (Poincaré, vorticidad; Alford et al. 2016 [verificar]).
- **0.5/T_M:** la mitad de la frecuencia de Merian del fiordo interior. Cubre una cuenca efectiva el doble de larga (fiordo más Bellsund hasta la plataforma) y la incertidumbre de la profundidad.
- **Nyquist:** 1.5 cph, con registros cada 1200 s.

**Advertencias:**
- Los armónicos n = 7–9 del circuito de Kelvin de verano (13.6, 11.9 y 10.6 h, con H máx.) caen en la banda cercana a la inercial. Las bandas clasifican frecuencias, no mecanismos.
- Como diagnóstico **secundario y descriptivo**, prerregistrado aquí, se reporta la fracción de energía de cada banda que está a menos de 2R₁ de la costa (atrapada).
- Con H medio, el modo 1 del circuito de Kelvin de verano (139.6 h) cae en la banda "baja" con 240 h y justo dentro de la banda de Kelvin con 280 h (límite: 140 h).

**Bahía de la tesis:**
- con f de la tesis y 240 h (R1k): baja < 0.0083; Kelvin 0.0083–0.0339; cercana a la inercial 0.0339–0.0479 (3.4 bins); superinercial 0.0479–0.149; seiche 0.149–2.0 cph (registros de 900 s);
- con las 120 h de R1, la banda cercana a la inercial tiene 1.7 bins y **no es resoluble**. Por eso H4d (R1k frente a R1) no usa bandas (§5);
- a 78°N (R2/R3): las mismas bandas que el fiordo, con el límite de seiche en 0.149 cph (Merian de la bahía, 3.35 h); con 280 h son 8.0 bins en la banda cercana a la inercial.

---

## 4. Validación V1–V4 (tolerancias congeladas)

Cada diferencia contra la referencia se clasifica como (a) error nuestro, (b) física ausente, (c) elección numérica o (d) error en la referencia. Un criterio que falla no se redefine.

- **V1, la copia reproduce la tesis (R1 frente a la corrida guardada):**
  - **V1-exacto:** máx |Δη| ≤ 1×10⁻¹⁰ m en todos los registros comunes. Esto es lo esperado con el mismo código, las mismas opciones y la misma precisión; solo cambia la E/S, de MNC a MDS.
  - **V1-espectral**, si no es exacto: (i) el pico de la seiche de la bahía en η, en la estación de la boca, cae en el mismo bin de frecuencia (Δf = 1/120 cph); y (ii) la DEP en el pico difiere ≤ 10 %.
  - **Justificación:** un cambio de compilador o de opciones produce diferencias de redondeo de crecimiento lento. Un error de configuración (f, perfil, forzamiento) cambia la frecuencia o la amplitud mucho más que un 10 %.
  - **De paso:** el f0 del `data` es 6.97×10⁻⁵ s⁻¹, con periodo inercial de 25.04 h. Los scripts de la tesis usan ~31 h. V1 lo resuelve con el espectro de R1.
  - **Nota para la Task 3 (profundidad de la columna):** la suma del `delR` de la tesis es 599.9 m, pero la plataforma de su batimetría mide 1000 m. MITgcm trunca la columna en el fondo de la malla:
    - `model/src/ini_masks_etc.F:112-113` acota la fracción de cada celda a [0, 1], y solo hay celdas k = 1…Nr;
    - `ini_masks_etc.F:126-141` recalcula el fondo como R_low = rF(1) − Σ drF·hFacC, que no puede bajar de 599.9 m;
    - `model/src/ini_depths.F:89-98` es la rama que pone R_low = rF(Nr+1) (el fondo de la malla) cuando no se da `bathyFile`; con `bathyFile`, R_low se lee en la línea 119 y lo trunca `ini_masks_etc.F`.

    La profundidad efectiva fuera de la bahía es, entonces, de 600 m. No afecta a V1 (misma configuración), pero la Task 3 lo declara (§0).
- **V2, frecuencias contra la teoría (F1):**
  - **V2a:** el pico de energía cinética horizontal dentro de la banda cercana a la inercial, **en la estación de la plataforma** (W ≫ R₁, donde la oscilación inercial no está restringida), está a ≤ 1Δf de f a 78°N. Δf = 1/T_fin, o el de la ventana analizada si es más corta.
    - Dentro del fiordo no se exige f: por §2.3, la teoría predice picos de 1.05–1.12 f (modos verticales altos) y de 1.36–1.69 f para el modo 1. Se reportan como prueba descriptiva de la teoría.
  - **V2b:** la velocidad de fase de la onda de Kelvin interna en F1 está dentro de [0.8 c₁(H̄), 1.2 c₁(H_máx)].
    - **Con la EOS ajustada: [0.28, 0.62] m s⁻¹.** Antes, con la EOS de la tesis, era [0.30, 0.69] m s⁻¹. La fórmula no cambia; solo se recalcula c₁ (`theory.json`, `v2b_range_m_s-1`).
    - Se mide por el tiempo de llegada del frente de desplazamiento de la picnoclina a ≥ 3 estaciones costeras sucesivas del fiordo interior, sobre la distancia a lo largo de la costa.
    - El rango cubre la incertidumbre de H del modelo plano, más ±20 % por la pendiente y la no linealidad. S&F14 obtuvieron 0.40–0.45 m s⁻¹ en su modelo, frente a 0.42 de la teoría de dos capas (p. 94), un acuerdo de ~7 %.
- **V3, M2 sola (contra S&F14):** pasa si se cumplen las tres:
  1. **presencia:** desplazamiento de la picnoclina a frecuencia M2 atrapado en la costa, con escala de decaimiento hacia fuera de la costa de 1–8 km (S&F14 ven un radio de 2–3 km, p. 93; R₁ = 2.5–3.6 km);
  2. **sentido ciclónico:** costa a la derecha, hacia el este por la costa sur del fiordo interior (S&F14 §4.2 y §4.3);
  3. **velocidad de fase** dentro de un factor 2 de 0.42 m s⁻¹, es decir 0.21–0.84 m s⁻¹. Es más estricto que el "orden de magnitud" de la spec.
- **V4, resolución (F1 a 400 m frente a 200 m):**
  - **V4a:** el flujo baroclínico neto hacia fuera, integrado en 0–T_fin y sumado sobre los dos estrechos, difiere ≤ 25 %.
  - **V4b:** la fracción del flujo de cada banda resuelta (Kelvin, cercana a la inercial, superinercial, seiche) difiere ≤ 0.10 en valor absoluto.
  - El flujo de cada estrecho se reporta por separado y es descriptivo. A 400 m, Mariasundet (0.5–0.7 km) mide 1–2 celdas, así que se espera que falle solo.
  - **Justificación:** R₁/Δx baja de 18 a 9 en verano. El 25 % es la tolerancia más estricta que se puede pedir sin exigir convergencia a una resolución que la Task 1b ya mostró que tiene errores de ancho de ~1 celda en los estrechos.
  - Si falla, se declara la limitación (spec §6).

---

## 5. Hipótesis H4a–H4d (umbrales congelados)

**Regla general:** un umbral de cociente debe superar la incertidumbre de resolución (V4a, 25 %). Por eso el **factor mínimo es 1.5**, el que se usa en H4c y H4d: un cociente entre 1/1.5 y 1.5 es "no distinguible". Todos los flujos son baroclínicos, por el método de Nash et al. 2005 [verificar] como en la tesis, y se integran sobre la corrida completa, 0–T_fin (T_fin fijado en el PD1), salvo R1k frente a R1, que va en 0–120 h. Las corridas que se comparan usan el mismo pulso, la misma duración, el mismo Δt y la misma EOS.

- **H4a, partición (F1):** φ_b es la fracción del flujo baroclínico neto hacia fuera por los dos estrechos que lleva la banda b, sobre la suma de las cuatro bandas resueltas. La banda baja se reporta aparte.
  - **Apoyada** (como la tesis) si φ_NI > 0.5.
  - **Alternativa** (domina Kelvin) si φ_K > 0.5.
  - Si no, es **mixta**.
  - Si la banda baja lleva más del 25 % del flujo total en valor absoluto, H4a se marca "no concluyente con T_fin" (§2.3, punto 2).
  - La spec dice que se acepta el resultado que salga; la predicción de §2.3, punto 3, favorece la alternativa.
- **H4b, hielo (F2 frente a F3):**
  - **(i)** Ē_in(F2)/Ē_in(F3) ≤ 0.5, con Ē_in la energía baroclínica (cinética más potencial disponible) del fiordo interior promediada en 0–T_fin, sin la banda baja.
    - **Relación con la regla del factor 1.5:** 0.5 equivale a un factor 2, es decir 2/1.5 = 1.33 veces el mínimo, un margen extra de 33 %. Cubre dos fuentes de incertidumbre propias de H4b que V4 no mide: (a) la resolución de invierno (R₁/Δx = 1.7–6.6 a 200 m, frente a ≥ 12 en el verano que prueba V4); y (b) que Ē_in no sale de un flujo por una sección sino del balance del fiordo, cuyo residuo se reporta.
    - Bajo hielo fijo el viento no trabaja dentro del fiordo, así que se espera una diferencia mucho mayor que un factor 2.
  - **(ii)** "entra por los estrechos": F_bc,in ≥ |F_bt,in|, con los flujos netos hacia dentro por los dos estrechos, integrados en 0–T_fin.
    - **Por qué este criterio:** en F2 el forzamiento local dentro del fiordo es nulo por construcción. Las fuentes de energía interna adentro son el flujo baroclínico que entra y la conversión, dentro del fiordo, del flujo barotrópico que entra. Si F_bc,in ≥ |F_bt,in|, más de la mitad llega ya como onda interna, aunque todo el barotrópico se convirtiera.
  - Se reportan (i) y (ii) por separado.
- **H4c, estación (F1 frente a F3, agua abierta):**
  - **(i)** W_NI(F1)/W_NI(F3) ≥ 1.5, con W_NI = ∫∫ τ·u_NI dA dt el trabajo del viento sobre la velocidad superficial filtrada en la banda cercana a la inercial, en el fiordo interior;
  - **(ii)** el cociente del flujo baroclínico radiado hacia fuera por los estrechos (todas las bandas resueltas), F1/F3, es ≥ 1.5.
  - **Apoyada** si se cumplen las dos; **refutada** si las dos son ≤ 1/1.5; **no distinguible** en otro caso.
  - **Advertencia:** el modo 1 de invierno está mal resuelto (§2.3, punto 1). Si H4c sale apoyada, puede deberse en parte a la resolución, y se declara (salvo que F3 corra a 100 m; decisión 1).
- **H4d, cierre:** η = F_bc,out/W_viento, la fracción del trabajo del viento que sale como flujo baroclínico (por la boca de la bahía en R1/R1k y por los estrechos en F1/F1c).
  - **Apoyada** en cada par por separado si η(KPP)/η(cierre de la tesis) ≥ 1.5; **refutada** si ≤ 1/1.5; **no distinguible** en otro caso.
  - Pares: R1k frente a R1 (0–120 h) y F1 frente a F1c (0–T_fin). No usa bandas (§3).
  - La tesis radiaba 0.3 % y 0.08 % (bitácora del 11 de noviembre); el factor 1.5 es deliberadamente modesto.
- **Puente:** descriptivo, sin umbral (spec §7).

---

## 6. Parámetros de mezcla y numéricos (aprobados)

| Parámetro | Valor | Fuente |
|---|---|---|
| EOS | lineal, ajustada a VMF (§1.3) en el fiordo y en R2/R3; la de la tesis en R1/R1k | decisión 3; `cases.yaml` `configs.fjord_vmf.eos`; `model/src/find_rho.F:104-107` |
| Viscosidad vertical de fondo `viscAr` (KPP la usa como fondo: `viscArNr(1)` en `pkg/kpp/kpp_routines.F:1197`) | 1.0×10⁻⁴ m² s⁻¹ | L94 ec. 29 (p. 373), ν_w = 1.0×10⁻⁴ m² s⁻¹ (onda interna) |
| Difusividad vertical de fondo `diffKrT = diffKrS` (fondo de KPP: `CALC_3D_DIFFUSIVITY`, `pkg/kpp/kpp_calc.F:515-526`; suma en `kpp_routines.F:1198-1199`) | 1.0×10⁻⁵ m² s⁻¹ | L94 ec. 29, κ_w = 0.1×10⁻⁴ m² s⁻¹; es el mismo valor de difusividad que la tesis (`data.stage1`), así que F1 y F1c solo difieren en la capa límite de KPP, la mezcla por cizalla y la viscosidad de fondo |
| Ri crítico de la capa límite `Ricr` | 0.3 (valor por omisión) | `pkg/kpp/kpp_readparms.F:112`; L94 §5, p. 377 ("Ric = 0.3") |
| Mezcla interior por cizalla: `Riinfty`, `difm0 = difs0 = dift0` | 0.7; 5×10⁻³ m² s⁻¹ (valores por omisión) | `kpp_readparms.F:129, 132-134`; L94 ec. 28 (p. 373): ν₀ = 50×10⁻⁴ m² s⁻¹, Ri₀ = 0.7, p₁ = 3 |
| Convección `difmcon = difscon = diftcon`, `BVSQcon` | 0.1 m² s⁻¹; −0.2×10⁻⁴ s⁻² (valores por omisión) | `kpp_readparms.F:130, 136-138` |
| Resto de KPP | valores por omisión de MITgcm; `KPP_OPTIONS.h` sin cambios (con `KPP_GHAT`, sin `EXCLUDE_KPP_SHEAR_MIX`) | `pkg/kpp/KPP_OPTIONS.h:62, 65` |
| Viscosidad horizontal | **Smagorinsky** `viscC2smag = 2.2` y biarmónica `viscC4smag = 1` | **S&F14 §4.1 (p. 92): la misma configuración de MITgcm en Van Mijenfjorden** ("Smagorinsky type with value 2.2 along with a small biharmonic viscosity factor (viscC4smag = 1) as suggested by Griffies and Hallberg (2000)"). La documentación de MITgcm da 2.2–4 para el océano (`doc/algorithm/algorithm.rst:2449-2451`) |
| Difusividad horizontal de T y S | 0 (la difusión numérica la pone el limitador) | F1c usa la misma, para aislar el cierre vertical |
| Arrastre de fondo `bottomDragQuadratic` | 2.5×10⁻³ | S&F14 §4.1, p. 92 ("quadratic bottom drag coefficient of 0.0025"); MITgcm: rango típico 0.001–0.003 (`doc/algorithm/algorithm.rst:1636-1638`); el valor por omisión es 0 (`model/src/set_defaults.F:137`) |
| Advección de T y S | `tempAdvScheme = saltAdvScheme = 77` (limitador de flujo de 2.º orden, Superbee) | `doc/algorithm/adv-schemes.rst:237-239, 455`. Se declara la mezcla numérica (Klingbeil et al. 2014, en la biblioteca de la tesis). El esquema 33 (DST3 con limitador) es menos compresivo; queda como alternativa si el PD1 muestra escalones en la picnoclina, y usarla exigiría una enmienda fechada |

**Leith o Smagorinsky.** Se usa Smagorinsky y no Leith. Leith responde solo a la vorticidad ("can only respond to buildup of vorticity at the grid scale", `doc/algorithm/algorithm.rst:2491, 2517-2518`); por eso la documentación añade `viscC2LeithD` para los modos divergentes (`algorithm.rst:2523-2545`), y las ondas internas son divergentes. Además, Smagorinsky 2.2 con `viscC4smag = 1` es la configuración con la que S&F14 reprodujeron las ondas de Kelvin internas en este mismo fiordo con MITgcm, lo que hace que V3 sea una comparación más limpia. La estabilidad se verifica en el cronometraje del PD1.

**Nota de S&F14 (descriptiva):** usaron Az = 1×10⁻³ m² s⁻¹, Δx = 100 m y 32 niveles (§4.1, p. 92). No se copia su Az, porque KPP lo calcula.

---

## 7. Viento (aprobado)

- **Dirección:** a lo largo del eje del fiordo interior, **hacia la boca** (down-fjord): viento **del ENE**, 078° de la cuadrícula UTM, hacia 258°. La convergencia de cuadrícula UTM 33 a ~16°E es ~1°, despreciable.
  - El eje es la componente principal de las celdas del fiordo interior: 12.4° al norte del este (`theory.json`, `wind`).
  - Vector de esfuerzo para el pulso de 0.1 N m⁻²: (τx, τy) = 0.1 × (−0.977, −0.215) N m⁻².
- **Justificación (F23, leído):**
  1. "near-surface winds in mountainous areas like Svalbard frequently get channeled along the valley and fjord axes", con diferencias de hasta 7 m s⁻¹ entre fiordos vecinos (§3.3.2, p. 9);
  2. los vientos sinópticos fuertes en invierno y primavera vienen de los sectores sur y este (§3.1, p. 4);
  3. "similar patterns can be expected for other locations in Svalbard" (p. 9).

  **Inferencia, no dato:** F23 estudia Isfjorden, no Van Mijenfjorden. Que el viento fuerte canalizado en Van Mijenfjorden sea del ENE (hacia la boca) es una extrapolación de su resultado, coherente con el plan ("de E a O, salvo que la climatología diga otra cosa").
- **Caveat de verano:** en JJA dominan los oestes sinópticos (16 % frente a 8 % de los estes; F23 §3.1, p. 5), con los vientos más débiles del año. Se usa la misma dirección en todas las F, para que los contrastes estacionales (F1 frente a F3) no mezclen dirección y estación. En teoría lineal la respuesta al pulso opuesto cambia de signo, no de amplitud. Las asimetrías no lineales no se estudian.
- **F5 (transversal):** el esfuerzo se gira 90° en sentido antihorario: viento **del NNO** (348°), hacia 168°. Es la dirección sinóptica más frecuente en MAM: "northerly directions dominating the total distribution (occurrence frequency 21%)" (F23 §3.1, p. 4–5).

---

## 8. Hielo: coeficiente de arrastre hielo–agua (aprobado)

- **Valor:** `SEAICE_waterDrag = 5.5×10⁻³`, el valor por omisión de MITgcm (`pkg/seaice/seaice_readparms.F:401`; usado en `pkg/seaice/lsr.F:221-223`). El mismo para el hielo fijo (F2) y a la deriva (F4).
- **Respaldo (leído):**
  - F22: "A representative range of drag coefficient for the MOSAiC sampling site was estimated to (4–6) × 10⁻³, which is a typical value for Arctic floe observations" (Key Points y resumen, p. 1; con la corriente relativa a 22 m, §4.3, p. 11–13). En su discusión (p. 17) recopilan otros valores: McPhee 2002 (SHEBA), 2–8×10⁻³ según el mes; Lu et al. 2011, 1–8×10⁻³ para témpanos lisos y 20–22×10⁻³ para rugosos; Cole et al. 2017, 1–4×10⁻³.
  - K24 (introducción): Shirasawa & Ingram (1997) proponen C_w = 5×10⁻³ a 1 m bajo el hielo en Barrow Strait. Es el análogo más cercano: **hielo fijo**, referencia somera, como nuestra celda superficial de 1.5 m. K24 mide C_w de 1–5×10⁻³ con capa estable (estaciones en fusión) y ≫ 10×10⁻³ con μ ≤ 0 (discusión), y advierte que un C_w constante es una simplificación.
- **Limitaciones:**
  - la referencia del modelo es la velocidad de la primera celda (centro a 0.75 m), no la de 22 m de F22;
  - el hielo de Van Mijenfjorden es de primer año y liso (Høyland 2009 [verificar]), así que lo esperable está en la parte baja del rango de F22;
  - no hay corrida de sensibilidad a C_w en el presupuesto, y se declara.
- **Verificación pendiente (PD1, Task 4):** qué esfuerzo recibe el océano bajo hielo inmóvil, τ = ρ₀ C_w |u_hielo − u|(u_hielo − u) con u_hielo = 0, y si KPP lo usa como u*. Se lee en el código y se cita en `notes/pd1_seaice.md`.

---

## 9. Lo que se decide después (PD1) y lo que queda declarado

1. **Malla de invierno** (decisión 1): F2 y F3 corren a 100 m si caben según el costo medido; si no, la limitación de R₁/Δx se declara.
2. **T_fin** (decisión 2): 280 h o 240 h. Con 240 h, 1.7 vueltas de Kelvin (H medio), declaradas.
3. **V2a se mide en la estación de la plataforma.** Dentro del fiordo los picos se reportan sin criterio (§2.3, punto 3).
4. **Task 3:** la profundidad efectiva fuera de la bahía de la tesis es de 600 m y no de 1000 m, por la suma del `delR` (§4, V1, con las citas de `ini_masks_etc.F` e `ini_depths.F`). Se declara en `notes/pd1_puente.md`.
5. **Task 3/4:** el `data` del fiordo y de R2/R3 usa `eosType = 'LINEAR'` con tAlpha, sBeta, rhoNil y los valores de referencia tRef y sRef de `cases.yaml` (`configs.fjord_vmf.eos`). R1/R1k conservan los valores de la tesis.
6. **Declarado *post hoc*:** la tolerancia de inversiones de 0.004 kg m⁻³ (§1.1).

## Enmiendas

### Enmienda 1 (2026-11-12)

**Qué cambia.**

1. Todos los flujos de H4a–H4d y el balance de energía usan la descomposición hidrostática de Kang y Fringer (2012, ecs. 13–17), la misma de la tesis en `Output_IWBay/balance_kang_fringer.py` (copia literal en `thesis_copy/Output_IWBay/`, con su sha256 en `thesis_copy/SOURCES.md`):
   - ū es la media vertical **instantánea** de la columna, u′ = u − ū;
   - ρ′ = ρ − ρ_b(z), con ρ_b el perfil inicial de reposo, el mismo en todos los términos;
   - el flujo por una cara es el trabajo de presión F (F_bt = H ū_n p_bt, F_bc = Σ u_n′ p′ dz) **más** la advección de energía A (ecs. 15–16), evaluados en las caras de la malla; el flujo de cada estrecho es la suma sobre sus caras;
   - la APE es la exacta (ec. 12); la APE lineal ½g²ρ′²/(ρ₀N_b²) se reporta solo como comparación;
   - C = Σ ρ′ g W_bt dV, con W_bt = −∇·((H + z) ū) en la forma discreta de la malla C;
   - P_bt = Σ τ·ū dA y P_bc = Σ τ·(u_sup − ū) dA;
   - D es el residuo; en el fiordo incluye el arrastre de fondo, Smagorinsky y KPP.

   Para el fiordo, la caja es el fiordo interior cerrado por secciones en escalera sobre las caras de la malla en Akselsundet y Mariasundet (`analysis/budget.py`, `box_from_sections`, con las líneas de sección de `figures/bathy_report.json`).
2. **H4d:** el umbral (factor 1.5, sin cambio) se evalúa con W_viento = P_bt + P_bc, es decir, η = F_bc,out / (P_bt + P_bc), con F_bc,out el flujo baroclínico de Kang y Fringer (F_bc + A_bc) que sale. El cociente con P_bc en el denominador, η′ = F_bc,out / P_bc, se reporta solo como descriptivo.
3. **Ningún umbral ni tolerancia cambia.** Los umbrales de H4a–H4d, V1–V4 y la regla del factor 1.5 quedan como están en §4–§5.
4. El flujo de Nash et al. (2005) de `energy_flux_IWBay.py` se conserva como contraste con la tesis: Task 5, paso 5, sección j = 323, i = 274–285, < 1 %. Fuera de eso es descriptivo.

**Por qué.** La tesis tiene dos definiciones de flujo: Nash et al. (2005), en `energy_flux_IWBay.py`, con medias temporales de ρ y u en la sección; y Kang y Fringer (2012), en `balance_kang_fringer.py`, con ū instantánea y ρ_b inicial. Solo la segunda entra en un balance que cierra, porque ρ′ y la separación barotrópica–baroclínica son las mismas en E, F, C y P. Un flujo de Nash, con otra referencia, no cierra el balance (`fondo_kang_fringer.py`). H4b(i) usa la energía del fiordo con su residuo, y H4d compara un flujo con el trabajo del viento. Por eso los dos necesitan las definiciones del balance. Lo decidió Esteban el 2026-11-12.

**Cuándo.** La decisión se tomó el 2026-11-12, antes de correr la campaña del fiordo (F1–F5, F1c, F1g, V3, V4 y R2/R3). R1 y R1k ya habían corrido (Task 3), pero no se había calculado ningún flujo sobre ellos. Esta enmienda se escribe en la Task 5. En esa misma tarea se calcula sobre R1 el contraste del método con los scripts de la tesis (`figures/flux_vs_thesis.json`, `figures/budget_vs_thesis.json`). Ese contraste reproduce números que ya estaban en la tesis, porque R1 es idéntico bit a bit a la corrida de la tesis. No se calculó H4d sobre R1/R1k: eso es de la Task 7.

### Enmienda 2 (2026-11-12)

**Qué cambia.**
- **Numerador de H4d.** F_bc,out es solo el trabajo de presión baroclínico de Kang y Fringer que sale de la caja: F_bc = Σ u_n′ p′ dz sobre las caras abiertas, hacia fuera. Es la contraparte de Kang y Fringer del flujo de Nash que pedía la H4d original. En el código es `F_bc_pressure_out` (`analysis/budget.py`, `kf_budget` y `f_bc_pressure_out`).
- **Cociente.** η = F_bc,out / (P_bt + P_bc); η′ = F_bc,out / P_bc sigue siendo descriptivo.
- **Advección.** La suma F_bc + A_bc (trabajo de presión más advección) se reporta solo como descriptiva.
- **Umbrales.** No cambia ninguno.

**Por qué.** Corrige el punto 2 de la Enmienda 1. El "(F_bc + A_bc)" que aparece ahí venía de una paráfrasis, no de una decisión explícita de Esteban. Esteban lo decidió el 2026-11-12. El texto de la Enmienda 1 no se edita.

**Qué se sabía al decidir.** El reparto en R1 entre 0 y 120 h ya se conocía por el contraste de la Task 5 (`figures/budget_vs_thesis.json`): el trabajo de presión baroclínico saca 0.28 GJ de la bahía y la advección A_bc mete 5.68 GJ. En cambio, no se había calculado ningún cociente de H4d en ninguna corrida, ni en R1 ni en R1k, y no se habían analizado las salidas de la campaña del fiordo.

### Enmienda 3 (2026-11-13)

**Qué cambia.**
- **Flujo que decide H4a, H4b(ii) y H4c(ii).** Es el trabajo de presión baroclínico de Kang y Fringer solo: F_bc = Σ u_n′ p′ dz sobre las caras de los estrechos (el mismo de `F_bc_pressure_out`, `analysis/budget.py`). La suma con la advección, F_bc + A_bc, pasa a ser descriptiva en las tres.
- **H4b(ii).** La contraparte barotrópica es el trabajo de presión barotrópico, F_bt = H ū_n p_bt. La comparación F_bc,in ≥ |F_bt,in| es, entonces, entre trabajos de presión. F_bt + A_bt es descriptivo.
- **H4d** queda como la dejó la Enmienda 2.
- **Umbrales.** No cambia ninguno: los de H4a–H4d, V1–V4 y la regla del factor 1.5 quedan como están en §4–§5.

**Por qué.** Hace que todos los flujos de las hipótesis sean la misma magnitud:
- el numerador de H4d (Enmienda 2);
- el flujo de V4, que se eligió como trabajo de presión porque V4a acota los cocientes de H4 (ruling R23 de la Task 6);
- la intención original del PD0, que pedía el flujo de Nash, un trabajo de presión.

Lo decidió Esteban el 2026-11-13. El texto de las Enmiendas 1 y 2 no se edita.

**Qué se sabía al decidir.**
- Los descriptivos de V4 y del PD2 (`figures/pd2_v4_diagnostics.json`, `pd2_v4_descriptive.json`), que son de trabajo de presión: en F1 el flujo neto por los estrechos es casi entero de la banda baja.
- El reparto F/A de R1 en 0–120 h (Task 5, `figures/budget_vs_thesis.json`): el trabajo de presión baroclínico saca 0.28 GJ de la bahía y la advección A_bc mete 5.68 GJ.
- **No se había calculado ningún número de H4a–H4d en ninguna corrida**, ni ninguna partición por bandas de flujo del fiordo para H4. Al decidir, la Task 7 estaba escribiendo las series de F1 en caché, sin ensamblarlas ni mirarlas. Los únicos números reales vistos en la Task 7 eran los de una corrida corta de cronometraje del PD1 (`F1_pd1final_48h`, 48 h): diferencias entre dos rutas de lectura (cero) y la energía de la caja por bandas, sin flujos.

#### Nota a la Enmienda 3 (2026-11-13)

La Enmienda 3 dice que, al decidirla, no se había calculado "ninguna partición por bandas de flujo del fiordo para H4". **Es incompleto.** `figures/pd2_v4_diagnostics.json` (Task 6, parte B, versionado el 2026-11-12) ya tenía, para F1, la partición por bandas del trabajo de presión baroclínico de Kang y Fringer en cada estrecho, la misma magnitud y el mismo método que la Enmienda 3 hace decisivos para H4a:
- Akselsundet: banda baja −0.997 GJ; Kelvin +0.052; cercana a la inercial −0.034; superinercial +0.011; seiche −0.009 (+ hacia fuera);
- Mariasundet: −0.147, −0.034, −0.009, −0.001 y −0.005 GJ.

La tabla del PD2 de `notes/pd2_validacion.md` los muestra. Se calcularon para V4, no para H4. Con ellos, la fracción de la banda baja de H4a y el signo de la suma de las bandas resueltas de F1 ya se podían deducir cuando se escribió la Enmienda 3.

Lo demás de la Enmienda 3 es correcto: no se había ensamblado ni evaluado ningún número de H4 de la Task 7. El texto de la Enmienda 3 no se edita.
