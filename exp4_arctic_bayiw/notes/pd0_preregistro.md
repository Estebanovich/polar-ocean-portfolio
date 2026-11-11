# PD0 — Preregistro del Exp. 4 (BORRADOR — pendiente de aprobación de Esteban (PD0))

Fecha: 2026-11-11. No hay ninguna salida de MITgcm del Exp. 4 al escribir esto: ni del fiordo ni de la bahía. Todo lo que sigue son **propuestas**. Esteban las aprueba o las cambia en el PD0. Hasta entonces no corre nada del fiordo (F*, V3, V4); R1 puede correr porque no depende del PD0 (ruling R3).

Los números salen de `analysis/theory_numbers.py` (`figures/theory.json`, `figures/profiles_vmf.json`) y la figura de control es `figures/profiles_vmf.{png,pdf}`. Lo que viene de la memoria y no se pudo comprobar lleva [verificar].

**Fuentes leídas para este borrador**
- UNIS HD (Skogseth, Ellingsen et al., doi:10.21334/unis-hydrography; CC BY 4.0). El resumen de cada archivo advierte que la calibración no está revisada en todos los datos, "salinity values in particular".
- Støylen & Fer 2014 (S&F14), *NPG* 21:87–100, versión en línea: §4.1 (p. 92, configuración de MITgcm) y §4.2.1 (p. 93–94, velocidad de la onda de Kelvin interna).
- Large, McWilliams & Doney 1994 (L94), *Rev. Geophys.* 32:363–403, en la biblioteca de la tesis: ec. 28–29 (p. 373) y la elección Ric = 0.3 (§5, p. 377).
- Código y documentación de MITgcm (`/Volumes/ESTEBAN_MIT/MITgcm`, solo lectura), citados como archivo:línea.
- Frank, Jonassen, Skogseth & Vihma 2023 (F23), *JGR Atmos.*, doi:10.1029/2022JD038011, copia de acceso abierto del repositorio de la U. de Bergen: §3.1 (p. 4–5) y §3.3.2 (p. 9).
- Fer et al. 2022 (F22), *JGR Oceans*, doi:10.1029/2022JC018751, copia de acceso abierto (BORA): Key Points (p. 1), §4.3 y §5 (p. 11–17), §6 (p. 18).
- Kawaguchi et al. 2024 (K24), *Sci. Rep.*, doi:10.1038/s41598-024-66124-8, texto completo en Europe PMC (PMC11224423): introducción y discusión.

Ningún PDF se guardó en el repositorio ni en `literatura/`.

---

## 1. Perfiles representativos (spec §3.3)

### 1.1 Criterio (propuesta; en `cases.yaml`, `configs.fjord_vmf.profiles`)

- **Región:** Van Mijenfjorden interior, al este de Akseløya. El polígono se definió en UTM 33N (`polygon_utm`) y se pasó a lon/lat (`polygon_lonlat`). Su borde oeste pasa por las secciones de Akselsundet y Mariasundet de la Task 1b; su borde sur (y = 8618 km) queda al norte de Van Keulenfjorden.
  - **Control:** el polígono contiene 13 167 de las 13 169 celdas del fiordo interior de la malla (`bathy.inner_fjord_mask`) y ninguna otra celda de mar fuera de los discos de cierre de los estrechos (`theory.json`, `fjord_geometry.polygon_control`).
  - Incluye los discos de cierre de los dos estrechos, así que entran los perfiles tomados en los estrechos.
- **Meses:** verano = agosto–septiembre; invierno tardío = marzo–abril (spec §3.3).
- **Años:** 1990–2019, en las dos estaciones. Es la era de CTD electrónico (excluye las estaciones de botella de 1912 y 1926 y los datos de 1981), dentro de las tres décadas del calentamiento descrito por Skogseth et al. 2020 [verificar]. La descarga de UNIS HD termina en 2019. El efecto de la ventana es pequeño: con todos los años habría 878 y 53 perfiles en lugar de 871 y 46.
- **Control de calidad** (`gendata/profiles.qc`), sobre el perfil completo:
  1. se quitan los niveles con valores faltantes; los 68 archivos sin `PSAL` (solo temperatura) se rechazan;
  2. se rechaza si quedan menos de 5 niveles;
  3. se rechaza si algún S está fuera de [0, 35.5] o algún T fuera de [−2.1, 10] °C;
  4. se rechaza si más del 10 % de los pasos de 1 m son inversiones de densidad, con σ0 de TEOS-10 (gsw) y una media móvil de 3 m. Una inversión es un paso con Δσ0 < −0.004 kg m⁻³, que es el efecto en σ0 de un error de salinidad de 0.005, del orden de la exactitud de un CTD.

  **Declaración:** la tolerancia de 0.004 kg m⁻³ se eligió **después de ver los datos**. Con 0.001 kg m⁻³ se rechazaban 31 de 46 perfiles de invierno (67 %), porque la columna de invierno es casi homogénea (Δσ0 de superficie a fondo: mediana 0.057 kg m⁻³) y el ruido del sensor cuenta como inversión. Es una elección de control de calidad de los datos, no un criterio sobre el modelo. La sensibilidad completa está en `theory.json` (`qc_pass_vs_inversion_tol`):

  | Tolerancia (kg m⁻³) | 0.001 | 0.002 | **0.004** | 0.01 |
  |---|---|---|---|---|
  | Perfiles de verano aceptados (de 871) | 821 | 831 | **842** | 854 |
  | Perfiles de invierno aceptados (de 46) | 15 | 22 | **36** | 44 |

- **Mínimos:** al menos 20 perfiles por estación tras el control de calidad. Un nivel entra en la mediana si lo alcanzan al menos 5 perfiles.
- **Perfil representativo** (`profiles.representative`):
  - mediana por nivel de la malla vertical del fiordo (`delR` de `cases.yaml`, 40 niveles, ruling R2);
  - cada perfil se interpola linealmente dentro de su rango de datos, y se extiende hasta la superficie si su primer dato está a menos de 5 m;
  - bajo el último nivel con al menos 5 perfiles se repite el último valor;
  - después se impone la estabilidad con la **EOS lineal del modelo**, con el ajuste mínimo: ajuste por mínimos cuadrados no decreciente de sβS − tαT, ponderado por el espesor de la celda (pool-adjacent-violators), aplicado solo a S.

### 1.2 Resultado

| | Verano (ago–sep) | Invierno (mar–abr) |
|---|---|---|
| Seleccionados (polígono, meses, 1990–2019) | 871 | 46 |
| Rechazados por control de calidad | 29 (5 sin salinidad, 23 por inversiones, 1 corto) | 10 (inversiones) |
| **Aceptados (n)** | **842** (agosto 521, septiembre 321) | **36** (marzo 9, abril 27) |
| Años | 18 años, 1993–2017; 2010 aporta 344 perfiles (41 %) | 2000 (12), 2007 (1), 2008 (2), 2016 (5), 2017 (16) |
| Instrumentos | CTD (C) y algunas botellas (B) | CTD |
| Último nivel con ≥ 5 perfiles | −107.9 m (se extienden 8 niveles) | −98.0 m (se extienden 9 niveles) |
| Niveles ajustados por estabilidad (EOS lineal) | 0 | 16 (\|ΔS\| máx. 0.0039) |

- **Figura de control:** `figures/profiles_vmf.png`: T, S, N (EOS del modelo y TEOS-10) y los modos w 1–3 a H = 119 m.
- **Comparación (descriptiva):** S&F14 usan c_i = 0.42 m s⁻¹ para agosto de 2010 en el fiordo exterior (§4.2.1, p. 94). Nuestro c₁ de verano está entre 0.37 (H medio) y 0.57 m s⁻¹ (H máximo); ver §2.
- **Limitación:** el invierno se apoya en 36 perfiles de 5 años, y dos años (2000 y 2017) aportan el 78 %. El c₁ de invierno depende mucho de unos pocos perfiles profundos (§2, sensibilidad "un voto por año").

### 1.3 La EOS lineal (decisión del plan, declarada)

El plan usa la EOS de la tesis, tα = 2×10⁻⁴ K⁻¹ y sβ = 7.4×10⁻⁴ (`thesis_copy/.../data.stage1:15-16`; `model/src/find_rho.F:104-107`), para que el puente sea limpio. A las temperaturas de Van Mijenfjorden, α de TEOS-10 es de 3.3×10⁻⁵ (−1.5 °C, S = 34.5) a 9.7×10⁻⁵ K⁻¹ (4.5 °C, S = 31) (gsw.alpha). La EOS del modelo exagera, entonces, el efecto de la temperatura de 2 a 6 veces.

- **Efecto medido:**
  - verano: c₁ = 0.574 m s⁻¹ con la EOS del modelo y 0.520 con TEOS-10 (+10 %);
  - invierno: 0.185 y 0.186 (la salinidad domina).
- **Propuesta:** se mantiene la EOS de la tesis y en el README se declara el sesgo de +10 % de c₁ en verano.
- **Alternativa para Esteban:** ajustar α y β a los perfiles de Van Mijenfjorden. Rompería el principio de que el puente R1→R2 solo cambia f y el perfil.

---

## 2. Números de la teoría (`figures/theory.json`)

**Ecuaciones** (`theory/modes.py`):
- modos w hidrostáticos con tapa rígida, W'' + (N²/c²)W = 0 con W(0) = W(−H) = 0 (Gill 1982 §6.11 [verificar]);
- R_n = c_n/f (m);
- Bu = R₁/L y R₁/W (adimensional);
- T_K = P/c₁ (s), con P el perímetro;
- Merian: T = 4L/((2n−1)√(gH)) para una cuenca cerrada en un extremo y abierta en el otro, y T = 2L/(n√(gH)) para una cuenca cerrada en ambos (Gill 1982 cap. 5 [verificar]);
- Poincaré transversal, ω = (f² + c_n²(π/W)²)^½ (Gill 1982 §10.4 [verificar]);
- f = 2Ω sin φ.

N² se calcula con la EOS lineal en las interfaces de la malla del modelo (`modes_column`).

### 2.1 Geometría

**Fiordo interior** (malla de la Task 1b, 200 m; `bathy.inner_fjord_mask`):
- área 526.8 km²;
- profundidad máxima 119.0 m y media 62.6 m;
- longitud a lo largo del eje L = 62.8 km. Es el camino más largo por agua (8-conexo, Dijkstra) desde las celdas junto a los estrechos;
- ancho medio W = área/L = 8.39 km;
- eje principal a 12.4° al norte del este (ENE–OSO);
- **perímetro P = 176.7 km**: lazo exterior del contorno 0.5 de la máscara (marching squares, `contourpy`), incluidas las bocas. La parte de costa es 170.0 km. El perímetro en escalera (aristas de celda) da 206.2 km (+17 %) y se descartó porque exagera las costas diagonales.

**Bahía de la tesis** (`bahia_01_expand_bat.bin` y `_dx/_dy.bin`, solo lectura): rectángulo de L = 121.0 km por W = 60.6 km y H = 164 m. Costa de tres lados 2L + W = 302.6 km; lazo cerrado 363.2 km.

### 2.2 Tabla

| Caso | f (s⁻¹) / T_i (h) | H modos (m) | c₁, c₂, c₃ (m s⁻¹) | R₁ (km) | Bu = R₁/L; R₁/W | T_K = P/c₁ (h) | ω_Poincaré,1/f | Merian (h): cerrada–abierta; cerrada–cerrada |
|---|---|---|---|---|---|---|---|---|
| VMF verano | 1.4266e-4 / 12.23 | 119.0 | 0.574, 0.335, 0.213 | 4.03 | 0.064; 0.48 | 85.5 | 1.81 | 2.82; 1.41 (H̄ = 62.6 m, L = 62.8 km) |
| VMF verano, H medio | ″ | 62.6 | 0.374, 0.192, 0.127 | 2.62 | 0.042; 0.31 | 131.2 | 1.40 | ″ |
| VMF invierno | ″ | 119.0 | 0.185, 0.075, 0.042 | 1.29 | 0.021; 0.15 | 265.9 | 1.11 | ″ |
| VMF invierno, H medio | ″ | 62.6 | 0.063, 0.031, 0.021 | 0.44 | 0.007; 0.05 | 777 | 1.01 | ″ |
| Bahía de la tesis, f de la tesis, perfil "realistic" | 6.97e-5 / 25.04 | 164 | 0.480, 0.215, 0.141 | 6.89 | 0.057; 0.11 | 174.9 (costa) / 210.0 (lazo) | 1.06 | 3.35; 1.68 |
| Bahía de la tesis a 78°N, perfil de verano VMF (R2) | 1.4266e-4 / 12.23 | 164 | 0.698, 0.358, 0.225 | 4.89 | 0.040; 0.081 | 120.5 / 144.6 | 1.03 | 3.35; 1.68 |

**Sensibilidades de c₁** (H = 119 m):
- verano: TEOS-10 0.520; "un voto por año" (mediana de las medianas anuales) 0.589;
- invierno: TEOS-10 0.186; "un voto por año" **0.063**.

El c₁ de invierno es incierto por un factor de ~3 (0.06–0.19 m s⁻¹): depende del escalón de salinidad cerca del fondo, que dan pocos perfiles.

### 2.3 Lo que estos números predicen (antes de correr)

1. **Resolución:** R₁/Δx = 20 (verano, H máx.) y 13 (verano, H medio) a 200 m, pero 6.5 y 2.2 en invierno. **En invierno el modo 1 está mal resuelto a 200 m** (y peor a 400 m). V4 solo prueba el verano: es una limitación de F2 y F3 que hay que declarar.
2. **Kelvin en 240 h:** en verano caben 2.8 vueltas (H máx.) o 1.8 (H medio). La spec pide "al menos dos", así que con H medio **no se cumple**. En invierno cabe menos de una vuelta (0.9 o 0.3): la onda de Kelvin interna de invierno no es resoluble espectralmente en 240 h.
3. **Poincaré y la banda inercial:** el fiordo es estrecho frente a πR₁ (W = 8.4 km). El primer modo de Poincaré transversal del modo vertical 1 está en 1.4–1.8 f (banda superinercial). Solo los modos verticales ≥ 3 (1.05–1.15 f) y la respuesta inercial forzada caen en la banda cercana a la inercial. En la bahía de la tesis (W = 60.6 km) el primer Poincaré está en 1.03–1.06 f, dentro de la banda cercana a la inercial.

   **Predicción del puente:** el paso de la bahía al fiordo real saca energía del modo 1 de la banda cercana a la inercial. Esto favorece la alternativa de H4a (Kelvin) por la geometría, no por f.
4. **Seiche de Merian** del fiordo interior: 1.4–2.8 h, muy por encima de f; no se mezcla con las bandas internas.

---

## 3. Bandas espectrales (propuesta; `spectral_bands` en `theory_numbers.py`)

Las bandas se definen relativas a f y a la seiche de Merian más lenta del fiordo interior (cuenca cerrada en un extremo y abierta en el otro, T_M = 2.82 h). Son contiguas y no se solapan:

| Banda | Límites | Fiordo a 78°N (cph; periodos) | Bins con 240 h (Δf = 1/240 cph) | Con 180 h de decaimiento libre |
|---|---|---|---|---|
| baja (no resuelta) | [0, 2Δf) | 0–0.0083; > 120 h | 2 | 2 |
| **Kelvin / subinercial** | [2Δf, 0.85 f) | 0.0083–0.0695; 14.4–120 h | 14.7 | 10.5 |
| **cercana a la inercial** | [0.85 f, 1.20 f) | 0.0695–0.0981; 10.2–14.4 h | 6.9 | 5.1 |
| superinercial | [1.20 f, 0.5/T_M) | 0.0981–0.1775; 5.6–10.2 h | 19.1 | 14.3 |
| **seiche** | [0.5/T_M, Nyquist] | 0.1775–1.5; 0.67–5.6 h | 317 | 238 |

**Justificación de cada límite:**
- **2Δf:** con una ventana Hann, la frecuencia más baja separable del cero.
- **[0.85, 1.20] f:** deja ≥ 5 bins en la banda cercana a la inercial incluso con 180 h, y deja un margen de ≥ 2Δf a cada lado de f para la fuga espectral. Es asimétrica hacia arriba porque las ondas cercanas a la inercial se desplazan al azul (Poincaré, vorticidad; Alford et al. 2016 [verificar]).
- **0.5/T_M:** la mitad de la frecuencia de Merian del fiordo interior. Cubre una cuenca efectiva el doble de larga (fiordo más Bellsund hasta la plataforma) y la incertidumbre de la profundidad.
- **Nyquist:** 1.5 cph, con registros cada 1200 s.

**Advertencias:**
- Los armónicos n ≥ 6 del circuito de Kelvin de verano (14.2, 12.2, 10.7 h) caen en la banda cercana a la inercial. Las bandas clasifican frecuencias, no mecanismos.
- Como diagnóstico **secundario y descriptivo**, prerregistrado aquí, se reporta la fracción de energía de cada banda que está a menos de 2R₁ de la costa (atrapada).
- Con H medio, el modo 1 del circuito de Kelvin de verano (131 h) cae en la banda "baja".

**Bahía de la tesis:**
- con f de la tesis y 240 h (R1k): baja < 0.0083; Kelvin 0.0083–0.0339; cercana a la inercial 0.0339–0.0479 (3.4 bins); superinercial 0.0479–0.149; seiche 0.149–2.0 cph (registros de 900 s);
- con las 120 h de R1, la banda cercana a la inercial tiene 1.7 bins y **no es resoluble**. Por eso H4d (R1k frente a R1) no usa bandas (§5);
- a 78°N (R2/R3): las mismas bandas que el fiordo, con el límite de seiche en 0.149 cph (Merian de la bahía, 3.35 h).

---

## 4. Validación V1–V4 (tolerancias propuestas)

Cada diferencia contra la referencia se clasifica como (a) error nuestro, (b) física ausente, (c) elección numérica o (d) error en la referencia. Un criterio que falla no se redefine.

- **V1, la copia reproduce la tesis (R1 frente a la corrida guardada):**
  - **V1-exacto:** máx |Δη| ≤ 1×10⁻¹⁰ m en todos los registros comunes. Esto es lo esperado con el mismo código, las mismas opciones y la misma precisión; solo cambia la E/S, de MNC a MDS.
  - **V1-espectral**, si no es exacto: (i) el pico de la seiche de la bahía en η, en la estación de la boca, cae en el mismo bin de frecuencia (Δf = 1/120 cph); y (ii) la DEP en el pico difiere ≤ 10 %.
  - **Justificación:** un cambio de compilador o de opciones produce diferencias de redondeo de crecimiento lento. Un error de configuración (f, perfil, forzamiento) cambia la frecuencia o la amplitud mucho más que un 10 %.
  - **De paso:** el f0 del `data` es 6.97×10⁻⁵ s⁻¹, con periodo inercial de 25.04 h. Los scripts de la tesis usan ~31 h. V1 lo resuelve con el espectro de R1.
  - **Nota para la Task 3:** la suma del `delR` de la tesis es 599.9 m, pero la plataforma de su batimetría es de 1000 m. MITgcm trunca la columna en el fondo de la malla, así que la profundidad efectiva fuera de la bahía es de 600 m. No afecta a V1 (misma configuración), pero hay que declararlo.
- **V2, frecuencias contra la teoría (F1):**
  - **V2a:** el pico de energía cinética horizontal dentro de la banda cercana a la inercial, **en la estación de la plataforma** (W ≫ R₁, donde la oscilación inercial no está restringida), está a ≤ 1Δf de f a 78°N. Δf es el de la ventana analizada.
    - Dentro del fiordo no se exige f: por §2.3, la teoría predice picos de 1.05–1.15 f (modos verticales altos) y de 1.4–1.8 f para el modo 1. Se reportan como prueba descriptiva de la teoría.
  - **V2b:** la velocidad de fase de la onda de Kelvin interna en F1 está dentro de [0.8 c₁(H̄), 1.2 c₁(H_máx)] = **[0.30, 0.69] m s⁻¹**.
    - Se mide por el tiempo de llegada del frente de desplazamiento de la picnoclina a ≥ 3 estaciones costeras sucesivas del fiordo interior, sobre la distancia a lo largo de la costa.
    - El rango cubre la incertidumbre de H del modelo plano, más ±20 % por la pendiente y la no linealidad. S&F14 obtuvieron 0.40–0.45 m s⁻¹ en su modelo, frente a 0.42 de la teoría de dos capas (p. 94), un acuerdo de ~7 %.
- **V3, M2 sola (contra S&F14):** pasa si se cumplen las tres:
  1. **presencia:** desplazamiento de la picnoclina a frecuencia M2 atrapado en la costa, con escala de decaimiento hacia fuera de la costa de 1–8 km (S&F14 ven un radio de 2–3 km, p. 93; R₁ = 2.6–4.0 km);
  2. **sentido ciclónico:** costa a la derecha, hacia el este por la costa sur del fiordo interior (S&F14 §4.2 y §4.3);
  3. **velocidad de fase** dentro de un factor 2 de 0.42 m s⁻¹, es decir 0.21–0.84 m s⁻¹. Es más estricto que el "orden de magnitud" de la spec.
- **V4, resolución (F1 a 400 m frente a 200 m):**
  - **V4a:** el flujo baroclínico neto hacia fuera, integrado en 0–240 h y sumado sobre los dos estrechos, difiere ≤ 25 %.
  - **V4b:** la fracción del flujo de cada banda resuelta (Kelvin, cercana a la inercial, superinercial, seiche) difiere ≤ 0.10 en valor absoluto.
  - El flujo de cada estrecho se reporta por separado y es descriptivo. A 400 m, Mariasundet (0.5–0.7 km) mide 1–2 celdas, así que se espera que falle solo.
  - **Justificación:** R₁/Δx baja de 20 a 10 en verano. El 25 % es la tolerancia más estricta que se puede pedir sin exigir convergencia a una resolución que la Task 1b ya mostró que tiene errores de ancho de ~1 celda en los estrechos.
  - Si falla, se declara la limitación (spec §6).

---

## 5. Hipótesis H4a–H4d (umbrales propuestos)

**Regla general:** un umbral de cociente debe superar la incertidumbre de resolución (V4a, 25 %). Por eso se usa un **factor 1.5** en H4c y H4d: un cociente entre 1/1.5 y 1.5 es "no distinguible". Todos los flujos son baroclínicos, por el método de Nash et al. 2005 [verificar] como en la tesis, y se integran en 0–240 h (0–120 h para R1k frente a R1). Las corridas que se comparan usan el mismo pulso, la misma duración y el mismo Δt.

- **H4a, partición (F1):** φ_b es la fracción del flujo baroclínico neto hacia fuera por los dos estrechos que lleva la banda b, sobre la suma de las cuatro bandas resueltas. La banda baja se reporta aparte.
  - **Apoyada** (como la tesis) si φ_NI > 0.5.
  - **Alternativa** (domina Kelvin) si φ_K > 0.5.
  - Si no, es **mixta**.
  - Si la banda baja lleva más del 25 % del flujo total en valor absoluto, H4a se marca "no concluyente con 240 h" (§2.3, punto 2).
  - La spec dice que se acepta el resultado que salga; la predicción de §2.3, punto 3, favorece la alternativa.
- **H4b, hielo (F2 frente a F3):**
  - **(i)** Ē_in(F2)/Ē_in(F3) ≤ 0.5, con Ē_in la energía baroclínica (cinética más potencial disponible) del fiordo interior promediada en 0–240 h, sin la banda baja.
    - **Justificación:** bajo hielo fijo, el viento no trabaja dentro del fiordo. El factor 2 es el menor que excede la incertidumbre de resolución y de residuo del balance.
  - **(ii)** "entra por los estrechos": F_bc,in ≥ |F_bt,in|, con los flujos netos hacia dentro por los dos estrechos, integrados en 0–240 h.
    - **Por qué este criterio:** en F2 el forzamiento local dentro del fiordo es nulo por construcción. Las fuentes de energía interna adentro son el flujo baroclínico que entra y la conversión, dentro del fiordo, del flujo barotrópico que entra. Si F_bc,in ≥ |F_bt,in|, más de la mitad llega ya como onda interna, aunque todo el barotrópico se convirtiera.
  - Se reportan (i) y (ii) por separado.
- **H4c, estación (F1 frente a F3, agua abierta):**
  - **(i)** W_NI(F1)/W_NI(F3) ≥ 1.5, con W_NI = ∫∫ τ·u_NI dA dt el trabajo del viento sobre la velocidad superficial filtrada en la banda cercana a la inercial, en el fiordo interior;
  - **(ii)** el cociente del flujo baroclínico radiado hacia fuera por los estrechos (todas las bandas resueltas), F1/F3, es ≥ 1.5.
  - **Apoyada** si se cumplen las dos; **refutada** si las dos son ≤ 1/1.5; **no distinguible** en otro caso.
  - **Advertencia:** el modo 1 de invierno está mal resuelto (§2.3, punto 1). Si H4c sale apoyada, puede deberse en parte a la resolución, y se declara.
- **H4d, cierre:** η = F_bc,out/W_viento, la fracción del trabajo del viento que sale como flujo baroclínico (por la boca de la bahía en R1/R1k y por los estrechos en F1/F1c).
  - **Apoyada** en cada par por separado si η(KPP)/η(cierre de la tesis) ≥ 1.5; **refutada** si ≤ 1/1.5; **no distinguible** en otro caso.
  - Pares: R1k frente a R1 (0–120 h) y F1 frente a F1c (0–240 h). No usa bandas (§3).
  - La tesis radiaba 0.3 % y 0.08 % (bitácora del 11 de noviembre); el factor 1.5 es deliberadamente modesto.
- **Puente:** descriptivo, sin umbral (spec §7).

---

## 6. Parámetros de mezcla y numéricos (propuesta)

| Parámetro | Valor | Fuente |
|---|---|---|
| Viscosidad vertical de fondo `viscAr` (KPP la usa como fondo: `viscArNr(1)` en `pkg/kpp/kpp_routines.F:1197`) | 1.0×10⁻⁴ m² s⁻¹ | L94 ec. 29 (p. 373), ν_w = 1.0×10⁻⁴ m² s⁻¹ (onda interna) |
| Difusividad vertical de fondo `diffKrT = diffKrS` (fondo de KPP: `CALC_3D_DIFFUSIVITY`, `pkg/kpp/kpp_calc.F:515-526`; suma en `kpp_routines.F:1198-1199`) | 1.0×10⁻⁵ m² s⁻¹ | L94 ec. 29, κ_w = 0.1×10⁻⁴ m² s⁻¹; es el mismo valor de difusividad que la tesis (`data.stage1`), así que F1 y F1c solo difieren en la capa límite de KPP, la mezcla por cizalla y la viscosidad de fondo |
| Ri crítico de la capa límite `Ricr` | 0.3 (valor por omisión) | `pkg/kpp/kpp_readparms.F:112`; L94 §5, p. 377 ("Ric = 0.3") |
| Mezcla interior por cizalla: `Riinfty`, `difm0 = difs0 = dift0` | 0.7; 5×10⁻³ m² s⁻¹ (valores por omisión) | `kpp_readparms.F:129, 132-134`; L94 ec. 28 (p. 373): ν₀ = 50×10⁻⁴ m² s⁻¹, Ri₀ = 0.7, p₁ = 3 |
| Convección `difmcon = difscon = diftcon`, `BVSQcon` | 0.1 m² s⁻¹; −0.2×10⁻⁴ s⁻² (valores por omisión) | `kpp_readparms.F:130, 136-138` |
| Resto de KPP | valores por omisión de MITgcm; `KPP_OPTIONS.h` sin cambios (con `KPP_GHAT`, sin `EXCLUDE_KPP_SHEAR_MIX`) | `pkg/kpp/KPP_OPTIONS.h:62, 65` |
| Viscosidad horizontal | **Smagorinsky** `viscC2smag = 2.2` y biarmónica `viscC4smag = 1` | **S&F14 §4.1 (p. 92): la misma configuración de MITgcm en Van Mijenfjorden** ("Smagorinsky type with value 2.2 along with a small biharmonic viscosity factor (viscC4smag = 1) as suggested by Griffies and Hallberg (2000)"). La documentación de MITgcm da 2.2–4 para el océano (`doc/algorithm/algorithm.rst:2449-2451`) |
| Difusividad horizontal de T y S | 0 (la difusión numérica la pone el limitador) | propuesta; F1c usa la misma, para aislar el cierre vertical |
| Arrastre de fondo `bottomDragQuadratic` | 2.5×10⁻³ | S&F14 §4.1, p. 92 ("quadratic bottom drag coefficient of 0.0025"); MITgcm: rango típico 0.001–0.003 (`doc/algorithm/algorithm.rst:1636-1638`); el valor por omisión es 0 (`model/src/set_defaults.F:137`) |
| Advección de T y S | `tempAdvScheme = saltAdvScheme = 77` (limitador de flujo de 2.º orden, Superbee) | `doc/algorithm/adv-schemes.rst:237-239, 455`. Se declara la mezcla numérica (Klingbeil et al. 2014, en la biblioteca de la tesis). El esquema 33 (DST3 con limitador) es menos compresivo; queda como alternativa si el PD1 muestra escalones en la picnoclina |

**Leith o Smagorinsky.** Se propone Smagorinsky y no Leith. Leith responde solo a la vorticidad ("can only respond to buildup of vorticity at the grid scale", `doc/algorithm/algorithm.rst:2491, 2517-2518`); por eso la documentación añade `viscC2LeithD` para los modos divergentes (`algorithm.rst:2523-2545`), y las ondas internas son divergentes. Además, Smagorinsky 2.2 con `viscC4smag = 1` es la configuración con la que S&F14 reprodujeron las ondas de Kelvin internas en este mismo fiordo con MITgcm, lo que hace que V3 sea una comparación más limpia. La estabilidad se verifica en el cronometraje del PD1.

**Nota de S&F14 (descriptiva):** usaron Az = 1×10⁻³ m² s⁻¹, Δx = 100 m y 32 niveles (§4.1, p. 92). No se copia su Az, porque KPP lo calcula.

---

## 7. Viento (propuesta)

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

## 8. Hielo: coeficiente de arrastre hielo–agua (propuesta)

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

## 9. Para que Esteban decida

1. **Perfiles:** polígono, años 1990–2019, tolerancia de inversión de 0.004 kg m⁻³ (elegida tras ver los datos; §1.1) y mínimos de 20 perfiles y 5 por nivel.
2. **EOS:** mantener la de la tesis (sesgo de +10 % en c₁ de verano) o ajustarla a Van Mijenfjorden (§1.3).
3. **H para la teoría:** las bandas y V2b usan el rango [H medio, H máximo]. ¿Basta un solo H?
4. **Duración:** con H medio caben 1.8 vueltas de Kelvin en 240 h (< 2, spec §3.4). ¿Se alarga a ~280 h si el cronometraje del PD1 lo permite?
5. **Invierno:** R₁/Δx = 2–6.5 a 200 m. ¿Se acepta como limitación declarada de F2 y F3?
6. **Tolerancias y umbrales:** V1 (10⁻¹⁰ m y 10 %), V2a (1Δf en la plataforma), V2b ([0.30, 0.69] m s⁻¹), V3 (factor 2), V4 (25 % y 0.10), H4a (0.5), H4b (0.5 y F_bc ≥ |F_bt|), H4c y H4d (factor 1.5).
7. **Mezcla:** L94 de fondo (1×10⁻⁴ y 1×10⁻⁵), KPP por omisión, Smagorinsky 2.2 con biarmónico 1 (como S&F14), arrastre 2.5×10⁻³, esquema 77.
8. **Viento:** del ENE (078°) hacia la boca; F5 del NNO.
9. **Hielo:** C_w = 5.5×10⁻³, sin sensibilidad.
