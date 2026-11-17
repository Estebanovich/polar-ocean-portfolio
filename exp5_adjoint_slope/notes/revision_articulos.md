# Revisión del diseño del Exp. 5 contra los artículos descargados

Fecha: 2026-11-16. Revisión de solo lectura: no se editó ningún archivo del repositorio.

**Fuentes leídas (texto completo, `portfolio/literatura/*.md`):**
- Stewart & Thompson (S&T): el archivo se llama "2014", pero el artículo es **GRL 42, 432–440, publicado en 2015**.
- Dettling et al. (2024), JPO 54, 1675–1690.
- Si, Stewart & Eisenman (2023), Sci. Adv. 9, eadd7049.
- Thompson et al. (2018), Rev. Geophys.
- Rodriguez, Mazloff & Gille (2016), JGR Oceans 121, 3337–3349.
- Mazloff (2012), J. Climate 25, 2279–2290.
- Mazloff, Heimbach & Wunsch (2010), JPO 40.
- Losch & Heimbach (2007), JPO 37.
- Heimbach & Losch (2012), Ann. Glaciol. El texto convertido no tiene las pp. 62–64, donde está la configuración (§2.3).
- Lea, Allen & Haine (2000), Tellus 52A, 523–532.
- McLay (2006), tesis, Southampton.
- Marotzke et al. (1999), JGR 104(C12).

**Convenciones:**
- "p." es la página de la revista.
- **[no verificado]** marca lo que no pude comprobar en el texto local.
- **[inferencia]** marca un juicio mío, no del artículo.
- **[cálculo propio]** marca una cuenta mía con los números del artículo o de nuestro código.

---

## 0. Resumen de hallazgos (por impacto)

1. **Nuestro perfil de viento hace casi lo contrario de lo que, según la literatura, sostiene el frente.**
   - Thompson et al. (2018, §2.5, p. 750) dicen que los vientos del este tienen "a tendency for the wind stress curl to induce downwelling velocities (Ekman pumping) over the continental slope".
   - S&T (2015, p. 434) y Dettling et al. (2024, p. 1677) usan solo viento del este, con el máximo 25 km mar afuera del centro del talud.
   - El nuestro es τ = −0.05 constante en y < 150 km, así que no hay rotacional ni bombeo de Ekman sobre la plataforma ni sobre el talud. La convergencia de Ekman cae toda en la pared costera.
   - Entre 150 y 250 km la rampa hacia +0.1 diverge: hay surgencia de ≈ 0.9 m/día, a 50–150 km del talud [cálculo propio].
   - [inferencia] Esto sube el CDW junto al talud en vez de hundir la picnoclina en el quiebre.
2. **Nuestro GM/Redi es 5–30 veces más difusivo que el de la literatura, y además reducimos Redi, cosa que Dettling no hace.**
   - Dettling usa κ_Redi = 15 m² s⁻¹ uniforme (p. 1682). Si et al. usan 100 m² s⁻¹ para GM y Redi (Métodos, p. 9).
   - Nosotros usamos 500 m² s⁻¹ en los dos. Al pie del talud bajan a 50, en la línea de J quedan en ≈ 117 y sobre la plataforma siguen en 500.
   - Dettling (p. 1684): con una difusividad alta, "WDW can directly access the continental shelf and erode the V-shaped isopycnal structure".
3. **Ningún artículo de proceso consigue una plataforma fría solo con relajación superficial.**
   - S&T y Dettling inyectan sal en la polinia sobre los primeros 50 km y usan un modelo termodinámico de hielo marino. Su régimen es de plataforma **densa**.
   - Si et al. relajan **toda la columna** en una esponja sur de 20 km: θ al punto de congelación y S fijada (Métodos, p. 8).
   - Nuestra plataforma solo se relaja en los 10 m de arriba, y el objetivo de S de la plataforma (34.3) es más salado que el de mar abierto (34.2). Eso no es una plataforma fresca en el sentido de Si et al. ni de Thompson et al.
4. **Plataforma fresca no quiere decir J ≈ 0 ni J pequeño: el criterio V4 está mal planteado.**
   - Si et al. (p. 3): el flujo de calor hacia la costa es casi nulo **en el caso de referencia, sin gradiente de salinidad**, y "both freshening and salinification of the shelf waters relative to the offshore waters lead to increased heat flux onto the continental shelf".
   - Sin marea, en la plataforma fresca el flujo de remolinos compensa por completo al medio, "leading to near-zero total onshore heat transport" (p. 6, fig. S9A). Nuestro modelo no tiene marea.
5. **El GM estándar falla justo en el régimen de plataforma fresca (Si et al. 2023), y GM sensible al talud solo está calibrado en plataforma densa (Dettling).**
   - Si et al. (p. 5): a 5–10 km con GM-Redi, "the standard GM-Redi eddy parameterization scheme is not enough to reproduce the positive feedback, leading to even worse performance".
   - Esto pega directo en H5b y en la validez física del adjunto grueso.
6. **El mecanismo de Rodriguez et al. (2016) no existe en nuestro canal, así que H5a está mal motivada.**
   - Rodriguez et al. cierran el balance de vorticidad con el término βv y la estructura a lo largo de la plataforma (Ec. 3, p. 3345).
   - Nuestro canal está en un plano f (β = 0) y su batimetría es uniforme en x.
   - El "rotacional manda donde la topografía guía la corriente" de Mazloff (2012) se refiere a contornos f/H obstruidos (p. 2289). Aquí no lo están.
7. **Errores de cita.**
   - "Mazloff (2011)" es **Mazloff (2012)**, J. Climate 25, 2279–2290, doi:10.1175/JCLI-D-11-00030.1.
   - "Stewart & Thompson (2014)" y "(2015)" son **el mismo artículo** (GRL 2015).
   - Las σ del PD0 no son las de Mazloff: él usa 0.1 N m⁻² (Fig. 3, p. 2283).
8. **El spec (§4) se contradice:** dice que "365 días … menos que el tiempo no lineal de McLay, unos 200 días". 365 es más que 200.
9. **La resolución del fino (2.5 km) está casi en la zona donde el mecanismo desaparece.**
   - S&T (p. 439): el transporte de CDW por remolinos "is completely absent even at … 5 km and is underestimated by close to 50% even at 2 km".
   - Si et al.: el radio de Rossby es ≈ 2.3 km sobre la plataforma.
10. **El diseño de dos modelos es justo el que propone McLay (2006, cap. 5, pp. 84–85).**
    - "a multi grid approach … a non eddy resolving ocean model and its adjoint are used to identify regions where the sensitivity was likely to be high in an eddy resolving ocean model. This could be tested using perturbed forward model experiments in both…"
    - Conviene citarlo de forma literal en el README.

---

## 1. Configuración: nosotros contra S&T (2015), Dettling et al. (2024) y Si et al. (2023)

### 1.1 Tabla comparativa

| Elemento | S&T 2015 (GRL) | Dettling 2024 (JPO) | Si 2023 (Sci. Adv.) | Exp. 5 (`chan10`/`chan2p5`) |
|---|---|---|---|---|
| Dominio (a través × a lo largo × prof.) | 450 × 400 km × 3000 m, reentrante (p. 434) | 400 × 450 km, periódico en x (p. 1677) | 450 × 400 km × 4 km (Métodos, p. 8) | 450 × 400 km × 3000 m ✓ |
| f₀ | −1.31×10⁻⁴, β = 0 (p. 434) | −1.31×10⁻⁴, β = 0 (Tabla 1) | [no verificado] | −1.4×10⁻⁴, β = 0 (≈) |
| Plataforma | 500 m (Tabla 1, rango 300–700) | 500 m; ancho de plataforma 125 km (Tabla 1) | 500 m (Métodos: "z = −500 m"), con fosas | 500 m, 0–75 km |
| Talud | centro y = 200 km, W = 150 km (p. 434; rango 50–250) | centro 200 km, W = 150 km; pendiente media ≈ 0.02 (Tabla 1; p. 1676) | variado ("WS"); valor [no verificado] | **centro 100 km, tanh w = 12.5 km: pendiente máx. 0.1 [cálculo propio]**. Fuera de rango |
| Puntos de malla gruesa sobre el talud | — | "15 grid points" a 10 km (p. 1677) | — | **≈ 3 a 10 km** [cálculo propio] |
| Viento | solo del este, sinusoidal por tramos; 0 en y = 50 y 450 km; τ_max 0.075 a 25 km mar afuera del centro del talud (p. 434; rango 0–0.1) | τ_max −0.075 en Y_w = 225 km (Tabla 1) | del este y meridional, máximos en la costa y decrecen linealmente a 0 en el norte (Métodos, p. 8) | −0.05 en y < 150 km; **+0.10 (del oeste)** en y > 250 km. **Los vientos del oeste no aparecen en ninguno de los tres** |
| Superficie | 1 m de hielo fijo al 100 %, flujos de calor y sal de dos ecuaciones (Schmidt 2004); sin forzamiento térmico en la polinia (p. 434) | igual: hielo termodinámico de dos ecuaciones, γ = 6×10⁻⁵ m s⁻¹ (p. 1677) | hielo dinámico y termodinámico, atmósfera fija que minimiza flujos netos (Métodos, p. 8) | **relajación de T y S en 30 días hacia el nivel 1 inicial**. Sin hielo |
| Forzamiento de flotabilidad en la plataforma | sal de polinia 2.5×10⁻³ g m⁻² s⁻¹ en los primeros 50 km (p. 434). "We found it necessary to impose this additional, explicit sink of buoyancy on the continental shelf in order to achieve a statistically steady state" | igual, 2.5 mg m⁻² s⁻¹ en 50 km (Tabla 1) | **esponja sur de 20 km: θ = congelación en toda la columna; S uniforme en la vertical (plataforma fresca) o creciente con la profundidad (densa)**; relajación de 10 a 0.5 días (p. 8) | **ninguno** salvo la relajación superficial; objetivo de SSS en la plataforma 34.3 > mar abierto 34.2 |
| Esponja norte | 50 km, T/S hacia un perfil de Weddell (Thompson & Heywood 2008) **sin AASW ni AABW**, de 0 a 56 días; **velocidad a 0 en 28 días** (p. 434) | 50 km; T/S en 54 días, **velocidad en 27 días** (p. 1677) | 20 km, Kapp Norvegia en invierno; contornos abiertos con marea (p. 8) | 30 km, T/S en 10 días hacia los perfiles iniciales; **sin relajar la velocidad** |
| Perfiles | Weddell observado, suavizado [valores no verificados] | Weddell (Thompson & Heywood 2008), Fig. 1c,d [valores no verificados] | Norte: Kapp Norvegia; superficie sur = superficie norte = **34.12 psu** en la referencia (p. 9) | mar afuera θ −1.8 → 0.8, S 34.2 → 34.7, termoclina tanh a 300 m; plataforma θ −1.85, S 34.3 → 34.5 |
| EOS | McDougall et al. (2003) | McDougall et al. (2003) | polinomial de alto orden (p. 8) | **lineal** (α 3×10⁻⁵, β 7.4×10⁻⁴) |
| Resolución | 1 km (referencia), rango 0.5–10 km; 53 niveles de 13–100 m; Δt 179 s (p. 435) | 1 km / 10 km; 77 niveles de 13–100 m; Δt 180/1200 s | 2 km (más 5 y 10 km); 70 niveles de 10–100 m | 10 km / 2.5 km; 35 niveles de 10–161 m |
| Mezcla y fricción | κ_z 5×10⁻⁶ + KPP; A_h 12; biarmónica y Leith; arrastre lineal 1×10⁻³ m s⁻¹ (p. 435) | igual (Tabla 1), en las dos resoluciones | Smagorinsky biarmónico (factor 4) | κ_z 1×10⁻⁵, sin KPP, `ivdc`; **viscAh 1000** (≈ 80 × A_h de Dettling a 10 km); arrastre cuadrático 2.5×10⁻³ |
| GM/Redi | ninguno (resuelve) | κ_GM diagnosticado o sensible al talud (Ec. 24, ε_c = 800); **κ_Redi = 15 m² s⁻¹ uniforme**; taper de Gerdes et al. 1991 (pp. 1678, 1682) | a 5–10 km: κ_GM = κ_Redi = **100** m² s⁻¹; pendiente máx. 0.025; DM95 con pendiente crítica 0.025 y ancho 0.0025 (p. 9) | κ_GM = κ_Redi = **500**, × (1 − 0.9 e^{−((y−100)/25)²}) **en los dos**; dm95, maxSlope 0.01 |
| Marea | no | no | sí; es esencial para la plataforma fresca | no |
| Spin-up y media | 40 años a 2 km + 5–20 años a 1 km; media de los últimos 5 años (p. 436) | 40 + 40 años a 10 km; los últimos 5 años, mensual (pp. 1677–1678) | 20 años + 7 de análisis (p. 9) | G0 30 años; **el fino F0 solo 6–12 meses** |

### 1.2 Qué valores nuestros quedan fuera de su rango

- **Pendiente y ancho del talud.**
  - El ancho efectivo es ≈ 25 km contra 50–250 km en S&T (Tabla 1). La pendiente máxima es ≈ 0.1 contra ≈ 0.02 en Dettling (p. 1676). Es 5 veces más empinado que el talud más empinado que probó S&T [cálculo propio].
  - S&T (p. 437): "Increasing W_slope reduces the topographic potential vorticity gradient, allowing the slope to support larger eddies that transport CDW shoreward more rapidly." Thompson et al. (2018, §2.5): "Increasing slope steepness increases the topographic PV gradient and stabilizes the along-slope flow."
  - Ojo: Si et al. encuentran que un talud **más empinado** aumenta el flujo hacia la costa en los dos regímenes (p. 5, fig. S2). El efecto del ancho depende del régimen.
  - A 10 km el talud cabe en unas 3 celdas. Dettling se aseguró de tener 15 (p. 1677).
- **Vientos del oeste mar afuera (+0.1 N m⁻²).** No hay nada así en los tres modelos. S&T sitúa "the continental shelf … to the left of the winds" en todo el dominio (p. 434).
- **Magnitud del viento del este (−0.05).** Está dentro del rango de S&T (0–0.1), pero por debajo de su referencia (0.075).
- **Viscosidad.** viscAh = 1000 m² s⁻¹ contra A_h = 12 m² s⁻¹ más biarmónica en el modelo de 10 km de Dettling, que tiene el mismo Δt de 1200 s.
- **Difusividades de GM y Redi** (ver §0.2).
- **Esponja.** Es más corta (30 contra 50 km) y más rápida (10 contra 54–56 días), y no relaja la velocidad. S&T relajan la velocidad "to avoid the development of a spurious jet at the northern wall" (p. 434).
- **Relajación superficial y ausencia de un sumidero de flotabilidad o de calor en la plataforma.** No hay nada análogo en ninguno de los tres.

### 1.3 ¿Qué mantiene el CDW fuera de la plataforma en sus corridas?

- **Viento (S&T 2015).**
  - "F_CDW decreases almost linearly with τ_max" (p. 436).
  - "Increasing τ_max steepens the isopycnals in the ASF, deepening the pycnocline at the shelf break until CDW can no longer mix across it" (p. 437).
  - "Shifting the surface wind stress maximum simply changes the effective wind stress over the continental slope, which in turn modifies F_CDW by changing the pycnocline depth at the shelf break" (p. 437).
  - Lo que cuenta es el esfuerzo **sobre el talud** y la convergencia de Ekman que lo acompaña.
- **Plataforma fresca (Thompson et al. 2018, §2.1, p. 748).**
  - "persistent and strong easterly winds … onshore Ekman transport … Cold shelf waters, weak cross-slope exchange due to a strong ASC, and a converging surface Ekman transport that produces downwelling velocities, all combine to establish a strong front."
  - El régimen se define así: "Fresh Shelf regions where neutral density layers γⁿ > 28.0 kg/m³ do not extend on to the shelf" (p. 747).
  - El agua dulce aparece como estado ("cold and fresh shelf water"), no como el mecanismo necesario. Además, "The influence of surface buoyancy fluxes on the structure of the ASF is less clear" (§2.5, p. 750).
- **Si et al. (2023): el frío de la plataforma lo impone la frontera.**
  - El calor que converge en la plataforma "is locally balanced by ocean-sea ice heat flux due to sea ice melting" (p. 4).
  - La esponja sur fija θ en el punto de congelación en toda la columna.
  - En la plataforma fresca, el agua dulce **aumenta** el transporte hacia la costa: hay una celda de vuelco somera conducida por remolinos que exporta agua superficial dulce y trae agua cálida por debajo (p. 5, Fig. 5A).
  - Que la plataforma quede fría no significa que el CDW esté aislado. Significa que hay un sumidero de calor fuerte.
- **Plataforma densa (S&T; Dettling).** El CDW **sí** cruza, por remolinos a lo largo de isopicnas que conectan con la plataforma. La sal de polinia es la fuente de energía potencial disponible (S&T, p. 437). Ninguno de los dos es un caso de plataforma fresca.

**¿Sus montajes dan una plataforma fresca con nuestros ingredientes?**
- No. Ninguno tiene nuestra combinación: relajación superficial sola, vientos del oeste, talud de 25 km y GM/Redi de 500.
- [inferencia, cálculo propio con nuestra EOS lineal] El agua del fondo de la plataforma inicial (−1.85 °C, 34.5) tiene ρ′ ≈ +2.6×10⁻⁵ (relativa a T = −1, S = 34.5). Esa es la densidad del perfil de mar abierto a ≈ 360 m, donde θ ≈ +0.2 °C.
- Así que **una isopicna de agua de termoclina cálida (≈ 0.2 °C) conecta con el fondo de la plataforma**. Por la definición de Thompson et al., eso no es una plataforma fresca estricta.
- Con Redi = 117–500 m² s⁻¹ a lo largo de esa isopicna y sin sumidero de calor en profundidad, lo esperable es que la plataforma se caliente hacia la temperatura de esa isopicna.
- El nuevo S0 sin el error de periodicidad en y dirá cuánto pesa este efecto.

### 1.4 Nuestro objetivo de salinidad superficial (34.3 en la plataforma contra 34.2 mar afuera)

- **Si et al.**
  - La referencia ("no hay gradiente") iguala la salinidad superficial sur y norte en 34.12 psu.
  - La plataforma fresca usa S uniforme en la vertical y **más dulce** que la referencia.
  - La densa usa 34.12 en superficie y S creciente con la profundidad (p. 9).
  - La palanca es "the shelf bottom salinity (or, equivalently, the maximum meridional salinity gradient at mid-depth)" (p. 9).
- **En nuestro caso.**
  - La SSS de la plataforma es más salada que la de mar abierto, y la S del fondo (34.5) solo es más dulce que el CDW (34.7). Por la clasificación de Si et al., estamos entre la referencia y una plataforma **densa** débil, no fresca.
  - Además, solo relajamos los 10 m de arriba. El fondo de la plataforma no está fijado por nada, salvo la convección y la mezcla.
- **[inferencia]** La variante SF (33.8 en superficie) va en la dirección correcta, pero solo actúa en la superficie. El análogo de Si et al. sería relajar toda la columna de la plataforma: una esponja sur con S uniforme y θ = θ_f.

---

## 2. La función de costo J

| Fuente | Qué mide | Referencia de θ | Dónde | Media/remolino | Promedio temporal |
|---|---|---|---|---|---|
| S&T 2015 | F_CDW: volumen por líneas de corriente del vuelco TEM en densidad potencial referida localmente, con S > 34.543, que conectan el quiebre (y = 100 km) con la esponja (p. 436). El calor se estima a posteriori | "relative to the surface freezing temperature of ∼ −1.8 °C"; con CDW a 0 °C, 3.5 Sv circumpolares ≈ 25 TW (p. 438) | quiebre, y = 100 km; isóbata de 1000 m en la comparación con observaciones (Fig. 3i) | la descomposición se infiere: "the shoreward CDW transport is entirely due to mesoscale eddy transfer" | últimos 5 años |
| Dettling 2024 | F_θ,mean = −c_p ρ₀ ∫∫ v̄ θ̄ dz dx; F_θ,eddy = −c_p ρ₀ ∫∫ v′θ′ (Ecs. 6–7); F_GM = −c_p ρ₀ ∫∫ v*θ; F_Redi = −c_p ρ₀ ∫∫ κ_Redi(∂θ/∂y + s ∂θ/∂z) (Ecs. 11–13) | **sin referencia**: θ absoluta en °C | función de y (Figs. 3, 6, 7); integrada al dominio sin la esponja (Fig. 9d) | la media es de 5 años más zonal. "the eddy heat flux is balanced by a mean offshore heat flux over the shelf and slope" (p. 1679). **Compara GM/Redi + remolinos resueltos contra el flujo de remolinos del fino, no el total** | 5 años (mensual) |
| Si 2023 | "Vertically and zonally integrated meridional advective heat flux" (Fig. 2A). F_shelf = media sobre la banda **y = 50–75 km**; F_slope = media sobre **y = 125–150 km** (p. 5). Descomposición en media, remolinos y marea (Ecs. 5–8), con medias diarias y de 7 años | "θ_ref = 0 °C" en la función de calor (Ec. 2): "The interpretation … is very similar for different choices of the reference temperature" (p. 9) | bandas en y, no una sola línea | además separa advección y difusión de remolinos en coordenadas σ₂ (Ecs. 13–14) | 7 años |
| Rodriguez 2016 | transporte advectivo de temperatura a través del contorno de **2600 m**, por cajas sobre y bajo 225 m; × C_v = 3985, ρ₀ = 1035 (pp. 3338, 3341) | θ absoluta; justificado porque "the net mass divergence in the modeled ASE is negligible" (p. 3341) | isóbata de 2600 m | u⊥θ = media geostrófica + remolino geostrófico + ageostrófico (Ec. 2, p. 3344) | 2005–2010, mensual y anual |
| Marotzke 1999 | Q = P⁻¹ c_p ρ₀ ∫∫∫ v T dx dz dt (Ec. 10) | T absoluta (cuenca cerrada) | 29°N | — | media anual |
| **Exp. 5** | −ρ₀c_p [Σ v (θ̄ − θ_f) + F_y^{GM/Redi}] en una cara v | θ_f = −1.9 °C | y = 90 km, isóbata de ~920 m (talud superior) | total: resuelto + GM/Redi | últimos 30 días de GA (360 d) o GB (180 d) |

**¿Es consistente?**
- **Referencia.** θ_f = −1.9 está en línea con S&T (−1.8). En el canal el transporte neto de volumen es casi nulo, así que la referencia casi no importa. Si et al. lo dicen explícitamente y usan 0 °C (p. 9). Bien.
- **Línea.**
  - 920 m se parece al estándar de 1000 m (S&T, Fig. 3; Thompson et al., Fig. 4). Pero cae **en el talud superior**, donde los remolinos están suprimidos.
  - Si et al. (p. 5): "Near the continental shelf break, the eddy heat transport is weaker … because the baroclinic eddies are suppressed there by strong topographic vorticity gradients."
  - Y (p. 5): "the cross-slope heat transport largely depends on the tidal and mean components; since eddies are suppressed over the slope, the cross-slope heat transport does not rely on the model's ability to resolve mesoscale eddies."
  - Sin marea y en la línea del talud, J depende del residuo entre la media y lo que quede de remolinos [inferencia].
- **Total contra remolinos.**
  - Nuestra J total es lo que entra al balance de calor de la plataforma, así que es una elección defendible.
  - En equilibrio, J ≈ el flujo de calor superficial al sur de la línea (la relajación) más dH/dt [inferencia].
  - Dettling y S&T reportan sobre todo la parte de remolinos o de CDW. Hay que reportar la descomposición (media, remolinos resueltos, GM, Redi) como en Dettling (Ecs. 6–7 y 11–13) y en Si et al. (Ecs. 5–8).
- **Redi dentro de J.** Dettling (p. 1683): "On the upper slope, the Redi scheme takes over and captures some of the shoreward heat flux across the shelf break." Con nuestro Redi ≈ 117 m² s⁻¹ en la línea (8 veces el de Dettling), es probable que J esté dominada por Redi [inferencia]. Hay que reportar la fracción de Redi en J.
- **Promedio de 30 días.** Coincide con Heimbach & Losch (2012, Ec. 1; p. 65: "J, is an average over the last 30 days of the integration"). Es mucho más corto que en todos los estudios de proceso (5–7 años) y que el de Mazloff (2012), que promedia toda la ventana de 2 años. Para H5b en el fino, 30 días son ruidosos (§4).
- **Rodriguez et al. (2016) no tiene adjunto.**
  - Es un diagnóstico de balances (calor y vorticidad) en SOSE, que es un modelo ajustado con el adjunto (p. 3338).
  - **No define una función objetivo adjunta.** Su "transporte" es el advectivo a través de la isóbata de 2600 m, sobre 2005–2010.
  - La pregunta "¿cómo define Rodriguez su J adjunta?" no tiene respuesta: no la hay.

**Órdenes de magnitud para V4** [cálculo propio]:
- S&T: F_CDW = 0.07 Sv × 1.8 K × ρc_p ≈ **0.5 TW por 400 km** (plataforma densa, τ = 0.075). F_AASW = 0.20 Sv ≈ Ekman.
- Nuestro Ekman: 0.05/(ρ₀|f|) × 400 km ≈ 0.14 Sv.
- Rodriguez: 5.85 TW bajo 225 m sobre ≈ 890 km de borde (θ absoluta, plataforma cálida).
- Si et al.: F_shelf ≈ 0 en la referencia y en la plataforma fresca sin marea.

---

## 3. Práctica adjunta

| | Ventana | Promedio de J | Escalado (σ) | Validación del gradiente | GM/Redi en el adjunto |
|---|---|---|---|---|---|
| Mazloff 2012 (SOSE 1/6°) | 2 años (2005–2006) | J = transporte medio por Drake sobre **toda** la ventana (Ec. 1); otro de 15 min para validar (pp. 2280–2281) | "0.1 N m⁻² for momentum flux, 10 W m⁻² for heat flux, and 10⁻⁸ PSU kg m⁻² s⁻¹ for salt flux" (Fig. 3, p. 2283) | perturbación directa (5 m s⁻¹, gaussiana de 1°) contra la convolución adjunta: "average difference … 8%, with a maximum difference of 13%"; "the error is less than 10%" (pp. 2281–2282). No reporta grdchk | SOSE usa GM-Redi (p. 2280); κ isopicnal 10 m² s⁻¹ (Mazloff 2010, Tabla 2). No dice nada de forma, taper ni estabilidad |
| Mazloff 2010 | 2 años; 26 iteraciones | misfit | pesos = desviaciones estándar | — | κ_iso 10 m² s⁻¹ (Tabla 2) |
| Losch & Heimbach 2007 (4°, 2°) | 1, 10, 100 y 500 años | último año de la integración (p. 385) | incertidumbres del 50 % para viento y flotabilidad, 100 m para topografía | "We performed various gradient checks … but will omit the technical details". Cuasi-linealidad con ΔH = 100 m: "<30%" a 100 años; "order 5% for 1-yr integrations" (Tabla 2, p. 385). Distinguen grdchk (ε pequeño, corto) de cuasi-linealidad (finito, largo) | "For technical reasons, the current model configuration does not include … Gent and McWilliams" (p. 390) |
| Heimbach & Losch 2012 | [no verificado: pp. 62–64 ausentes] | "averaged over a period Δτ preceding the final time" (Ec. 1); Δτ = 30 días (p. 65) | no escalan ("we are not considering responses", p. 65) | [no verificado] | [no verificado] |
| McLay 2006 | 278 y 690 días | media temporal de toda la integración | — | FD con τx 1×10⁻⁵ N m⁻² en 5 celdas: 2 % (0.48°), 12 % (0.1° amortiguado) y 31 % (0.1° sin amortiguar) a 690 días (§3.2.2, p. 42). TLM contra par no lineal para el tiempo no lineal | a 0.48° **sin** GM: "More sophisticated schemes for isopycnal mixing are not used as they tend to be highly nonlinear and were found to cause the adjoint model to become unstable after a few timesteps" (p. 23) |
| Marotzke 1999 | 1 año | media anual (Ec. 10) | — | "We have confirmed the linearity of our perturbation run by scaling the salinity perturbation down by an order of magnitude" | — |
| **Exp. 5** | GA 360 d, GB 180 d | últimos 30 días | τx 0.02 N m⁻² por mes; T₀ 0.1 °C; S₀ 0.02; κ_v 1×10⁻⁵ | V1 grdchk < 1 % en ≥ 10 puntos; V2 \|Δ⁺+Δ⁻\|/\|Δ⁺−Δ⁻\| < 0.1; V3 < 30 % | skew-flux en el adjunto; ver abajo |

**Comparación.**
- **V1.** El < 1 % es más estricto que todo lo local, pero corresponde a grdchk con ε pequeño (la distinción de Losch & Heimbach). La fuente es Gaikwad et al. (2024), que no está en `literatura/` [no verificado aquí].
- **V2.** Es una versión cuantitativa del test de escalado de Marotzke et al. Bien.
- **V3.**
  - El umbral del 30 % viene de Losch & Heimbach, pero para **100 años**. Para 1 año ellos obtienen ~5 %, y Mazloff < 10 % en un modelo con remolinos.
  - Un revisor dirá que el 30 % es laxo para 360 días en un modelo no caótico.
  - El umbral está congelado y no debe cambiarse. Hay que reportar el error obtenido y compararlo con el ~5 % de L&H y el ~8 % de Mazloff.
- **Ventanas.**
  - 360 días está entre el año de L&H (≈ 5 %) y los 690 días que McLay valida en el modelo no caótico de 0.48° (FD 2 %). Es razonable para un modelo grueso **si no es caótico** (§4).
  - Mazloff promedia J sobre **toda** la ventana. Nosotros medimos la sensibilidad de los últimos 30 días, lo que da sensibilidades más ruidosas y menos "climáticas".
- **H5c.**
  - Los "~100 días" de Mazloff no son una semivida de anticipación.
  - El resumen dice "time scale of influence … on the order of 100 days". En el cuerpo, la respuesta a un impulso "is negligible after a month", con e-folding de ≈ 3–4 días (p. 2282). Los ~100 días son el tiempo en que la sensibilidad de la media "reach a quasi-equilibrium state on the order of 100 days" (p. 2283).
  - Además es una respuesta **barotrópica** ("This primary response, however, is barotropic", p. 2287).
  - La comparación natural para un J baroclínico es Heimbach & Losch: "dominant timescales of 30–60 days" (resumen, p. 59).

**Forma de GM, GM_Small_Number, taper y estabilidad (nuevo contexto: skew-flux inestable en el directo, con θ de −7 a +2.4 °C en 2 años; la forma advectiva es estable en el directo; su adjunto solo es estable con GM_Small_Number = 1e-6).**
- **Ningún artículo local trata la forma skew contra la advectiva, `GM_Small_Number` ni la estabilidad del adjunto de GM.** Busqué "skew", "small number", "advective form" y "bolus velocity" en los 12 artículos: sin resultados pertinentes.
- **Lo más cercano:**
  - **McLay (2006, p. 23):** los esquemas isopicnos "were found to cause the adjoint model to become unstable after a few timesteps". Es el único precedente local de inestabilidad del adjunto con mezcla isopicna, y la evitó quitándola.
  - **Si et al. (2023, p. 9):** a 10 km, κ_GM = κ_Redi = 100 m² s⁻¹, "maximum isopycnal slope is 0.025", DM95 "with a DM95 critical slope of 0.025 and a DM95 tapering width of 0.0025". No hablan de inestabilidad. Es el conjunto de parámetros más parecido al nuestro que funcionó en el directo.
  - **Dettling (p. 1678):** "To avoid numerical instability in the presence of large isopycnal slopes, we use the tapering scheme of Gerdes et al. (1991). No major differences were observed when testing other tapering schemes."
    - Describen GM con la función de corriente de bolo y la divergencia advectiva F_GM = −∇·(φu*) (Ecs. 1–3).
    - No dicen qué forma de MITgcm usaron. Su código está en github.com/nicolasdettling/weddell_gm [no verificado].
  - Losch & Heimbach no tenían GM. Mazloff no lo discute.
- **[inferencia]**
  - El rango de θ de −7 a +2.4 °C con skew-flux en el directo es una sobrepasada no física, muy por debajo del punto de congelación. Es coherente con tres cosas a la vez:
    1. κ = 500 m² s⁻¹, 5 veces el de Si et al. y 33 veces el Redi de Dettling;
    2. isopicnas empinadas en un talud de pendiente 0.1, con taper dm95 a maxSlope 0.01, que es más bajo que el 0.025 de Si;
    3. DST3 sin limitador, que no es monótono.
  - En la forma skew, la parte diagonal de Redi con flujos cruzados puede producir extremos espurios cuando κ Δt/Δx² y las pendientes son grandes.
  - Antes de aceptar la forma advectiva con GM_Small_Number = 1e-6, probaría bajar κ a 100 (Si) y subir maxSlope a 0.025 con DM95 (Si).
  - Si se usa GM_Small_Number = 1e-6, **tiene que ser el mismo valor en el directo**, y V1 tiene que pasar con ese valor. Si no, el adjunto no es el gradiente del modelo que se reporta.
  - Todo esto es post hoc y debe declararse así.

---

## 4. Caos: Lea et al. (2000) y McLay (2006)

**Lo que dicen:**
- **Lea et al. (2000):**
  - En el resumen: aplicado a trayectorias largas, "cumulative error growth means that adjoint results diverge exponentially from the 'macroscopic climate sensitivity'"; aplicado a trayectorias cortas, "the adjoint provides an incorrect estimate of the sensitivity".
  - En una ventana intermedia, "an ensemble of adjoint gradients can give a reasonably accurate (O(10%)) estimate".
  - Usan ensambles de 199 a 1314 segmentos (Fig. 7, p. 530), y advierten que "this ensemble-adjoint approach is unlikely to be reliable for more complex systems".
  - Para el método directo, la estimación converge si "(a) we allow sufficient time for the system to display its true sensitivity, (b) Δr is large enough that the error due to local extrema is small, and (c) Δr is small enough that the underlying curvature … is weak" (p. 528).
  - Los extremos espurios escalan como 1/√t con la longitud del promedio (Fig. 4).
  - Conclusión: "the direct perturbation method is the most reliable" (p. 531).
- **McLay (2006):**
  - El canal es de **fondo plano a 45°N**, a 0.1° y 0.48° (§2.3.1, p. 22). No tiene talud y no es polar, así que no es "el mismo marco" del Exp. 5 salvo por ser un canal reentrante de MITgcm.
  - Tiempo no lineal ≈ 200 días; información útil en la estructura espacial a 278 días.
  - A 0.48° (no caótico, sin GM) el adjunto funciona a 690 días.
  - "Finite difference gradient checks do not provide a good indication of whether useful information remains in gradients … in eddy resolving ocean models. Perturbed forward model experiments provide an alternative check" (§3.5, p. 53).
  - El modelo "amortiguado" usa el esquema de acoplamiento C-D, no una parametrización de remolinos tipo GM (pp. 22–23).
  - Sus experimentos directos son corridas **únicas**, sin ensamble, con perturbaciones de gran escala que duran toda la integración (p. 45).

**(a) ¿Amenaza al modelo grueso con GM en 360 días?**
- No por construcción: el análogo de McLay (0.48°, sin remolinos) es lineal a 690 días.
- Pero hay que **comprobar** que el grueso no es caótico. Dettling (p. 1683): "even at coarse resolution, there is a contribution of resolved eddies at the lower continental slope where the GM coefficient is low … the Rossby radius of deformation is about 30 km so that some eddies are resolved".
- [cálculo propio] Nuestro radio de deformación mar afuera es ≈ 18 km, del orden de 2Δx a 10 km. Hay que reducir GM en el talud (factor 0.1) y vigilar que no aparezcan remolinos resueltos.
- [inferencia] viscAh = 1000 probablemente los mata. Una prueba barata, siguiendo a McLay (§3.2.1, que sigue a Kleist & Morgan):
  - dos corridas directas de 360 días, una de ellas con 10⁻³ °C en T₀;
  - si ΔJ(t) crece exponencialmente, el modelo grueso es caótico y V1–V3 no bastan.
- Es una prueba post hoc y se declara como diagnóstico, no como criterio.

**(b) ¿Amenaza a H5b (fino, N = 3, |ΔJ| > 2σ)?**
- Sí, de tres maneras [inferencia apoyada en las citas]:
  1. **Potencia estadística.**
     - Con N = 3 por grupo, s_C y s_P tienen 2 grados de libertad cada uno, y la regla "2σ" equivale a un t con ≈ 4 grados de libertad: p ≈ 0.1, no 0.05.
     - Lea et al. piden promedios largos (extremos ∝ 1/√t) y Δ "large enough". 30 días de promedio en el fino es corto.
  2. **Equilibrio del fino.**
     - S&T (p. 433): "It is also important that each simulation reaches a statistically steady state, as otherwise the results could be strongly dependent on arbitrary choices made in setting the initial conditions." Ellos necesitaron 5–20 años tras interpolar de 2 a 1 km (p. 436). Si et al. usaron 20 años.
     - F0 de 6–12 meses no llega al equilibrio. Como FC y FP comparten la deriva, la diferencia sigue siendo interpretable, pero hay que decirlo.
  3. **Resolución.**
     - S&T (p. 439): el mecanismo es "underestimated by close to 50% even at 2 km". Si et al.: a 5–10 km "they do not adequately resolve baroclinic eddies over the continental shelf" (p. 5).
     - El fino de 2.5 km no es la "verdad", y una discrepancia entre los dos modelos no falsa por sí sola al adjunto.
- **Lo que recomiendan:**
  - Lea et al.: ensamble de adjuntos de longitud intermedia (cientos), o el método directo con perturbaciones finitas.
  - McLay: perturbaciones directas de gran escala más el enfoque de malla múltiple.
  - Ninguno da un N mínimo para ensambles directos [no verificado: Lea et al. 2002 no está local].

---

## 5. Verificación de citas

| Afirmación (archivo) | Veredicto | Evidencia |
|---|---|---|
| "Mazloff (2011)", J. Climate, sensibilidad de Drake (spec §1, §4, §7; PD0; rev §1.5) | **Año incorrecto** | Mazloff, M. R., 2012: *J. Climate*, 25, 2279–2290, doi:10.1175/JCLI-D-11-00030.1 ("1 APRIL 2012"; recibido en ene. 2011, versión final oct. 2011) |
| "Stewart & Thompson (2014)" y "(2015), citado por Bull et al. 2025, ~1 km" (rev §1.2; spec §3.2) | **Es el mismo artículo, GRL 2015** | Su línea de cita dice "(2015) … Geophys. Res. Lett., 42, 432–440". Thompson et al. (2018) lo listan como 2015a. El "1 km" está en p. 433 y p. 437 ("no further increase in F_CDW for Δx ≤ 1 km"). Que Bull et al. citen este artículo: [no verificado] |
| Rev §1.5: Mazloff, "Drake responde al viento en ~100 días; sensibilidad en escalas de ~100 km en topografía compleja; manda el rotacional" | **Parcial** | La parte de 100 km y del rotacional está apoyada (resumen; p. 2289: dipolo "indicative of a sensitivity to wind stress curl"). Los "~100 días" son el equilibrio de la sensibilidad de la media; el impulso decae en semanas (p. 2282) |
| Rev §1.5 y spec §7: Rodriguez, "rotacional sobre el quiebre sostiene el transporte; un adjunto puede probarlo directamente" | **Lo primero está apoyado; la transferencia a nuestro canal, no** | Resumen: "the cross-shelf transport is sustained by vorticity input from the localized wind-stress curl over the shelf break". El balance es βv contra el rotacional y el estiramiento (Ec. 3; −0.31, −1.50, +29.34, −47.35, +16.64, +3.18 Sv, p. 3346). Con β = 0 y batimetría uniforme en x, el término no existe en nuestro canal. Rodriguez no usa adjunto |
| Rev §1.6: L&H, "sensibilidad a la batimetría tan grande como a las condiciones de superficie; esfuerzo de forma" | **Apoyado** | Resumen: "sensitive to bottom topography as much as they are to surface boundary conditions"; "interpreted in terms of bottom form stress" |
| Rev §4: H&L 2012, "30–60 días por advección; rutas sobre un umbral; regiones de sombra" | **Apoyado** | Resumen: "dominant timescales of 30–60 days … via advective processes"; "Dominant sensitivity pathways are found over a sill, as well as 'shadow regions'" |
| Rev §5: Lea, "trayectorias largas divergen, cortas subestiman, ensamble intermedio ~10 %" | **Apoyado** | Resumen. Matiz: Lea et al. dicen que el ensamble "is unlikely to be reliable for more complex systems" |
| Rev §1.2 y §5: McLay, "200 d; útil a 278 d; FD no sirve; parametrización que reduce EKE coincide a 690 d"; "el mismo marco del Exp. 5" | **Apoyado con matices** | Resumen. La "parametrización" es el esquema C-D, no GM (pp. 22–23). El canal es de fondo plano a 45°N (p. 22), no "el mismo marco" |
| Spec §4: "365 días … menos que el tiempo no lineal de McLay, unos 200 días" | **Error lógico** | 365 > 200. Lo defendible: "el grueso no caótico, como el 0.48° de McLay, es útil a 690 días" |
| Spec §3.2 y PD1: GM "~500 m² s⁻¹ baja un orden de magnitud sobre el talud, como en Dettling", aplicado a GM **y** Redi | **Parcial** | Dettling apoya el factor 10 en κ_GM (resumen). No da 500 m² s⁻¹ (el valor sobre el talud es O(15), p. 1680). κ_Redi = 15 uniforme, **sin** reducir (p. 1682) |
| Spec §3.1: "arquetipo de S&T 2014, 2016" | **Parcial** | El dominio coincide. El talud (100 km, 25 km de ancho) y el viento (con westerlies) no |
| Spec §6, V3: "< 30 % (Lösch & Heimbach 2007)" | **Apoyado, pero para 100 años** | "<30%" a 100 años; "order 5% for 1-yr integrations" (p. 385) |
| PD0: σ "Mazloff 2011; Jones et al. 2019" (τx 0.02 N m⁻²) | **No apoyado por Mazloff** | Mazloff usa 0.1 N m⁻² (Fig. 3) y no escala T₀, S₀ ni κ_v. Lo más cercano a 0.02 es el "50 % de incertidumbre" de L&H (p. 390). Jones et al. (2019): [no verificado] |
| v4: "más viento del este … hunde las isopicnas y enfría la plataforma (S&T 2012; Morrison 2023)" | **Apoyado por S&T 2015** (los otros dos no están locales) | S&T 2015, p. 437: "Increasing τ_max steepens the isopycnals … deepening the pycnocline at the shelf break". Matiz: en S&T es el esfuerzo **sobre el talud**. Nuestra variante SW sube τ uniforme en y < 150 km y no cambia el rotacional sobre el talud [inferencia] |
| v4: "la plataforma fresca necesita agua superficial dulce (Thompson 2018)" | **Débil** | Thompson et al. atribuyen el frente al viento y al bombeo de Ekman (p. 748); los flujos de flotabilidad son "less clear" (p. 750). Si et al. sí fijan la frescura de la plataforma, pero con toda la columna en la frontera sur, y encuentran **más** calor hacia la costa al refrescarla |
| Rev §2: "el CDW cruza solo por remolinos; la superficie por Ekman" (S&T) | **Apoyado** (régimen denso) | Resumen; p. 436 |
| Rev §2: "más deshielo da más flujo hacia la costa" (Si et al.) | **Apoyado** | Resumen; p. 3 |
| Spec §3.3: "la salinidad importa en el frente (Si et al.)" | **Apoyado** | Título; p. 3: "largely controlled by the magnitude of the cross-slope salinity gradient" |

---

## 6. Recomendaciones, por impacto

### A. Apoyadas por los artículos (con cita)

1. **Poner el viento con convergencia de Ekman sobre el talud y quitar o alejar los vientos del oeste.**
   - Un perfil solo del este, que se anule cerca de la costa y tenga el máximo cerca de 25 km mar afuera del talud, como en S&T (p. 434) y Dettling (Tabla 1). O lineal desde la costa, como en Si et al. (p. 8).
   - Cita: "with a tendency for the wind stress curl to induce downwelling velocities (Ekman pumping) over the continental slope" (Thompson et al. 2018, p. 750); "the local along-shore wind stress controls the pycnocline depth at the shelf break" (S&T, p. 438).
2. **Bajar κ_Redi a O(15–100) m² s⁻¹ y aplicar la reducción del talud solo a κ_GM.**
   - Dettling: κ_Redi = 15 uniforme (p. 1682); con una difusividad alta, "the onshore heat fluxes are strongly overestimated" (p. 1684). Si et al.: 100 m² s⁻¹ con pendiente máx. 0.025 y DM95 (p. 9).
   - Probablemente cura también la inestabilidad del directo en skew-flux (§3) [inferencia].
3. **Darle a la plataforma un sumidero de calor o de flotabilidad en toda la columna.** Hay dos opciones:
   - (a) plataforma densa como S&T y Dettling: sal de polinia de 2.5×10⁻³ g m⁻² s⁻¹ en 50 km. Es el único régimen donde un GM sensible al talud está validado (Dettling);
   - (b) plataforma fresca como Si et al.: esponja sur (≈ 20 km) con θ = θ_f y S uniforme y más dulce que la superficie de mar afuera.
   - Cita de S&T (p. 434): "necessary to impose this additional, explicit sink of buoyancy on the continental shelf in order to achieve a statistically steady state".
4. **Revisar el criterio V4 con Si et al.**
   - En la plataforma fresca sin marea, J total ≈ 0 (p. 6, fig. S9A).
   - "J > 0" y "plataforma fresca" pueden estar en tensión. Al ser un criterio congelado, no se cambia, pero el README debe citar este resultado al interpretar un fallo.
5. **Advertir que GM falla en la plataforma fresca:** "The gray lines in fig. S2 … show that the straightforward application of the GM-Redi parameterization cannot correctly capture the onshore heat transport" (Si et al., p. 9). Si se queda en plataforma fresca, H5b es una prueba de esa limitación.
6. **Ensanchar el talud** (W ≥ 50 km, para tener ≥ 5–15 celdas a 10 km). Rango de S&T: 50–250 km. Dettling: 15 celdas.
7. **Esponja norte de 50 km que relaje también la velocidad** (S&T, p. 434; Dettling, p. 1677).
8. **J: reportar la descomposición** en media, remolinos resueltos, GM y Redi (Dettling, Ecs. 6–7 y 11–13; Si et al., Ecs. 5–8), y además J relativa a −1.8 °C (S&T, p. 438). Usar bandas de y como F_shelf (50–75 km) y F_slope (125–150 km) de Si et al. como diagnóstico.
9. **Corregir las citas:** Mazloff (2012); S&T (2015) como un solo artículo; σ de Mazloff = 0.1 N m⁻²; L&H 30 % a 100 años y ~5 % a 1 año.
10. **Citar a McLay (pp. 84–85) como precedente del diseño de dos modelos**, con la advertencia "processes governing the sensitivity … are not necessarily the same in models of different resolution" (p. 81).

### B. Inferencias mías

1. **H5a.** Replantearla como bombeo de Ekman sobre el talud, que cambia la profundidad de la picnoclina en el quiebre. Es el mecanismo "L_wind" de S&T (p. 437). En nuestro canal f/β = 0, Rodriguez y Mazloff no aplican de forma literal.
2. **Predictibilidad del grueso.** Hacer un par de corridas gemelas (10⁻³ °C) antes de dar por buenos V1–V3 en 360 días.
3. **H5b.**
   - Usar pares de miembros con la misma semilla para control y perturbado, y evaluar la diferencia por pares.
   - Promediar J en el fino sobre más de 30 días (p. ej., los últimos 90).
   - Declarar que con N = 3 la potencia es baja.
   - Esto no cambia los umbrales; añade diagnósticos post hoc declarados.
4. **GM_Small_Number = 1e-6.** Usarlo igual en el directo y en el adjunto, y repetir V1 con ese valor.
5. **Viscosidad.** viscAh = 1000 deforma el ASC. Dettling corre a 10 km y Δt = 1200 s con A_h 12 más biarmónica. Si el directo es estable tras arreglar viento y Redi, conviene bajarla y declararlo.
6. **Lo que objetaría Mazloff:**
   - (i) la cita;
   - (ii) que las σ no son las suyas;
   - (iii) J de 30 días en vez de la media de la ventana;
   - (iv) "the more nonlinear the solution, the shorter time window the adjoint model will be relevant" (p. 2281);
   - (v) que la resolución cambia las sensibilidades: "increasing to 1/6° resolution results in the sensitivities becoming much more localized, and even, in some cases, changing signs" (p. 2280). Esto justifica H5b;
   - (vi) un canal de plano f no puede probar el mecanismo de Rodriguez;
   - (vii) "rbcs", la relajación superficial y la ausencia de hielo hacen el forzamiento poco realista para el Amundsen, que es plataforma cálida, fuera del régimen de Si et al. (p. 9, limitación 5).
7. **Opción de máximo crédito con poco tiempo.** Reproducir la configuración pública de Dettling (zenodo 10.5281/zenodo.10033249) en el grueso. Daría un V4 cuantitativo contra un artículo publicado, en régimen de plataforma densa. Su viabilidad en los 2 días que quedan: [no verificado].

---

## 7. Huecos por artículos ausentes

- **S&T 2012 y Morrison et al. 2023:** la sensibilidad de la plataforma al viento del este está cubierta por S&T 2015, pero no la de la plataforma fresca al viento.
- **S&T 2016:** es la configuración base de Dettling ("closely resembles the one described in Stewart and Thompson (2016)", p. 1677). Sin él no se verifican los perfiles exactos de T y S.
- **Si et al. 2022 (ref. 52 de Si 2023):** tiene los perfiles, el viento en N m⁻² y la tabla S1. Sin él, los valores numéricos de salinidad de la plataforma fresca quedan [no verificado].
- **Stewart et al. 2018:** la isóbata de 1000 m como línea estándar y la compensación entre marea y media (citado vía Thompson et al. y Si et al.).
- **Bull et al. 2025 y Jones et al. 2019:** la "respuesta lineal de primer orden" y las σ de Jones quedan [no verificado].
- **Gaikwad et al. 2024 y Forget et al. 2015:** el criterio del 1 % y la práctica de ECCO con GM en el adjunto (forma, taper, viscosidad del adjunto) quedan [no verificado]. Ningún artículo local cubre GM_Small_Number ni la forma skew contra la advectiva.
- **Heimbach & Losch 2012, pp. 62–64:** la ventana, la resolución y las pruebas de gradiente quedan [no verificado].
