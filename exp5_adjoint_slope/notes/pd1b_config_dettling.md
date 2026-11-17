# PD1b: el canal pasa a la configuración de Dettling et al. (2024) y Stewart & Thompson (2015)

Fecha: 2026-11-16. Decisión de Esteban, tomada después de la revisión de los artículos (`revision_articulos.md`). **Post hoc**: se decide después de ver fallar V4 y el adjunto, y el README lo declara. Reemplaza la configuración del PD1 (`pd1_configuracion.md`) y la revisión de forzamiento de `v4_revision_posthoc.md` §2: las variantes S0/SW/SF/SWF ya no se corren.

## Por qué
La configuración del PD1 se apartaba de la literatura en todos los ejes que controlan el intercambio a través del talud:
- talud de 25 km de ancho, cinco veces más empinado que el rango de S&T;
- vientos del oeste mar afuera, con surgencia junto al talud;
- ningún sumidero de flotabilidad en la plataforma, que S&T declaran necesario para un estado estacionario;
- κ_Redi = 500 contra 15–100 m²/s.

Dettling et al. (2024) es el único trabajo que valida un GM grueso de 10 km, sensible al talud, contra una referencia de 1 km en este tipo de canal. Copiar su configuración da un V4 comparable con un artículo publicado.

## Valores (fuente entre paréntesis; Tabla 1 = Dettling et al. 2024, Tabla 1)

| Elemento | Valor | Fuente |
|---|---|---|
| Dominio | 400 × 450 km × 3000 m, periódico en x, **cerrado en y** | Tabla 1; Dettling p. 1677 |
| f₀, β | −1.31×10⁻⁴ s⁻¹, 0 | Tabla 1 |
| Batimetría | H = 500 + 1250·[1 + tanh((y − 200 km)/75 km)] | Tabla 1 (H_s, H, Y_s = 200 km, W_slope = 150 km, W_shelf = 125 km). La forma exacta del tanh ("tanh-like", S&T p. 434) **[verificar]**: escala W/2 = 75 km, coherente con la plataforma de 125 km y con la pendiente media ≈ 0.02 (Dettling p. 1676) |
| Viento | solo del este, sinusoidal por tramos: 0 en y = 50 y 450 km, máximo −0.075 N m⁻² en y = 225 km | S&T p. 434; Tabla 1. Forma de cada tramo (cuarto de seno) **[verificar]** |
| Sal de polinia | 2.5×10⁻³ g m⁻² s⁻¹ en y < 50 km (`saltFlux` = −2.5×10⁻³: positivo es hacia arriba, `FFIELDS.h`) | Tabla 1; S&T p. 434 |
| Esponja norte | 50 km; T y S hacia los perfiles iniciales en 54 días; velocidad hacia 0 en 27 días; tasa lineal de 0 en y = 400 km al máximo en la pared | Tabla 1; S&T p. 434 |
| GM | κ_GM = 130/(1 + 800·s_topo) m²/s (Ec. 24, ε_c = 800; κ^{const,high} ≈ 130, leído de la Fig. 5); κ_Redi = 15 uniforme; taper de Gerdes et al. (1991) (`gkw91`) | Dettling Ec. 24, Fig. 5, pp. 1678 y 1682 |
| Viscosidad y mezcla | A_h = 12, biarmónica `viscA4Grid` = 0.1, Leith biarmónica 1.0 (vorticidad y divergencia), A_z = 3×10⁻⁴, κ_z = 5×10⁻⁶, arrastre lineal 1×10⁻³ m/s | Tabla 1 |
| Δt | 1200 s | Tabla 1 (10 km) |

## Diferencias que se mantienen (declaradas)
- **Hielo marino.** S&T y Dettling usan 1 m de hielo fijo con flujos de dos ecuaciones (Schmidt et al., 2004; γ = 6×10⁻⁵ m/s). Aquí:
  - el flujo de calor se aproxima relajando la temperatura superficial a θ_f = −1.9 °C con τ = Δz₁/γ = 10 m / 6×10⁻⁵ m s⁻¹ ≈ 1.9 días, que es la misma conductancia;
  - el agua dulce de fusión se aproxima relajando la SSS a 34.2 en 30 días, salvo en la zona de polinia (y < 50 km), donde no se relaja.
- **Ecuación de estado lineal** (α = 3×10⁻⁵, β = 7.4×10⁻⁴) en vez de McDougall et al. (2003): sin cabbeling ni termobaricidad.
- **Sin KPP:** `ivdc` (κ = 10) para la convección. **35 niveles** de 10 a 161 m, en vez de 77 de 13 a 100 m.
- **Perfiles.** Thompson & Heywood (2008), que usan S&T y Dettling, no está disponible. Se usan perfiles idealizados uniformes en y:
  - θ de −1.8 a 0.8 °C y S de 34.2 a 34.7, con una termoclina tanh en 300 m (100 m de ancho) **[verificar]**;
  - a diferencia de S&T, el perfil de restauración conserva la capa superficial (AASW).
- **GM en forma advectiva** con `GM_Small_Number = 1e-6` (opción A, `v4_revision_posthoc.md` §4).

## Lo que cambia en el experimento
- **Régimen:** plataforma **densa** (sal de polinia), no plataforma fresca.
  - El criterio congelado V4 pide un "frente de plataforma fresca presente y J > 0". La primera parte **falla por construcción** y se declara así. J > 0 se evalúa igual.
  - Si et al. (2023, p. 9) muestran que un GM estándar a 5–10 km no captura el transporte hacia la costa en plataforma fresca. Dettling valida el GM sensible al talud en plataforma densa. Por eso se eligió este régimen.
- **Línea de J:** y = 150 km, en la isóbata de ≈ 1020 m (cara v de las dos mallas). Conserva el "~1000 m" del spec.
- **H5a:** el dipolo "centrado en el quiebre" se centra en el quiebre de la nueva geometría, y₀ = 125 km (antes 100 km). La motivación se reformula como bombeo de Ekman sobre el talud (S&T p. 437), no como el mecanismo de vorticidad de Rodriguez et al. (2016), que necesita β y contornos f/H obstruidos.
- **El fino (`chan2p5`)** usa la misma geometría y el mismo forzamiento, sin GM.
