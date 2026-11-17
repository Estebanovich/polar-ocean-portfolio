# PD1 — configuración del canal (Exp. 5)

Estado: **provisional**, para revisión de Esteban. Figura de control: `figures/pd1_section.png`.

## Malla y batimetría
- Canal reentrante en x, Lx = 400 km, Ly = 450 km, plano f con f₀ = −1.4×10⁻⁴ s⁻¹.
- `chan10`: 40 × 45 columnas (Δx = 10 km). `chan2p5`: 160 × 180 (Δx = 2.5 km).
- Vertical: Nr = 35, con espesores que crecen linealmente de 10 m a 161 m (suma 3000 m).
- H(y) = 500 + 1250·(1 + tanh((y − 100 km)/12.5 km)) m: plataforma de 500 m, talud de ~25 km centrado en 100 km, fondo de 3000 m. Igual en todas las columnas x.

## Línea de J
- y_qb = 90 km: es cara v en las dos mallas (j_qb = 9 en `chan10` y 36 en `chan2p5`, base 0; en Fortran, fila j_qb + 1).
- La profundidad de H(y) en y = 90 km es 920 m; la isóbata de 1000 m cae en y ≈ 91.3 km. Se declara que la línea está "en la isóbata de ~920 m, sobre el quiebre", no exactamente en la de 1000 m.

## Perfiles iniciales [verificar]
- **Fuera de la plataforma:** θ de −1.8 °C en superficie a 0.8 °C (CDW) bajo una termoclina tanh en 300 m (ancho 100 m); S de 34.2 a 34.7 con la misma forma.
- **Plataforma:** θ = −1.85 °C uniforme; S de 34.3 a 34.5 con la misma forma.
- La transición en y usa la misma tanh del talud. Cada columna es una combinación convexa de dos perfiles estables con EOS lineal (α = 3×10⁻⁵ K⁻¹, β = 7.4×10⁻⁴), así que es estable.
- **Ruling:** los artículos que el spec pide citar para estos valores (Stewart & Thompson 2014/2015; Si et al. 2023) no están en la biblioteca local. Para no frenar el día 1 se usan estos valores idealizados, que siguen la estructura de esos modelos ("plataforma fresca", Thompson et al. 2018). Esteban los descarga y en el PD1 se confirman o se cambian antes de GA/GB. Cambiarlos obliga a repetir G0.

## Forzamiento
- τx(y): −0.05 N m⁻² en y < 150 km, +0.1 N m⁻² en y > 250 km, rampa de coseno entre ambos. τy = 0.
- Relajación superficial de T y S hacia el nivel superior inicial (30 días).
- Esponja `rbcs` en las últimas 30 km: máscara de 0 a 1 hacia la pared norte, hacia los perfiles iniciales.
- GM: factor 1 − 0.9·exp(−((y − 100 km)/25 km)²): 0.1 en el centro del talud, 1 fuera.

## `chan10`: numérica y cronometría (Tarea 2)
- Δt = 1200 s (al inicio 1800 s; ver abajo). Con 1800 s y 30 días, CFL advectivo máximo: u 0.035, v 0.022, w 0.015; |u| ≤ 0.2 m/s.
- Advección de T y S con DST3 sin limitador (`tempAdvScheme = 30`): es lineal y diferenciable, y la misma en el directo y en el adjunto. Los limitadores de flujo (p. ej. el esquema 7 del tutorial) no son diferenciables en los extremos. Precio: pequeños rebasamientos (θ máx. 0.85 °C contra 0.80 °C inicial).
- viscAh = 1000 m² s⁻¹ (al inicio 200), viscAr = 1×10⁻⁴ m² s⁻¹ implícita, arrastre de fondo cuadrático 2.5×10⁻³, paredes con deslizamiento libre.
- **Cambio en la Tarea 3:** con viscAh = 200 y Δt = 1800 s, G0 divergió en el año 2 (chorro de talud de 0.44 a 1.5 m/s, CFL w ~0.9, luego NaN). Con viscAh = 1000 y Δt = 1200 s, T2Y y los 30 años de G0 son estables (|u| ≤ 0.15 m/s en T2Y).
- **GM en forma skew-flux** (`GM_AdvForm = .FALSE.`, post hoc, 2026-11-16): en forma advectiva el adjunto de Tapenade diverge (×10¹³ por paso con el `GM_Small_Number` por defecto; ×300 con 1e-8). J incluye entonces el flujo meridional completo de GM/Redi. Ver `v4_revision_posthoc.md`.
- **Cronometría:** 30 días de modelo (T10) en 9 s de reloj con 5 rangos. G0 (30 años, 30 etapas de un año) ≈ 1 h.
- `DIAGNOSTICS_SIZE.h` con `numDiags = 12*Nr` (el valor por defecto, 1*Nr, no alcanza para los 6 campos 3-D).

## Pendiente (Tarea 4)
- Resolución final del fino.
