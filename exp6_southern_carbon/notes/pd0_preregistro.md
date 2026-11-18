# Exp. 6: prerregistro (PD0)

**Estado: BORRADOR (18 nov 2026).**
- Se congela después de G0/G0s y antes de correr GT, con el visto bueno de Esteban. Ahí se anotan la fecha y el hash del commit.
- Lo marcado *[se confirma con G0]* se revisa con el estado de equilibrio. Si cambia, cambia **antes** de congelar.

## Configuración que entra al prerregistro

- **Canal:** `tutorial_reentrant_channel` con `pkg/dic`, 6 trazadores (DIC, Alk, PO4, DOP, O2, PO4pre), sin hierro.
- **Ciclo estacional:** ±1.5 °C en la SST de relajación, con el mínimo en agosto. La luz es la de 55°S.
- **Intercambio de gas:** viento fijo de 10 m s⁻¹. Salinidad uniforme de 34.5.
- **Pasos de tiempo:** `g50` a 2000 s; `f10` a 1000 s. G0 es acelerado (×4) y le siguen 48 años síncronos (G0s). Ver `puerta_dia1.md`.
- **CO₂:** GC en el valor CMIP6 de 1850. GT con el histórico CMIP6 y luego SSP3-7.0, desde el archivo de D1.
- **Esponja norte:**
  - caso principal: preindustrial, con los perfiles `PTRACERS_ref` de `tutorial_global_oce_biogeo`;
  - GTa: con C_ant y P(z) = exp(−z/1000 m).
- **Ramas:** viento ×1.2 como escalón, en 1900, 1980, 2020 y 2060, sobre GT y sobre GC. Las gruesas duran 25 años; las finas, 15.
- **Miembros del fino:** 3. La regla de recorte de `puerta_dia1.md` no se activa con la cronometría del día 1.

## Métricas

- **Región R:** las filas j = 1–19 de `g50` (y < 1000 km), sin la pared.
  - En T10 (año 10), el 99.97 % del crecimiento de la celda superior de la MOC residual cae en R. *[se confirma con G0]*
- **Flujo F:** media de `DICCFLX` sobre el canal, en mol C m⁻² año⁻¹, positivo hacia el océano.
- **Sensibilidad de una rama:** S(t_b) = media de F_rama − F_padre en los años 6–15 de la rama.
  - S_nat = S_GC y S_ant = S_GT − S_GC.
- **Gradientes**, en la década centrada en t_b del padre y sobre R:
  - G₂₀₀ = [DIC(200 m) − DIC(nivel 1)] / 200 m;
  - ΔDIC₁₀₀₀ = DIC(1000 m) − DIC(nivel 1);
  - ambos en µmol kg⁻¹ (ρ₀ = 1035).
- **Capa de mezcla de invierno de R:** el máximo mensual de la climatología de `MXLDEPTH` en el padre.
  - En T10 es de 671 m (septiembre). *[se confirma con G0: tiene que ser > 100 m]*
- **MOC de R:** máx_l ψ(l, j = 20), con ψ la MOC residual en capas de temperatura (`LaVH1TH`, incluye el bolo de GM).

## Validaciones

- **V6a:** en el último siglo de G0, |Δ flujo medio decenal| < 0.01 mol C m⁻² año⁻¹ y |Δ DIC medio bajo 1000 m| < 1 µmol kg⁻¹, con ≥ 1000 años de trazador. Se reporta además la deriva en G0s.
- **V6b:** en GC y GT, en cada año de rama, el residuo del presupuesto anual de DIC del canal es < 1 % del término de superficie.
  - En T2 (año 2) es de 3×10⁻⁵.
  - En R y en la capa de mezcla se reporta, sin criterio.
- **V6c:** Δ_GC = S_GC(2060) − S_GC(1900). Entra en H6a.
- **V6d:** la EKE mensual del fino cumple la regla del Exp. 5 (los últimos 3 meses dentro del 10 % de los 3 anteriores), con al menos 3 años de spin-up.
  - Tope de 8 años. Si no la cumple, H6c se marca "con remolinos sin equilibrar".
- **V6e:** S_GTa(1900) y S_GTa(2060) junto a las de GT. Descriptiva.

## Hipótesis

- **H6a:** pasa si ΔS = S_GT(2060) − S_GT(1900) > 0 y ΔS > 2|Δ_GC|.
  - Se reportan ΔS en mol C m⁻² año⁻¹ y en % de |S_GT(1900)|, S_nat y S_ant por rama, y el año de rama en que S_GT cambia de signo, si cambia.
  - No hay umbral de magnitud.
- **H6b:** en cada una de las cuatro ramas de GT, la suma *offline* (ec. 1 + ec. 4 + término biológico con PO4; `analysis/decomp.py`) debe quedar dentro de ±25 % de la suma medida de circulación más biología.
  - La suma medida es el cambio de rama menos padre en los términos `adv` + `mix` + `bio` del presupuesto de DIC en la capa de mezcla de R, años 6–15, en mol C año⁻¹.
  - Pasa solo si las cuatro ramas cumplen.
  - Se reporta cada término. El término de temperatura se calcula y se reporta aparte, fuera de la suma.
- **H6c:** con ΔS_fino = media de miembros de S(FT_2060) − media de S(FT_1900), pasa si:
  1. |ΔS_fino| > 2σ, con σ = √(s₁₉₀₀²/N + s₂₀₆₀²/N);
  2. tiene el mismo signo que ΔS del grueso;
  3. el cociente ΔS_fino / ΔS_grueso está en [0.5, 2].
  - Con N < 3, H6c es "indicativa".

## Reglas

- Ningún umbral, ventana ni región cambia después de congelar.
- Los fallos se clasifican: (a) nuestro error, (b) física ausente, (c) elección numérica, (d) error de la referencia.
- Lo que se añada después de congelar se rotula *post hoc*.
