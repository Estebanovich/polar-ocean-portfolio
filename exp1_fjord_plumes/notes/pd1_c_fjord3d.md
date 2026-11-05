# PD1 — configuración (c) `c_fjord3d` (Task 7, Step 4; sin correr)

Fiordo 3D hidrostático: pluma subglacial parametrizada (teoría de plumas →
sumideros y fuentes de `addMass` + RBCS, opción (iii) de PD3 aprobada), pulso de
viento a lo largo del fiordo y oscilación de la picnoclina en la boca. Cuatro
casos: `c_fjord_plume`, `c_fjord_wind`, `c_fjord_shelf`, `c_fjord_all`, de 30
días cada uno. Build MPI (5 rangos):
`/Volumes/ESTEBAN_MIT/portfolio_builds/exp1/c_fjord3d/`. **No se ha corrido
MITgcm** (la campaña de (a) sigue en curso). El Step 5 (corrida cronometrada de
1 día) viene después.

## Parámetros

| Parámetro | Valor | Fuente / razón |
|---|---|---|
| Dominio | 40 km × 6 km × 500 m; 200 × 30 × 50 celdas | `cases.yaml` (Task 3) |
| Celdas | 200 m × 200 m × 10 m | Cowton et al. 2015 §3.2: la misma malla (200 m a lo ancho, 10 m en la vertical) |
| Paredes | filas de tierra j = 0 y Ny−1 (ancho mojado 5.6 km); columna Nx−1 de tierra, OB en Nx−1 (`OB_Ieast=-2`) | lecciones de (a)/(b) |
| Descomposición | `nPx = 5`, `sNx = 40`, `sNy = 30` (np = 5) | prueba `test_c_grid_matches_size_h_and_5_ranks` |
| Hidrostático | `nonHydrostatic=.FALSE.`, `#undef ALLOW_NONHYDROSTATIC` | 200 m ≫ H/aspecto; pluma no resuelta |
| f | 1.37e-4 s⁻¹ | Cowton et al. 2015, Tabla 1 |
| EOS | lineal; tAlpha 3.87e-5, sBeta 7.86e-4, rhoNil = rhoConst = 1027, Cp 3974 | paridad con (a)/(b) |
| viscAh, diffKh(T,S) | 30 m² s⁻¹ (Laplaciano) | Cowton et al. 2015, Tabla 1 |
| viscAr, diffKr(T,S) | 1e-5 m² s⁻¹ (implícitos) | Cowton et al. 2015, Tabla 1 |
| Convección | `ivdc_kappa = 10` m² s⁻¹ | hidrostático: mezcla implícita en columnas inestables (`calc_3d_diffusivity.F:96-110`) |
| Paredes / fondo | sin roce lateral (free-slip); arrastre cuadrático 2.5e-3 en el fondo | 200 m no resuelve capas límite laterales; Cd = `MeltParams.Cd` |
| Advección | `tempAdvScheme = saltAdvScheme = 33`, `staggerTimeStep` | como (a)/(b) |
| dt | 60 s (provisional; lo confirma el Step 5) | ver estabilidad |
| Duración | 30 d = 43 200 pasos | plan |
| Salida | MDS float64 big-endian; medias de 3 h de UVEL y THETA (Q_ex, calor); instantáneas diarias de U, V, W, T, S; ETAN horario | Fig. 3 |
| Pickups | cada 5 d (y al final de la etapa 1 de los casos con viento) | — |
| icefront | **no** | ver abajo |

## Pluma (PD3, opción iii)

- `plume_source_sink` (semicono, modo `mitgcm`, α = 0.1, Q0 = 100 m³/s,
  T0 = S0 = 0) en el perfil Sermilik de (c) da:
  - z_neutral = −181 m;
  - arrastre total 3625 m³/s en las celdas k ≥ 19 (todas las celdas bajo la salida);
  - salida Q(z_lo) = 3725 m³/s a 3.42 °C y 33.53 psu.
- `addMass` (kg/s) va en la columna de la cara (i = 0).
- **Reparto (regla R18):** j = 12–17 (6 celdas a lo ancho, 1200 m), con un
  sexto cada una. La salida va en **un solo nivel**, k = 18 (la celda de
  z_neutral, como en Cowton et al. 2015). `n_plume_cols = 6`,
  `n_out_levels = 1` en `cases.yaml`.
- **Neto:** Σ addMass = ρ0·Q0 exacto (C5). Queda en `case.json`
  (`net_addmass_over_rho0_Q0 = 1.0`).
- **Agua local en sumideros y fuentes:** `temp_addMass` y `salt_addMass` sin
  definir, y `salt_EvPrRn = 1.234567E+05` (UNSET_RL; si no, S de `addMass` = 0).
- **T/S de la salida:** RBCS solo en las 6 celdas de salida (k = 18, j = 12–17).
  - `tauRelaxT = tauRelaxS = 600` s y máscara = q·τ/V = 0.93 (≤ 1).
  - Objetivos `rbcs_T` y `rbcs_S` = 3.42 °C / 33.53 psu ahí; perfil inicial
    fuera (sin efecto, máscara 0).
  - La tendencia RBCS es igual a q(T_out − T)/V (spike PD3: 7e-5 °C).
- Casos sin pluma (`wind`, `shelf`): addMass = 0 y máscara 0. Los namelists
  son idénticos.

### Recambio de la celda de salida: elección

- **Una sola celda** (200 m × 200 m × 10 m, V = 4e5 m³) recibiría 3725 m³/s:
  - recambio V/q = 107 s, es decir λ·dt = 0.56 con dt = 60 s, cerca del límite
    de AB2 para un término de relajación explícito (λ·dt < 1);
  - salida de ~1.9 m/s por una sola cara de 2000 m² (CFL 0.56).
  - Habría que bajar a dt ≈ 20–30 s, lo que duplica o triplica el costo.
- **Elegido (R18): 6 celdas a lo ancho × 1 nivel.**
  - Por celda: q = 621 m³/s, recambio 644 s, λ·dt = 0.093.
  - Salida por las 6 caras en x (12 000 m²): ≈ 0.31 m/s (CFL 0.09).
  - La primera versión (Step 4) usaba 3 × 2 celdas, con los mismos números.
    La revisión encontró que la mitad puesta en k = 17 es **0.068 kg/m³ más
    densa** que el ambiente ahí: una inversión. En k = 18 es 0.049 kg/m³ más
    ligera. Cowton et al. 2015 terminan la pluma en la flotabilidad neutra
    justamente para evitarlo. Ahora solo se reparte lateralmente; la prueba
    `test_c_outflow_one_level_six_columns_not_denser_than_ambient` lo vigila.
  - Justificación: el radio del semicono en z_neutral es b ≈ 45 m (teoría),
    menor que una celda, pero la intrusión se esparce lateralmente por gravedad.
    Repartir en 1200 m a lo ancho (21 % del ancho mojado) es una elección
    numérica (clase (c)), no una resolución física de la intrusión.
  - Los sumideros se reparten en las mismas 6 columnas. Así el transporte
    vertical de la pluma ocurre en una sola columna de 1200 m de ancho y no se
    crea una recirculación horizontal artificial entre columnas de sumidero y de
    salida.

### Diferencias conocidas, clase (b) (física faltante; regla R19)

1. **Forzamiento de la pluma fijo.** z_neutral, T_out, S_out y los sumideros se
   calculan una vez, con el perfil inicial. Cowton et al. 2015 los recalculan
   en cada paso con el perfil del modelo. Consecuencias:
   - En `shelf` y `all`, la salida **no responde** a la oscilación de la
     picnoclina: no cambia su profundidad ni sus propiedades.
   - **Calor y sal no se conservan estrictamente.** Los sumideros retiran agua a
     la T/S local, que evoluciona, y la fuente la añade a la T_p fija. El
     desbalance es ∝ Q_arrastre·(T_local(t) − T_local(0)). Es cero mientras la
     columna de la cara conserve su perfil inicial.
2. **ΣQ = Q0 incluye el derretimiento.** La teoría añade el agua de deshielo de
   la pared del semicono (≈ 0.39 m³/s bajo la salida) al flujo de la pluma. Para
   que Σ addMass = ρ0·Q0 sea exacto, los sumideros retiran del ambiente
   arrastre + deshielo, y el deshielo no es agua del ambiente. Es 1e-4 del
   arrastre. Queda en `case.json` (`melt_volume_in_sinks_m3_s`, `sum_rule_note`).

## Forzamientos

- **Viento** (`wind`, `all`): τx = +0.2 N/m², a lo largo del fiordo y hacia
  afuera (down-fjord: Jackson et al. 2018 §1, "stronger down-fjord winds from
  the ice sheet"), durante 12 h.
  - Uniforme en el dominio, en `windx.bin`: registros horarios centrados
    (`get_periodic_interval.F`), 12 registros con τ y el resto en 0. El ciclo es
    de 6 d, compartido con la OB.
  - Rampas de media hora al inicio y al final por la interpolación lineal.
  - **Protocolo de 2 etapas de la tesis** (`verification/BayIW_Rectan/README.md`,
    "Two-stage run strategy"):
    - etapa 1 (`data`): 0–24 h con `zonalWindFile`, pickup a las 24 h
      (iteración 1440);
    - etapa 2 (`data.stage2`): `nIter0 = 1440` hasta los 30 d, sin
      `zonalWindFile`. Sin esto, el pulso se repetiría cada ciclo (6 d).
  - El runner de la etapa 2 está en `scripts/run_campaign.py`
    (`run_stage2`). Si la etapa 1 está OK, copia (nunca enlaza) las entradas y
    todos los `pickup*.<n1>*`, incluido `pickup_orlanskiE` (lección C3 de (a)),
    a `<run>/input_stage2` y corre en `<run>/stage2`. Cada etapa tiene su
    `run_status.json`, y `<run>/run_status_final.json` está OK solo si ambas
    lo están y last_iter = nIter0 + nTimeSteps de la etapa 2.
  - Nota: en `BayIW_Rectan` el protocolo está en su `README.md`.
- **Plataforma** (`shelf`, `all`): picnoclina en la boca con
  η(t) = A sin(2πt/P), A = 30 m y P = 6 d. Jackson et al. 2018: "hS of 30 m
  (typical for Sermilik Fjord)" y su control de 6 días.
  - `OBEt.bin` y `OBEs.bin`: perfil Sermilik con centro en 200 m − η(t_n), 144
    registros horarios.
  - **OBCS:**
    - u y v por Orlanski (como (a)/(b));
    - T y S **prescritas** (`useOBCSprescribe`, `ALLOW_OBCS_PRESCRIBE`).
      `OBCS_PRESCRIBE_READ` corre después de `ORLANSKI_EAST`
      (`obcs_calc.F:99, 413`) y solo sobrescribe los campos con archivo
      (`obcs_fields_load.F:159-174, 518-524`);
    - esponja de T/S de 10 celdas (2 km), τ de 600 s en la OB a 6 h en el borde
      interno (`obcs_sponge.F:553-559`). Es la "relaxation zone" de Cowton et
      al. 2015 §3.2, aquí dependiente del tiempo.
  - **Velocidades no prescritas:** fijar u en la boca bloquearía el
    intercambio baroclínico de la pluma (~3700 m³/s de salida) y del viento, y
    el neto Q0. `useOBCSbalance` sigue apagado: ignora `addMass`. Q0 sale por
    la radiación.
  - En `plume` y `wind`, los registros son el perfil inicial (constante): la
    misma zona de relajación que Cowton.

## Por qué sin pkg/icefront

- El derretimiento propio de la pluma ya está en la teoría (0.4 m³/s, incluido en ΔQ).
- Fuera de la pluma, con γT en su piso (1e-4 m/s, `icefront_thermodynamics.F:94-95`)
  y ΔT ≈ 6 K en la AW, el derretimiento sería ≈ 7e-6 m/s sobre ~1.7e6 m², es
  decir ≈ 10–15 m³/s de agua dulce. Es ≈ 0.3–0.4 % del transporte de salida de
  la pluma y mucho menor que el intercambio forzado por viento o plataforma.
- A 200 m las velocidades en la cara no representan la capa límite, así que ese
  derretimiento sería sobre todo una elección del piso de γT.
- En `wind` y `shelf` (Q0 = 0) sería la única fuente de flotabilidad y
  ensuciaría la partición del intercambio (Fig. 3).
- Se escriben `R_icefront.bin` e `icefrontlength.bin`, pero ningún namelist los usa.

## Sección para Q_ex (Fig. 3)

La esponja de U/V y T/S de OBCS cubre las 10 columnas antes de la OB: OB en
i = 198 (Python, `OB_Ieast = -2` → columna Nx−1 de Fortran), esponja en
**i = 188–197**. Ahí u y v se relajan a los valores de Orlanski y T/S a las de
la OB, así que el intercambio está amortiguado artificialmente. **La sección de
Q_ex y del calor hacia el frente debe estar a media distancia del fiordo**
(p. ej. i = 100, x = 20 km), fuera de la esponja y lejos de la columna de la
pluma (i = 0). También conviene evitar la zona de recirculación de la salida
cerca de la cara.

## Números adimensionales (para la Fig. 3)

- Δρ (PW/AW, EOS lineal) = 1.154 kg/m³, g' = 0.0110 m/s².
- c = √(g' h1 h2/H) = 1.15 m/s (h1 = 200 m, h2 = 300 m).
- **Rd = c/f = 8.4 km**, **W/Rd = 5.6/8.4 = 0.67**, ωL/c = 0.42 (L ≈ 39.8 km).
- Jackson et al. 2018: Sermilik ≈ 0.9 con W = 7 km. Con ny = 37 (35 filas
  mojadas, 7.0 km) quedaría W/Rd = 0.83. **No se cambió `ny`** sin visto bueno.

## Estabilidad con dt = 60 s

| Proceso | Número |
|---|---|
| Onda interna, modo 1: c·dt/dx | 1.15·60/200 = 0.35 |
| Salida de la pluma: u·dt/dx | 0.31·60/200 = 0.09 |
| Recambio de la celda de salida (RBCS/addMass): λ·dt | 60/644 = 0.093 |
| Viscosidad lateral: ν·dt/dx² | 30·60/200² = 0.045 |
| Esponja OBCS: dt/τ en la celda más externa | 60/2700 = 0.022 |
| Superficie libre | implícita |

## Costo a priori (escalado de (a): 2.8e-6 s por celda·paso, serial)

| | Valor |
|---|---|
| Celdas | 200·30·50 = 3.0e5 |
| Pasos (30 d) | 43 200 → 1.30e10 celda·pasos |
| Serial | 36 300 s ≈ **10.1 h por caso** |
| 5 rangos (factor 3.3–5, el de (b)) | **2.0–3.1 h por caso**; **8–12 h los 4 casos** |
| Corrida corta del Step 5 (1 d = 1440 pasos) | 1210 s serial ≈ 4–6 min con 5 rangos |
| Disco | ≈ 1.7 GB por caso, ≈ 7 GB en total en LaCie |

- **Cota superior:** el 2.8e-6 de (a) es no hidrostático, dominado por cg3d
  (300–500 iteraciones). En (c) no hay cg3d, así que se espera 3–10 veces menos.
- **Presupuesto del Exp. 1 (≤ 24 h):** (a) 2.6 h + (b) 8.5–12 h + (c) 8–12 h
  = 19–27 h en el peor caso. **Si el Step 5 confirma > 2 h por caso, se para en
  PD1.**
  - Opciones: dt = 90 s (onda interna 0.52); 20 d en vez de 30; o un solo caso
    `all` más los de un solo forzamiento acortados.

## Comandos (Steps 5–7)

```bash
cd ~/Phd_project/portfolio
PY=/opt/anaconda3/envs/pyocean/bin/python
# Step 5: corrida cronometrada de 1 día (ejercita pluma, viento, OB y el
# reinicio de la etapa 2: etapa 1 de 12 h con viento, etapa 2 de 12 h)
$PY exp1_fjord_plumes/scripts/run_campaign.py --config c_fjord3d --cases c_fjord_all \
    --np 5 --end-time 86400 --set wind_stage1_s=43200 --run-suffix _short1d \
    > /Volumes/LaCie/portfolio_runs/exp1/c_short1d_campaign.log 2>&1
# Step 6: campaña de 4 casos de 30 d (tras PD1; --skip-ok para reanudar)
$PY exp1_fjord_plumes/scripts/run_campaign.py --config c_fjord3d --cases 'c_fjord_*' \
    --np 5 --skip-ok > /Volumes/LaCie/portfolio_runs/exp1/campaign_c.log 2>&1
# Step 7: Fig. 3 (+ C5') y GIF
$PY exp1_fjord_plumes/analysis/fig3_exchange_partition.py
$PY exp1_fjord_plumes/analysis/animate_c.py /Volumes/LaCie/portfolio_runs/exp1/c_fjord_all --check-scales
for c in plume wind shelf all; do
  $PY exp1_fjord_plumes/analysis/animate_c.py /Volumes/LaCie/portfolio_runs/exp1/c_fjord_$c
done
```

- Disco adicional por la animación: la lista `animV` (VVEL en el nivel 19,
  medias de 3 h) suma 240 × 48 kB ≈ 12 MB por caso, ≈ 46 MB en total. Los
  cuadros reutilizan `mean3h` (medias de 3 h de U y T 3D), que ya existía para
  la Fig. 3.
- GIF: ≈ 0.07 MB por cuadro con datos sintéticos ruidosos a tamaño real. Con
  240 cuadros la zancada sube sola a 2–4 para quedar ≤ 8 MB.

## Pendiente (Steps 5–6)

- Corrida cronometrada de 1 día (`c_fjord_all`, que ejercita todo: pluma,
  viento etapa 1, OB). Revisar en el log:
  - `OBCS_CHECK` OK;
  - lectura de RBCS;
  - `salt_addMass = 1.234567E+05`;
  - CFL;
  - monitor de T/S (C5': T y S dentro del rango inicial y de la fuente).
