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
| Salida | MDS float64 big-endian. Medias de 3 h de UVEL, THETA y **ADVx_TH** (calor); instantáneas **horarias** de UVEL y THETA 3D (Q_ex); instantáneas diarias de U, V, W, T, S; ETAN horario; VVEL en el nivel 19 (GIF) | Fig. 3, regla R31 |
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
| Disco (regla R31) | ≈ 5.8 GB por caso, ≈ 23 GB en total en LaCie (3.6 TB libres) |

- **Extrapolación del Step 5 (regla de la revisión):** se usa la tasa por paso
  **sin los arranques**. El temporizador `MAIN LOOP [THE_MAIN_LOOP]` de MITgcm
  en cada etapa excluye la inicialización y la lectura de entradas;
  `c_runs.step_rate(run)` lo suma sobre las etapas y lo divide entre los pasos.
  - Costo de un caso de 30 d = s/paso × 43 200 + arranques (≈ 1–2 por caso).
  - No se usa `wall_s` de `run_status.json`, que incluye el arranque y, en la
    corrida de 1 d, dos arranques.
  ```bash
  $PY -c "from exp1_fjord_plumes.analysis.c_runs import step_rate; \
  r = step_rate('/Volumes/LaCie/portfolio_runs/exp1/c_fjord_all_short1d'); \
  print(r, 'h por caso de 30 d:', r['s_per_step'] * 43200 / 3600)"
  ```

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
# escalas: --check-scales en LOS CUATRO casos antes de fijarlas (iguales para todos)
for c in plume wind shelf all; do
  $PY exp1_fjord_plumes/analysis/animate_c.py /Volumes/LaCie/portfolio_runs/exp1/c_fjord_$c --check-scales
done
for c in plume wind shelf all; do
  $PY exp1_fjord_plumes/analysis/animate_c.py /Volumes/LaCie/portfolio_runs/exp1/c_fjord_$c
done
```

- **Diagnósticos (regla R31).**
  - **Calor hacia el frente:** media de 3 h de `ADVx_TH`, el flujo advectivo
    zonal de θ en la cara U, en degC·m³/s (`gad_advection.F:841`;
    available_diagnostics.log). Positivo hacia +x, lejos del hielo. Se
    convierte a W con H = ρ0 c_p Σ[−ADVx_TH + T_ref·ū·dA], con
    T_ref = −1.90 °C. Es la media exacta de ρ0 c_p Σ(−u)(T − T_ref)dA con la T
    de cara del propio esquema, incluido u'T' de escala menor a 3 h.
  - **Q_ex:** instantáneas horarias de UVEL (|u| no es lineal; las medias de
    3 h lo sesgan hacia abajo).
  - **Comprobaciones cruzadas:** los productos de medias de 3 h de u y T quedan
    en `cross_check` de `fig3_metrics.json`, con su diferencia relativa.
- **Disco, por caso de 30 d:**

  | Lista | Cálculo | Tamaño |
  |---|---|---|
  | `snap1h` (UVEL + THETA horarios) | 720 × 2 × 2.4 MB | 3.5 GB |
  | `mean3h` (UVEL, THETA, ADVx_TH) | 240 × 3 × 2.4 MB | 1.7 GB |
  | `snap1d` | — | 0.36 GB |
  | `animV` (VVEL en el nivel 19, 3 h) | 240 × 48 kB | 12 MB |
  | `ETAN` y pickups | — | ≈ 0.2 GB |
  | **Total** | | **≈ 5.8 GB por caso, ≈ 23 GB en total** |

  El buffer de diagnósticos usa 502 niveles de los 1000 de `numDiags`
  (`DIAGNOSTICS_SIZE.h`).
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

## Corrida cronometrada de 1 día (2026-11-05): falla de la frontera de Orlanski

`c_fjord_all_short1d` (2 etapas de 12 h, np = 5, dt = 60 s) terminó OK y es muy barata: 13 s + 15 s de pared. **Pero el dominio se vacía:**

| t (h) | η medio (m) | obc_E_uVel_Int (m³/s) |
|---|---|---|
| 0 | 0 | 0 |
| 12 | −3.03 | +28 900 (sale) |
| 24 | −11.38 | +58 300 (sale) |

- **Lo esperado:** la fuente neta es Q0 = 100 m³/s (Σ addMass/ρ0 = Q0, `net_addmass_over_rho0_Q0` = 1.0). η medio debería *subir* ~Q0·t/A ≈ 0.036 m/día si nada saliera, y quedar ≈ 0 si la frontera sacara Q0.
- **Lo observado:** sale 300–600 veces Q0. Es la misma inestabilidad de Orlanski que en (b), aquí en sentido de vaciado.
- **Por qué Flather tampoco sirve en (c):** c = √(9.81·500) ≈ 70 m/s, así que c·dt/dx = 70·60/200 ≈ 21. Forma explícita inestable.
- **Propuesta:** Stevens, igual que en (b), con u_ob media en profundidad = Q0_total/A_ob. El forzamiento de plataforma (OBEt/OBEs variables + esponja T/S) se mantiene. **Costo:** el modo barotrópico queda cerrado en la frontera (no hay intercambio de setup de viento con la plataforma). El intercambio baroclínico (Jackson 2018), que es el objetivo de (c), sigue libre.

## Stevens + zona RBCS: criterios de la prueba de 1 día, fijados antes de correr (2026-11-05)

**Configuración:**
- Stevens en el este, con u_ob = Q0/A_ob.
- La esponja T/S de OBCS se reemplaza por una zona RBCS equivalente (R33; commits 78a191d, 5a6d27e, b5a0588).
- La esponja U/V se mantiene.

**Caso:** `c_fjord_all_stevens_short1d`. Son 2 etapas de 12 h, con `--set wind_stage1_s=43200` y np = 5.

**Criterios (todos deben pasar):**
1. Termina OK en las dos etapas, y la etapa 2 arranca desde `pickup_stevensE`.
2. `obc_E_uVel_Int` = +100 m³/s, con error relativo < 1e-6, en todas las salidas del monitor.
3. |η medio del dominio| < 0.01 m en todo momento. Lo esperado es ≈ 0: la frontera saca exactamente Q0 = Σ addMass/ρ0, y ni pkg/icefront ni RBCS agregan volumen.
4. T y S quedan dentro del rango físico de las condiciones iniciales y de la plataforma, con un margen de 0.5 °C y 0.2 psu:
   - mín(T) ≥ mín(T_init, T_shelf) − 0.5; máx(T) ≤ máx(T_init, T_shelf, T_out) + 0.5;
   - lo mismo para S;
   - sin NaN.
5. **Informativo:** se reporta el desplazamiento de la picnoclina en la zona RBCS a las 24 h frente a la serie de plataforma (debe seguir el forzamiento de ±30 m).

**Si falla cualquiera de 1–4,** no se lanza la campaña de (c) y se analiza el fallo.

## Resultado de la prueba de 1 día con Stevens + zona RBCS (2026-11-05, 22:35)

**Corrida:** `/Volumes/LaCie/portfolio_runs/exp1/c_fjord_all_stevens_short1d`. El tiempo de pared fue 13 s + 14 s con np = 5.

| # | Criterio | Resultado | Veredicto |
|---|---|---|---|
| 1 | Las 2 etapas terminan OK; la etapa 2 arranca desde `pickup_stevensE` | `run_status_final` OK (720 + 720 pasos). STDOUT de la etapa 2: `MDS_READ_SEC_YZ: opening file: pickup_stevensE.0000000720.*` | PASA |
| 2 | `obc_E_uVel_Int` = +100 m³/s, con error relativo < 1e-6 | Error relativo máximo de 6.9e-13 (26 salidas del monitor) | PASA |
| 3 | \|η medio\| < 0.01 m | \|η medio\| máximo de 3.6e-16 m | PASA |
| 4 | T y S dentro del rango físico ± el margen, sin NaN | T ∈ [−1.49983, 4.00000] °C y S ∈ [32.90005, 34.60000] psu. Las entradas tienen T ∈ [−1.49997, 4.0] y S ∈ [32.9000, 34.6]. No hay NaN | PASA |
| 5 | Informativo: la picnoclina de la zona RBCS frente a la de la plataforma | Ver la tabla siguiente | **No sigue el forzamiento** |

**Criterio 5.** La picnoclina se mide como la profundidad del valor medio de T de cada columna, interpolada linealmente. Se promedia en y sobre las columnas mojadas.

| | 0 h | 12 h | 24 h |
|---|---|---|---|
| Objetivo de plataforma (OBEt = rbcs_T), analítico 200 − 30 sin(2πt/6 d) | 200.0 | 185.0 | 174.0 |
| Zona, isl 1–3 | 200.0 | 198.5 | 192.9 |
| Zona, isl 4–10 | 200.0 | 198.8 | 193.3 |
| Al oeste de la zona, isl 11–20 | 200.0 | 198.8 | 193.2 |
| Centro del fiordo, x = nx/2 | 200.0 | 197.7 | 194.2 |

**Lectura.**
- La zona sube solo ~1.5 m a las 12 h y ~7 m a las 24 h, mientras que el objetivo sube 15 m y 26 m. La zona queda prácticamente igual que el interior. Incluso isl = 1, con una tasa de 0.2/600 s⁻¹ (e-folding ≈ 50 min), queda en 192.5 m.
- El código RBCS es el estándar: `rbcs_add_tendency.F:88-97`, con la máscara × maskC de `rbcs_init_fixed.F:240`.
  - Los registros temporales se cargan bien: a las 12 h se leen los registros 12 y 13 con pesos 0.5/0.5, según el log `RBCS_FIELDS_LOAD`.
  - rbcs_T coincide con OBEt.
- **Hipótesis (sin verificar):** el levantamiento que impone la zona se exporta como onda interna, hacia el fiordo y hacia afuera por la frontera radiante de Stevens, más rápido de lo que la relajación lo repone.
  - La zona mide 2 km. La relajación es fuerte solo en las 2–3 columnas exteriores: e-folding ≥ 2.7 h desde isl = 3.
  - El interior sí se levanta (~6–7 m a las 24 h, en todo el fiordo), así que la señal entra, pero atenuada a ~25 % del forzamiento nominal.
- Un día no basta para separar la amplitud del desfase, porque el periodo es de 6 d. Las 4 corridas de 30 d (5 ciclos) lo miden directamente.
- **Clase (c), elección numérica:** la intensidad de la zona heredada de los parámetros `sponge_*` de cases.yaml.
- **Consecuencia para el análisis:** el forzamiento de plataforma efectivo dentro del fiordo no es ±30 m. Cualquier resultado de (c) debe reportarse frente a la amplitud de picnoclina **realizada** en el fiordo, medida en la corrida, y no frente a la nominal. Esto vale igual que con los parámetros reales de case.json.
- Cambiar la intensidad de la zona después de ver este resultado sería ajustar el diseño a un criterio informativo. Solo se haría con una nueva prerregistración y con aprobación de Esteban (PD1(c)).

## Campaña de 4 × 30 d (opción A, 2026-11-05): dos casos divergen

`run_campaign` reporta OK en los 4 casos, pero **`c_fjord_shelf` da NaN desde el día 12.8 y `c_fjord_all` desde el día 20.3** (etapa 2). Los casos `c_fjord_plume` y `c_fjord_wind` terminan finitos.
- Antes del NaN, max|u| oscila entre 0.4 y 1 m/s; en `plume` es ~0.3 m/s.
- Es la misma firma que en (b): velocidades espurias junto a la OB de Stevens (ver pd1_b_twoplume3d.md, "Hallazgo crítico"). Aquí se dispara con el forzamiento de plataforma.
- **Las corridas de (c) no se usan** hasta corregir la frontera. `require_ok` debe rechazar campos no finitos.

## Verificación cruzada de la OB de (c) con pyOCEAN: criterios fijados antes de correr (2026-11-06, aprobado por Esteban)

**Pregunta.** La divergencia de (c) con forzamiento de plataforma, ¿es propia de la OB de MITgcm (Stevens, la esponja autorreferente, el balance que ignora `addMass`) o del diseño?

**Configuración.**
- El modelo es pyOCEAN, hidrostático, en modo ALE z*.
- La OB este usa Flather implícito (forma θ), Orlanski y una esponja de 10 celdas (1200/60 s) con trazadores relajados hacia un exterior variable en el tiempo.
- La malla y la estratificación son las de (c): 200 × 30 celdas de 200 m, 50 niveles de 10 m, f0 = 1.37e-4, perfil Sermilik y EOS lineal con tAlpha/sBeta de (c).
- Viscosidad y difusividad: Ah = Kh = 30 y Ar = Kr = 1e-5.
- El exterior es la serie de plataforma de (c), con la picnoclina en 200 − 30 sin(2πt/6 d).
- La fuente es una descarga de Q = 100 m³/s como runoff en la superficie de la cabeza del fiordo. Es un análogo estuarino, porque pyOCEAN no tiene fuente en profundidad: diferencia de clase (b) frente a la pluma de (c), declarada.
- **Casos:**
  - `pyo_shelf`: solo plataforma, Q = 0. Es el análogo de `c_fjord_shelf`, que en MITgcm diverge en el día 12.8.
  - `pyo_shelf_runoff`: plataforma más Q = 100.
- Duración: 30 días.
- El código vive en `exp1_fjord_plumes/crosscheck/pyocean_c/`, solo importa pyOCEAN y no modifica `~/pyOCEAN`.

**Criterios (todos deben pasar):**
1. Llega a los 30 d con campos finitos.
2. Balance de volumen. Con Q = 0, |η̄| < 0.01 m en todo momento. Con Q = 100, el transporte neto saliente por la OB, promediado en los días 20–30, vale Q ± 2 %, y la deriva de η̄ en los días 20–30 es < 0.01 m.
3. En cada salida diaria, max|u| en las 3 columnas junto a la OB ≤ 1.5 × max|u| en el interior (columnas 10–189). Es el mismo criterio que en MITgcm.
4. **Informativo:** la amplitud de la picnoclina realizada a mitad del fiordo (ajuste senoidal de periodo 6 d en los días 6–30) frente a los 30 m nominales, y su fase.

**Interpretación pre-registrada.**
- Si pyOCEAN pasa 1–3, la falla de (c) se atribuye a la OB de MITgcm, como clase (c), elección numérica.
- Si falla, se atribuye al diseño del forzamiento o de la zona, y se reconsidera el diseño de (c).

## Verificación de V3 (Stevens + u baroclínica interior) en `c_fjord_shelf` a 30 d (2026-11-06)

**Corrida:** `/Volumes/LaCie/portfolio_runs/exp1/c_fjord_shelf_stv_iv`. Antes, con Stevens por defecto, este caso divergió en el día 12.8.

**Criterios pre-registrados (pd1_b_twoplume3d.md, "Diagnóstico de la OB"): todos PASAN.**
1. `require_ok` pasa: 30 d, campos finitos.
2. En el volcado de U (solo existe el del día 30), la razón OB/interior es 1.20.
3. |η̄| máx 4.3e-16 m; transporte de la OB, error máximo 3e-10 m³/s.

**Chequeo adicional, no pre-registrado; se reporta tal cual.** Se aplicó el mismo criterio 2 a los 720 snapshots horarios de UVEL (`snap1h`):
- 1 de 720 supera 1.5. Es el de la hora 1, con razón 1.54, cuando |u| ≈ 1e-3 m/s durante el arranque.
- Después del día 1, el máximo es 1.44 y la mediana 1.11.
- El máximo de |u| junto a la OB es 0.23 m/s, frente a 1–2 m/s y la divergencia que había antes.

La frontera sigue con velocidades 10–40 % mayores que el interior, pero sin crecimiento. Esto queda como limitación conocida (clase (c)).

**Decisión (R38):** se adopta V3 en las configuraciones de producción (b) y (c) y se vuelven a correr los 14 + 4 casos, como en el plan aprobado. Las corridas anteriores se conservan renombradas con el sufijo `_stv0`, como evidencia de la falla.

**Precisiones fijadas antes de evaluar la verificación con pyOCEAN** (commit 32a21ed, antes de cualquier corrida de 30 d):
- "En todo momento", para Q = 0, significa cada paso del modelo.
- "Deriva", para Q = 100, es el máximo de |η̄(t) − η̄(20 d)| en los días 20–30, la lectura más estricta.
- La picnoclina es la isoterma de 1.25 °C, que reproduce los 200.0 m de la tabla de MITgcm en t = 0.
- La "OB ± 3 columnas" son las caras u 198–200, porque en pyOCEAN la OB está en el borde del dominio.
- El runoff de 100 m³/s se reparte en las 6 celdas de la pluma de (c), en j = 12–17.

**Diferencias declaradas (clase b/c).** El runoff entra en superficie, no en profundidad. La esponja es de 1200/60 s, no la zona de (c) de 600/21600 s. No hay ajuste convectivo. Los esquemas de advección son distintos. `pyo_shelf` es el análogo limpio de `c_fjord_shelf`.

## Caso de control sin forzamiento `c_fjord_control` (aprobado por Esteban, 6 nov; fijado antes de correr)

**Configuración.** Es (c) con Q0 = 0, sin viento y con shelf_amp = 0: el exterior OBEt/OBEs y rbcs_T/S son el perfil inicial constante, lo que verifica `test_c_control_has_no_forcing`. Usa la misma OB (V3) y 30 d.

**Para qué sirve.** Un fiordo estratificado en reposo, con EOS lineal y fondo plano, debería quedarse en reposo. Cualquier movimiento es espurio: viene de la OB, de la zona RBCS o de errores numéricos. Ese movimiento fija el **piso de ruido** del intercambio y de las velocidades de los otros cuatro casos.

**Criterios (deben pasar):**
1. `require_ok` pasa: 30 d con campos finitos.
2. |η̄| < 0.01 m y |`obc_E_uVel_Int`| < 1e-6 m³/s en todas las salidas del monitor.

**Expectativa (informativa; si no se cumple se reporta y se investiga):**
3. El máximo de |u| en el dominio es < 0.01 m/s durante los 30 d (serie `dynstat_uvel_max`).

**Uso en el análisis.** El transporte de intercambio de cada caso de (c) se reporta junto al del control, calculado con la misma métrica y en la misma sección. Una señal menor que 3× la del control no se interpreta como física.

### Re-corrida de los 4 casos (c) con V3 (2026-11-06)

Todos pasan los criterios 1–3. Razón OB/interior en los volcados de U:

| caso | razón |
|---|---|
| plume | **1.49** |
| wind (etapa 1) | 0.93 |
| wind (etapa 2) | 1.39 |
| shelf | 1.20 |
| all (etapa 1) | 1.00 |
| all (etapa 2) | 1.35 |

- |η̄| ≤ 1.6e-14 m y el error del transporte de la OB es ≤ 3.5e-10 m³/s.
- **`plume` queda justo en el límite (1.49 frente a 1.5).** Pasa, pero la frontera de (c) acelera el flujo de salida un 35–50 % sobre el interior en varios casos. Es una limitación que se declara (clase c). El control sin forzamiento dirá cuánto de esto es espurio.

## Resultado de la verificación cruzada con pyOCEAN (2026-11-06)

Corridas: `/Volumes/LaCie/portfolio_runs/exp1/crosscheck_pyo/{pyo_shelf,pyo_shelf_runoff}`, 30 d cada una (~68 min de reloj, 5 hilos). Evaluador: `crosscheck/pyocean_c/eval_pyo_c.py` → `figures/crosscheck_pyo_c.json`.

**Veredictos pre-registrados, tal como salieron:**

| Criterio | `pyo_shelf` (Q = 0) | `pyo_shelf_runoff` (Q = 100) |
|---|---|---|
| 1. 30 d, campos finitos | PASA | PASA |
| 2. Balance de volumen | PASA: max\|η̄\| = 8.5e-3 m | **FALLA**: Q_out(d20–30) = 97.48 m³/s (−2.5 %, límite ±2 %); deriva de η̄ d20–30 = 1.05e-2 m (límite 0.01) |
| 3. OB/interior ≤ 1.5 | PASA: máx 0.51 (día 12) | PASA: máx 0.32 (día 7) |
| 4. Informativo: picnoclina a mitad del fiordo | A = 29.6 m (0.99 de 30), desfase 27° (0.45 d) | A = 29.7 m (0.99), desfase 27° |

**Veredicto global pre-registrado: FALLA** (por el criterio 2 de `pyo_shelf_runoff`, en ambas partes y por poco).

**Diagnóstico post-hoc (no cambia el veredicto):**
- η̄ y el caudal saliente por la OB oscilan con el periodo de 6 d del forzamiento de plataforma (η̄ entre 0.005 y 0.019 m; Q_out horario entre 84 y 124 m³/s).
- La ventana d20–30 cubre 1.67 periodos, no un número entero; su promedio y su "deriva" (η̄(30) − η̄(20)) dependen de la fase de la oscilación. El criterio se diseñó pensando en un estado estacionario y no previó esto.
- Volumen acumulado en 30 d (5 periodos completos): salida 2.551e8 m³ frente a runoff 2.592e8 m³, es decir 98.4 % (−1.6 %). La diferencia (4.1e6 m³) es la que queda guardada en η̄(30 d). No hay fuga ni crecimiento.
- En ninguno de los dos casos aparece inestabilidad en la OB: con el mismo forzamiento de plataforma que hizo divergir a MITgcm con Stevens por defecto (día 12.8 y 20.3), pyOCEAN corre 30 d con la OB por debajo del interior (razón ≤ 0.51) y transmite al fiordo el 99 % de la amplitud de la picnoclina.

**Interpretación.** La regla pre-registrada dice que si pyOCEAN falla, la causa se atribuye al diseño. La falla ocurrió en el balance de volumen, no en la estabilidad, que es lo que motivó la verificación. Cómo leer esto es decisión de Esteban; no se redefine el criterio.

**Decisión de Esteban (2026-11-06): opción (b).** La FALLA se reporta tal cual, con el diagnóstico post-hoc al lado. No se usa como evidencia sobre el diseño de (c) ni se redefine el criterio. Lo que la verificación sí aporta (estabilidad de 30 d de la OB de pyOCEAN con el mismo forzamiento) se cita como evidencia de apoyo, de clase (c), a que la divergencia original era propia de la OB Stevens por defecto de MITgcm.
