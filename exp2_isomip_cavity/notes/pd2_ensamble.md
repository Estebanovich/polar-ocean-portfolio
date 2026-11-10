# PD2 — Ensamble a 4 km: H1, H2 y N1–N3

Fecha: 2026-11-09. Cola de producción (PID 23111): el ensamble terminó a las 22:03 CST (`END ens rc=0`). Ocean1 sigue corriendo y no entra aquí.

Definiciones: `pd0_preregistro.md`, incluida la aclaración R22, que se fijó a las 16:41 CST (commit 968410f), antes de evaluar cualquier salida del ensamble. Todo se calcula sobre los últimos 6 meses (registros 19–24) de corridas de 2 años a 4 km, con ΓT_cal,4km = 0.029356.

Código y salidas:
- código: `analysis/ensemble.py`, `fig2_melt_sensitivity.py` y `fig3_meltwater_export.py`;
- tabla: `figures/ensemble_table.json`;
- veredictos: `figures/fig2_h1.json` y `figures/fig3_h2.json`.

## N1–N3 (29 corridas: referencia y 28 `ens_*`)

- **N3:** pasan las 29; no hay huecos.
- **N2:** pasan las 29 (máx. |residuo| = 3.1e-7 < 1 %).
- **N1:** fallan 5, con deriva > 3 %:

  | Corrida | N1 |
  |---|---|
  | `ens_T-0.5_zc250_s1.0` | 3.7 % |
  | `ens_T-0.5_zc550_s1.0` | 6.5 % |
  | `ens_T1.5_zc550_s1.0` | 3.4 % |
  | `ens_T-0.42_zc402_s0.8` | 8.8 % |
  | `ens_T-0.39_zc503_s1.2` | 4.9 % |

  Son, sobre todo, casos fríos o de termoclina profunda, con melt pequeño (0.3–4.5 m a⁻¹), en los que la misma deriva absoluta pesa más en términos relativos. Según la aclaración R22, siguen en el ajuste y se señalan. Ninguna corrida se alargó.
- **Verificación cruzada:** `ensemble.collect` reproduce exactamente los N1 y N2 de `ensemble_n_check.py`, el paso `ncheck` de la cola (diferencia máxima 0.0).

## H1 — forma cuadrática (criterio: ⟨m⟩ de toda la plataforma, R22)

T_f,GL = −2.3341 °C (Sbot = 34.7; base mojada más profunda 565.2 m).

| z_c (m) | R² cuadrática | n (ley de potencia) | k | Veredicto | Corridas con N1 fallido |
|---|---|---|---|---|---|
| 250 | 0.9998 ✔ | 2.43 ✘ | 1.04 | **FALLA** | T −0.5 |
| 400 | 0.9995 ✔ | 2.78 ✘ | 0.29 | **FALLA** | — |
| 550 | 1.0000 ✔ | 2.70 ✘ | 0.069 | **FALLA** | T −0.5, T 1.5 |

**H1 falla en los tres z_c.** La forma cuadrática ajusta muy bien, pero el exponente queda por encima del rango preregistrado [1.3, 2.0]: el melt crece más rápido que (Tbot − T_f,GL)².

**Secundario, no criterio (R22):** con el melt profundo (base < −300 m) sale n = 1.99 (z_c 250, dentro del rango), 2.49 (400) y 2.62 (550), con R² ≥ 0.9998. No cambia el veredicto.

**Robustez (*post hoc*):** se calcularon tres variantes, todas sobre el melt de toda la plataforma.

| z_c (m) | n sin un punto (5 ajustes) | n sin las corridas con N1 fallido | pendientes locales d ln m / d ln ΔT |
|---|---|---|---|
| 250 | 2.40–2.48 | 2.40 | 2.19–2.55 |
| 400 | 2.75–2.82 | 2.78 | 2.68–2.99 |
| 550 | 2.60–2.76 | 2.67 | 2.43–2.88 |

En ninguna variante n baja de 2.19. El fallo no depende de un punto aislado ni de las corridas que fallan N1.

**Lectura (hipótesis, sin verificar):** en la familia "rampa", la base del hielo poco profunda está en la termoclina o por encima de ella. Al subir Tbot, crecen a la vez la temperatura del agua que llega a la base y el área de base expuesta a agua cálida. Una variable explicativa que es Tbot, y no la temperatura en la base, puede dar un exponente efectivo mayor que 2. Clasificación provisional: (c), por elección del diseño (la variable Tbot de la rampa), no (a). Para probarlo haría falta un ajuste contra la temperatura media en la base, que sería *post hoc*.

## H2 — profundidad de salida del agua de deshielo

| z_c (m) | D_out (m) para Tbot = −0.5, 0.0, 0.5, 1.0, 1.5 | Spearman ρ | Veredicto |
|---|---|---|---|
| 250 | 250.6, 199.3, 141.8, 127.9, 117.6 | −1.0 | **PASA** |
| 400 | 374.2, 322.4, 181.2, 133.4, 102.6 | −1.0 | **PASA** |
| 550 | 460.4, 334.5, 113.7, 97.4, 90.7 | −1.0 | **PASA** |

**H2 pasa en los tres z_c:** D_out baja de forma monótona con Tbot. Las corridas de prueba (Tbot = 2.0: 113, 99 y 97 m) siguen la tendencia, pero no entran en H2. Ocean0: D_out = 119 m.

## Figuras

- `figures/fig2_melt_sensitivity.{png,pdf}`, con tres paneles:
  - (a) ⟨m⟩ contra Tbot, con el ajuste cuadrático, el punto de Ocean0 (10.19 m a⁻¹) y Tbot = 2.0 como prueba;
  - (b) mapa de calor;
  - (c) efecto de s, solo descriptivo.
- `figures/fig3_meltwater_export.{png,pdf}`: (a) secciones en x_front para Tbot −0.5 y 1.5 con z_c = 400; (b) D_out contra Tbot; (c) mapa de ∫C dz para Tbot 1.5 y z_c 400.
- `figures/fig3_anim_ocean0_4km_r21.gif`: 0.71 MB, 24 meses.

**Resolución.** El pie de la Figura 2 remite a la prueba de resolución, que es descriptiva. Ahí, el melt de toda la plataforma de Ocean0 difiere en −0.05 % entre 4 y 2 km. Esa coincidencia es, en gran parte, el objetivo de la calibración (ΓT se calibró por separado a cada resolución, K1), no evidencia independiente. Las métricas que no se calibraron difieren más: el máx. de |ψ_bt| en la cavidad, en +125 %, y los mapas de melt correlacionan con r = 0.91. El ensamble es a 4 km y sus números de melt heredan esa limitación.

ψ_bt no se usa en ningún resultado: a 4 km su máximo en la cavidad es 2.2 veces el de 2 km.

## H1 post hoc: forzamiento en la base del hielo (R23, aprobado 2026-11-10)

**POST HOC.** Este análisis no es una prueba preregistrada y no cambia el veredicto de H1, que sigue siendo **FALLA** en los tres z_c, como se preregistró. Código: `analysis/h1_posthoc_tf.py`. Salidas: `figures/fig2_posthoc_tf.{png,pdf,json}`.

**Hipótesis.** El exponente n > 2 aparece porque Tbot, la temperatura de fondo lejos del hielo, no es la temperatura que ve el hielo. En la familia "rampa", el hielo poco profundo queda en la termoclina o por encima de ella.

**Forzamientos** (15 corridas del núcleo; ⟨m⟩ = melt medio de toda la plataforma, últimos 6 meses):
- **TF_eff, a priori, solo con entradas:** media ponderada por área, sobre las columnas de plataforma, de T_rampa(z_d) − T_f(S_rampa(z_d), ρ_ref g |z_d|) (AD16 Ec. 25). z_d es la base del hielo de `icetopo.bin`; T y S salen de `profiles.ramp` con Tbot, z_c y d de case.json. Hay una variante para toda la plataforma y otra para la región profunda (< −300 m).
- **TD_diag, diagnosticado:** es la estimación de `jpl_compare.py`, T_w − T_b = −⟨SHIfwFlx⟩ L / (ρ_ref c_w ⟨SHIgammT⟩) (AD16 Ec. 24, κ_i = 0), con la misma advertencia: es una media ponderada por γT. Además **es en parte circular**, porque sale del propio melt. Como m = c·γT·TD_diag, su exponente vale 1 + d ln⟨γT⟩ / d ln TD_diag y mide la respuesta de u*; no es una relación forzamiento–melt independiente.
- **Ningún TF ≤ 0;** no se excluyó ninguna corrida.

**Ajustes ⟨m⟩ = k·TFⁿ** (MCO en log–log; R² en espacio logarítmico; entre corchetes, el rango de n al dejar fuera una corrida):

| Forzamiento | z_c = 250 | z_c = 400 | z_c = 550 | Conjunto (15 corridas) |
|---|---|---|---|---|
| Tbot − T_f,GL (variable de pd0, referencia) | 2.43 [2.40–2.48] | 2.78 [2.75–2.82] | 2.70 [2.60–2.76] | 2.64, R² = 0.32 |
| TF_eff, toda la plataforma (a priori) | 2.36 [2.34–2.40] | 3.44 [3.42–3.45] | 6.49 [5.88–6.94] | **1.55**, R² = 0.85 [1.41–1.61] |
| TF_eff, < −300 m (a priori) | 2.35 [2.33–2.39] | 2.89 [2.86–2.92] | 4.26 [3.95–4.48] | 1.89, R² = 0.94 [1.78–1.96] |
| TD_diag, toda la plataforma (en parte circular) | 1.52 [1.49–1.55] | 1.44 [1.42–1.45] | 1.41 [1.39–1.41] | 1.45, R² = 0.9995 [1.44–1.45] |
| TD_diag, < −300 m (en parte circular) | 1.99 [1.98–2.01] | 1.76 [1.66–1.96] | 1.49 [1.42–1.67] | 1.74, R² = 0.987 [1.71–1.84] |

**Lectura honesta.**
- **Con el forzamiento a priori, que es la prueba independiente, la hipótesis no se sostiene dentro de cada z_c.** A z_c fijo, n queda aún más lejos de [1.3, 2.0] que con Tbot: 2.4, 3.4 y 6.5. El melt crece 5.9, 7.8 y 7.3 veces (z_c = 250, 400, 550) mientras la media de TF_eff en la base cambia poco, sobre todo con z_c = 550, donde TF_eff va de 0.18 a 0.24 °C. La sensibilidad supercuadrática no se explica con la temperatura media que la rampa pone en la base: es una **respuesta de la cavidad**. El melt responde a la parte profunda y cálida de la base, cerca de la línea de apoyo, y a la circulación que esa agua impulsa. La media de área de la base la diluye.
- **El conjunto a priori da n = 1.55**, dentro del rango y cerca del 1.61 de H08. El colapso, sin embargo, es pobre (R² = 0.85 en log): no hay una sola ley de potencia en TF_eff para todos los z_c. La variante profunda colapsa mejor (n = 1.89, R² = 0.94), lo que concuerda con que domine la base profunda. Esta pendiente del conjunto mezcla la variación entre z_c con la que ocurre dentro de cada z_c, así que no se toma como evidencia a favor de la hipótesis.
- **Contra TD_diag hay un colapso casi perfecto** (n = 1.45, R² = 0.9995), pero es en parte circular. Lo que sí informa es esto: si m = c·γT·TD con γT = ΓT u*, entonces n − 1 ≈ 0.45 dice que u* (la velocidad en la capa límite) crece de forma sublineal con T_w − T_b. H08 §4 (`literatura/Holland_2008*.md` L650–657) encuentra en cambio que |u| y T − T_b crecen *ambos linealmente* con T_D, y por eso su producto es cuadrático. Aquí, la diferencia T_w − T_b efectiva crece más rápido que linealmente con Tbot, con TD_diag de 0.17 a 0.56 °C en z_c = 250, es decir, ×3.2 mientras Tbot − T_f,GL va ×2.1.
- **Referencia de H08 §6** (`literatura/Holland_2008*.md` L917–929): ajustando una ley de potencia a su modelo obtienen "an exponent of 1.61 (r² = 0.994)". En la figura se dibuja solo como guía, no como ajuste. Su T_D es la temperatura del agua profunda, uniforme desde 200 m (L501); no es la misma variable que nuestra Tbot de rampa ni que TF_eff.
- **Conclusión:** el exceso de n sobre 2 no se debe solo a haber elegido Tbot como variable. Con el forzamiento a priori en la base, el exponente dentro de cada z_c sigue fuera de [1.3, 2.0]. La sensibilidad supercuadrática pertenece a la respuesta de la cavidad (dónde está el agua cálida respecto a la base profunda, y cómo responde la circulación), no a la variable de forzamiento. Clasificación (c)/(b) [verificar]: puede ser un efecto del diseño de rampa, con el gradiente vertical concentrado cerca de la línea de apoyo, o de la física de la cavidad a 4 km. **H1 sigue siendo FALLA, tal como se preregistró.**

## R24 — prueba post hoc de exposición (perfil tipo Holland 2008), predicción escrita antes de correr

**POST HOC** (aprobado por Esteban, 2026-11-10). No cambia el veredicto preregistrado de H1, que sigue siendo **FALLA**. Esta sección se escribe y se sube **antes** de lanzar las corridas; la hora del commit es la evidencia.

**Contexto.** La banda de H1, n ∈ [1.3, 2.0], viene de Holland et al. 2008. Su agua profunda es uniforme desde 200 m (`literatura/Holland_2008*.md` L501, §3a), así que el agua cálida cubre la base profunda del hielo. En nuestro núcleo, la termoclina de la rampa (z_c = 250, 400 y 550 m; d = 300 m) deja cálida solo una fracción pequeña de la base. Con R23 vimos que el forzamiento medio en la base no explica n > 2 dentro de cada z_c.

**Corridas.**
- Casos `hol_T{−0.5, 0.0, 0.5, 1.0, 1.5}` en `cases.yaml`: rampa con z_c = 150 m y d = 100 m. Es decir, −1.9 °C y S0 por encima de 100 m, rampa lineal de 100 a 200 m, y Tbot y Sbot por debajo de 200 m. Se usa `profiles.ramp` sin cambios.
- Todo lo demás es igual al núcleo del ensamble: 4 km, ΓT_cal,4km = 0.029356, 2 años, `melt_law`, `checkIniTemp`, s = 1.
- **No** son miembros del ensamble: no tienen `ensemble_group`.
- **Verificación en seco:** `make_case` de `hol_T1.0` contra `ens_T1.0_zc250_s1.0`.
  - Los binarios de geometría (batimetría, `icetopo.bin`, máscara rbcs) son idénticos.
  - Los namelists solo difieren en la línea de comentario con el nombre del caso.
  - En case.json solo difieren `case`, `z_c`, `d` y `ensemble_group`.
  - Difieren, como debe ser, T/S iniciales y de restauración y `shelficeLoadAnomaly.bin`.
  - Perfil comprobado: −1.90 °C a 90 m, −0.45 °C a 150 m, 1.00 °C a 210 m.

**Exposición de la base** (fracción del área de plataforma, s = 1, geometría a 4 km):

| Perfil | Base en agua a Tbot plena | Base en agua fría de superficie | Fracción de rampa r media (área) |
|---|---|---|---|
| `hol` (z_c 150, d 100) | 0.497 | 0.062 | 0.824 |
| núcleo z_c = 250 | 0.126 | 0.062 | 0.435 |
| núcleo z_c = 400 | 0.011 | 0.711 | 0.132 |
| núcleo z_c = 550 | 0.000 | 0.875 | 0.032 |

**Advertencia sobre la premisa.** El **50.3 %** del área de la base está a menos de 200 m de profundidad (6.2 % a menos de 100 m; la base llega hasta −5.6 m). El perfil tipo Holland **no** pone toda la base en agua cálida, y ningún perfil de rampa podría hacerlo, salvo uno uniforme. La prueba es, por tanto, de **exposición mucho mayor**: ~50 % de la base en agua a Tbot plena, contra ≤ 13 % en el núcleo. No es de exposición total. Esto se reporta tal cual y la regla de lectura no se cambia.

**Hipótesis.** Si n > 2 se debe a la exposición parcial de la base a la termoclina, con la base mucho más expuesta el exponente contra (Tbot − T_f,GL) debería caer en [1.3, 2.0].

**Regla de lectura (fijada ahora).**
- Mismos ajustes que H1, con `ensemble.quadratic_fit` y `ensemble.power_law_fit` sobre las 5 corridas `hol_T*`.
- ⟨m⟩ = melt medio de toda la plataforma (R22), registros 19–24.
- T_f,GL = `metrics.freezing_point_gl`, la misma geometría que el núcleo: −2.3341 °C.
- Se reporta n con su rango al dejar fuera una corrida, y R² de la cuadrática.
- **"Apoya"** si n ∈ [1.3, 2.0]; **"no apoya"** en otro caso.
- Se reportan N1–N3 de las 5 corridas. Un fallo de N3 deja la prueba "no evaluable": no se imputa ni se vuelve a correr. Un fallo de N1 o N2 se señala y la corrida se queda en el ajuste (como R22).

### R24, parte 2 (sin corridas): predictores cuadráticos local y no local (POST HOC)

Código: `h1_posthoc_tf.py --quad`. Salidas: `figures/fig2_posthoc_quad.{png,pdf,json}`. TF_loc es el forzamiento lejano a priori en la base de cada columna (el mismo de R23). ⟨·⟩ es la media ponderada por área sobre la plataforma. Se ajusta m = k·P^q (MCO en log–log); la familia cuadrática predice q = 1.

**Predictores:**
- **P_nonlocal = ⟨TF_loc·⟨TF⟩⟩.** Favier et al. 2019, Ec. 5 (`literatura/Favier_2019*.md` L491–497): "quadratic, local and non-local". Es la "quadratic-semilocal" de Burgard et al. 2022, Ec. 18 (`literatura/Burgard_2022*.md` L687–702). **Identidad:** ⟨TF⟩ es un solo número por corrida, así que P_nonlocal = ⟨TF⟩² = TF_eff² exactamente y q = n(TF_eff)/2. Este predictor **no aporta información nueva** sobre R23; se reporta porque se pidió.
- **P_local = ⟨TF_loc·|TF_loc|⟩.** Favier 2019, Ec. 4 (L447–455, "quadratic, local"), en la forma que conserva el signo de Burgard 2022, Ec. 14 (L648–654, "quadratic-local"). No se reduce a TF_eff.

**Resultados.** No hay ningún P ≤ 0, así que no se excluyó ninguna corrida. Entre corchetes, el rango de q al dejar fuera una corrida.

| Predictor | z_c = 250 | z_c = 400 | z_c = 550 | Conjunto (15 corridas) |
|---|---|---|---|---|
| P_nonlocal (= TF_eff²) | 1.18 [1.17–1.20] | 1.72 [1.71–1.73] | 3.24 [2.94–3.47] | 0.78, R² = 0.85 [0.71–0.80] |
| P_local | 1.18 [1.17–1.20] | 1.44 [1.43–1.46] | 1.83 [1.70–1.92] | **1.02**, R² = 0.96 [0.97–1.05] |

**Lectura.**
- En el conjunto de las 15 corridas, el cuadrático local da q = 1.02 (R² en log = 0.96, rango LOO 0.97–1.05), con mejor colapso que TF_eff (R² = 0.85). **Esa pendiente no se toma como evidencia:** igual que en R23, mezcla la variación entre z_c con la que ocurre dentro de cada z_c, y dentro de cada z_c q vale 1.18–1.83. Que q del conjunto quede cerca de 1 no muestra que la familia cuadrática describa la respuesta a Tbot. Además, las corridas `hol` de R24 no caen sobre esa misma recta (`fig2_posthoc_r24`, panel b). Pesar cada columna por su propio forzamiento, que es lo que hace ⟨TF_loc|TF_loc|⟩, da más peso a la base profunda y cálida. Esto concuerda con la lectura de R23: el melt lo controla la parte expuesta de la base, no la media.
- **Dentro de cada z_c, q sigue por encima de 1** (1.18, 1.44 y 1.83) y crece con z_c, es decir, con menos exposición. Ningún predictor a priori de esta familia reproduce la sensibilidad a Tbot a z_c fijo. El residuo concuerda con una respuesta de la circulación de la cavidad que no está en TF_loc.
- **H1 sigue siendo FALLA**, como se preregistró.

## R24 — lectura de la prueba de exposición (POST HOC)

Lectura mecánica de la regla fijada antes de correr (commit 874205e). Cola `queue_r24.sh`: 11:01–11:34 CST, rc = 0; corridas de unos 33 min cada una. Código: `h1_posthoc_tf.py --r24`. Salidas: `figures/fig2_posthoc_r24.{png,pdf,json}`.

**N1–N3 de las 5 corridas `hol_T*`:**
- N3: pasan las 5, sin huecos.
- N1: máximo 0.51 % (`hol_T-0.5`); pasan las 5.
- N2: máximo 2.0e-7; pasan las 5.

**Ajustes de H1 sobre `hol`** (⟨m⟩ de toda la plataforma, registros 19–24; T_f,GL = −2.3341 °C). ⟨m⟩ = 8.17, 14.12, 21.46, 30.21 y 40.31 m a⁻¹ para Tbot = −0.5 … 1.5.
- Cuadrática: R² = 0.9999995.
- Ley de potencia contra (Tbot − T_f,GL): **n = 2.16**, rango LOO 2.11–2.19 (R² en log = 0.9996).

**Veredicto (regla): NO APOYA.** n = 2.16 queda fuera de [1.3, 2.0], también en todo el rango LOO.

**Contexto (no es regla):**

| Variable | hol | núcleo z_c = 250 | 400 | 550 |
|---|---|---|---|---|
| n contra Tbot − T_f,GL | **2.16** [2.11–2.19] | 2.43 | 2.78 | 2.70 |
| Base en agua a Tbot plena (fracción de área) | 0.497 | 0.126 | 0.011 | 0.000 |
| n contra TF_eff (plataforma, a priori) | 1.98 [1.96–1.99] | 2.36 | 3.44 | 6.49 |
| q contra P_local = ⟨TF_loc·\|TF_loc\|⟩ | 0.99 [0.98–1.00] | 1.18 | 1.44 | 1.83 |

**Lectura honesta.**
- La regla dice "no apoya". Aun así, con mucha más exposición el exponente baja (2.16 frente a 2.43–2.78) y queda cerca del borde de la banda. El orden no es estrictamente monótono en la exposición: z_c = 400 (2.78) está por encima de z_c = 550 (2.70).
- **Advertencia de la premisa:** solo ~50 % del área de la base está por debajo de 200 m (el 50.3 % está por encima). No es la exposición total del supuesto de la hipótesis. Una exposición completa podría acercar n a la banda, pero eso no se probó y la geometría no lo permite con un perfil de rampa.
- Dentro de la familia `hol`, el melt sigue casi exactamente las formas de forzamiento local: n = 1.98 contra TF_eff y q = 0.99 contra el cuadrático local. Con la base expuesta, la sensibilidad a z_c fijo coincide con la familia cuadrática local. En el núcleo, con menos exposición, no coincide (q = 1.18–1.83). Esto concuerda con que el exceso de n en el núcleo venga de la exposición parcial a la termoclina. Pero es una lectura **post hoc** y la regla preregistrada para R24 dio "no apoya".
- **H1 sigue siendo FALLA,** como se preregistró.
