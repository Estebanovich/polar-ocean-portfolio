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
- **Con el forzamiento a priori, que es la prueba independiente, la hipótesis no se sostiene dentro de cada z_c.** A z_c fijo, n queda aún más lejos de [1.3, 2.0] que con Tbot: 2.4, 3.4 y 6.5. El melt crece de 4 a 15 veces mientras la media de TF_eff en la base cambia poco, sobre todo con z_c = 550, donde TF_eff va de 0.18 a 0.24 °C. La sensibilidad supercuadrática no se explica con la temperatura media que la rampa pone en la base: es una **respuesta de la cavidad**. El melt responde a la parte profunda y cálida de la base, cerca de la línea de apoyo, y a la circulación que esa agua impulsa. La media de área de la base la diluye.
- **El conjunto a priori da n = 1.55**, dentro del rango y cerca del 1.61 de H08. El colapso, sin embargo, es pobre (R² = 0.85 en log): no hay una sola ley de potencia en TF_eff para todos los z_c. La variante profunda colapsa mejor (n = 1.89, R² = 0.94), lo que concuerda con que domine la base profunda. Esta pendiente del conjunto mezcla la variación entre z_c con la que ocurre dentro de cada z_c, así que no se toma como evidencia a favor de la hipótesis.
- **Contra TD_diag hay un colapso casi perfecto** (n = 1.45, R² = 0.9995), pero es en parte circular. Lo que sí informa es esto: si m = c·γT·TD con γT = ΓT u*, entonces n − 1 ≈ 0.45 dice que u* (la velocidad en la capa límite) crece de forma sublineal con T_w − T_b. H08 §4 (`literatura/Holland_2008*.md` L650–657) encuentra en cambio que |u| y T − T_b crecen *ambos linealmente* con T_D, y por eso su producto es cuadrático. Aquí, la diferencia T_w − T_b efectiva crece más rápido que linealmente con Tbot, con TD_diag de 0.17 a 0.56 °C en z_c = 250, es decir, ×3.2 mientras Tbot − T_f,GL va ×2.1.
- **Referencia de H08 §6** (`literatura/Holland_2008*.md` L917–929): ajustando una ley de potencia a su modelo obtienen "an exponent of 1.61 (r² = 0.994)". En la figura se dibuja solo como guía, no como ajuste. Su T_D es la temperatura del agua profunda, uniforme desde 200 m (L501); no es la misma variable que nuestra Tbot de rampa ni que TF_eff.
- **Conclusión:** el exceso de n sobre 2 no se debe solo a haber elegido Tbot como variable. Con el forzamiento a priori en la base, el exponente dentro de cada z_c sigue fuera de [1.3, 2.0]. La sensibilidad supercuadrática pertenece a la respuesta de la cavidad (dónde está el agua cálida respecto a la base profunda, y cómo responde la circulación), no a la variable de forzamiento. Clasificación (c)/(b) [verificar]: puede ser un efecto del diseño de rampa, con el gradiente vertical concentrado cerca de la línea de apoyo, o de la física de la cavidad a 4 km. **H1 sigue siendo FALLA, tal como se preregistró.**
