# PD0 — Preregistro del Exp. 3 (aprobado por Esteban el 2026-11-10)

Fecha: 2026-11-10. **Ninguna corrida del protocolo** (`run_protocol.py --cv` o `--final`) se ha hecho. No existe ninguna predicción de un emulador para las 3 corridas a Tbot = 2.0 ni para `ocean0_4km`. Lo que sí se vio antes de fijar los umbrales está en "Divulgación". Esteban aprobó el PD0 con las decisiones de abajo. La configuración quedó congelada con `run_protocol.py --freeze-config` en `figures/protocol_config.json`: el protocolo la lee, y nada se cambia después.

**Fuentes**
- Spec `2026-11-10-exp3-melt-emulator-design.md` y nota `rev_literatura.md` (§1–8).
- Burgard et al. 2022 (B22), Ec. 14, 32–35. Favier et al. 2019 (F19), Ec. 2, 4, 5. Jourdain et al. 2020 (J20), Ec. 1.
- Myren & Lawrence 2021 (ML21), §3.1. Lakshminarayanan et al. 2017 (L17).
- Ockenden et al. 2026 (O26), *GMD* 19:7979: dispersión = desviación estándar de los miembros del ensamble, p. 7991 (Fig. 8b); umbral sugerido "> 10 %" del melt integrado, en las conclusiones.
- Lambert & Burgard 2025 (LB25), Ec. 1–3. El PD0 del Exp. 2 define las ventanas y N1.

## Definiciones comunes

- **Datos:** `/Volumes/LaCie/portfolio_runs/exp2/exp3_dataset.npz`, sha256 `3d2eceb1…6200f`. Se verifica antes de cada uso. Hay 29 corridas: 25 de entrenamiento (15 del factorial y 10 de geometría), la referencia `ocean0_4km` y 3 de prueba (Tbot = 2.0; z_c = 250, 400, 550).
  - El spec dice "26 corridas"; eso quedó superado cuando se decidió dejar `ocean0_4km` fuera de todo entrenamiento.
- **Melt:** m a⁻¹ de agua dulce, positivo = fusión, la media de los meses 19–24. Los valores negativos (rehielo, hasta −1.36 m a⁻¹) no se recortan.
- **Plataforma:** el `shelf_mask` propio de cada corrida. Tiene 726, 696 y 620 celdas con s = 0.8, 1.0 y 1.2. Cada método predice en la malla unión, pero cada métrica usa solo la máscara de la corrida evaluada.
- **Ponderación por área:** la malla es uniforme (4 km), así que la ponderación es la media simple sobre la máscara.
- **TF** = T − T_f(S, p(z_draft)), con el liquidus de AD16 Tabla 4 (constantes importadas del Exp. 2). Solo M7 usa TF₊ = max(TF, 0); M1–M3 tienen signo (R1).
- **Región GL:** la fracción **0.10** de las celdas de plataforma más cercanas a la línea de apoyo (`dist_gl`). El número de celdas es `round(0.10·n)`, y los empates se rompen por índice.
- **N1:** las 5 corridas que fallan N1 (`T−0.5_zc250`, `T−0.5_zc550`, `T1.5_zc550`, `T−0.42_zc402_s0.8`, `T−0.39_zc503_s1.2`) entran al entrenamiento y se señalan en cada tabla por corrida.

## Pliegues (`folds.py`, Task 2)

Todo el ajuste se hace dentro de cada pliegue: la estandarización, las EOF, el núcleo, K y n, el grado de M4 y la elección de `tf_kind`.

| Tipo | Pliegues | Entrenamiento / prueba |
|---|---|---|
| `loro` | 25 | 24 / 1, una vez por cada corrida de entrenamiento |
| `zc_out` | 3 | tercios por z_c con cortes en 325 y 475 m: z_c < 325 (17/8), 325 ≤ z_c ≤ 475 (16/9), z_c > 475 (17/8). Incluyen las corridas de geometría |
| `warm_out` | 1 | se entrena con Tbot ≤ 1.0 (20) y se prueba con Tbot ≥ 1.4 (5: tres a 1.5, una a 1.42 y una a 1.43; `T1.5_zc550` falla N1) |
| `profile_test` | 1 | 25 / `ocean0_4km`. **Solo para M0 y los métodos por celda** (M1–M3, M6, M7). M0 entra porque la destreza se mide contra él. M4 y M5 no lo ven, porque solo reciben (Tbot, z_c, s) y no pueden representar d = 720 m |
| `final` | 1 | 25 / las 3 a Tbot = 2.0. **Se evalúa una sola vez**, con `--final`, que se niega a correr si existe `figures/final_test.json` |

**Repetición de sensibilidad N1** (`drop_n1`): `--cv` repite todos los pliegues salvo `final` sin las 5 corridas que fallan N1 en el entrenamiento (la prueba no cambia). Se etiqueta "sensibilidad". La prueba final no se repite (Decisión 6). **H3a–H3c usan solo los pliegues principales**, nunca los de `drop_n1`.

## Métodos e hiperparámetros congelados (valores por defecto del código)

| # | Método | Configuración |
|---|---|---|
| M0 | Mapa medio | media por celda de los mapas de entrenamiento; una celda no vista recibe la media global |
| M1 | γ c TF | F19 Ec. 2. Un coeficiente, mínimos cuadrados en forma cerrada sobre todas las celdas de plataforma |
| M2 | K c² \|TF\| TF | B22 Ec. 14 (la misma forma funcional que F19 Ec. 4, salvo la constante, que K absorbe). Igual que M1 |
| M3 | K c² TF \|⟨TF⟩\| | J20 Ec. 1 (F19 Ec. 5). Igual que M1 |
| M4 | Regresión por celda | sobre (Tbot, z_c, s) estandarizados; grado 1 o 2 (10 términos), elegido por LORO interno con la media de `rmse_map`; un empate va al lineal |
| M5 | GP sobre EOF | ML21 §3.1. Estandarización por celda, SVD y el menor k con ≥ 99 % de la varianza, con tope k ≤ m − 1. Núcleo: Constant(1; 1e-3–1e3) × RBF ARD(1; 1e-2–1e3) + White(1e-3; **1e-8–1e-2**; un *nugget* pequeño y acotado, Decisión 4), `normalize_y`, 5 reinicios, `random_state` 0. σ se propaga por las EOF; el error de truncamiento no se suma |
| M6 | Ensamble de redes | entradas: las 10 de `FEATURES`. MLP 32-64-32 ReLU con salidas (μ, σ²), NLL (L17 Ec. 1), sin *bagging*, mezcla gaussiana uniforme. M = 5 semillas, 0–4 |
| M7 | Híbrido | m = K·b(TF)ⁿ·g_θ(x_g), con x_g = (draft, draft − z_c, dist_gl, dist_front, slope, s): sin ninguna temperatura. K y n se entrenan por miembro (n₀ = 2; K₀ sale de M2/M3 ajustado a los runs de entrenamiento del miembro). b = TF₊ ("local") o (TF₊\|⟨TF⟩\|)^½ ("nonlocal", R3). La elección es automática, por **LORO interno completo** (`tf_kind_folds = None`, Decisión 3) con la media de `rmse_map`; un empate va a local. La varianza sigue R11 (Decisión 2) |

**Comunes a M6 y M7:**
- Adam con lr = 1e-3, lotes de 512 celdas, float64 en CPU, `torch.use_deterministic_algorithms(True)` y 1 hilo.
- `max_epochs = 5000` como tope que no debe activarse (R13; se registra `hit_max_epochs`) y `patience = 250`.
- **Parada temprana:** cada miembro reserva 3 corridas de validación, estratificadas por Tbot (una al azar en el nivel más alto, una en la mediana y una en el más bajo, con la semilla del miembro; R10, R12). El criterio es la media de `rmse_map` sobre esas 3. Las corridas de validación no entran en los gradientes, ni en el estandarizador, ni en `y_scale`, ni en K₀.

## Métricas (Task 3)

- `rmse_map`: B22 Ec. 32.
- `integrated_rel_error` = (Σpred − Σmelt)/Σmelt: B22 Ec. 33–34, para una plataforma.
- `rmse_gl`: B22 Ec. 35. Por corrida es \|media_GL(pred) − media_GL(melt)\| (R8); la RMS sobre corridas se toma al agregar.
- `skill` = 1 − err/err_M0.
- Para M5, M6 y M7: cobertura al 68 % y al 95 %, y CRPS.
- En la prueba final, como descripción: dM/dT, la no linealidad y la amplificación profunda de LB25 Ec. 1–3.
- Todo se reporta por corrida, con las corridas N1 señaladas y con intervalos *bootstrap* sobre las corridas.

**Piso de ruido** (`analysis/noise_floor.py`, `figures/noise_floor.json`; npz sha256 `2146fc72…3cd035`):
- Para cada una de las 25 corridas: la RMS ponderada por área, sobre su `shelf_mask`, de (mapa de los meses 19–24 − mapa de los meses 13–18).
- **Mediana = 0.157 m a⁻¹** (rango 0.029–0.366). Para el melt integrado, la mediana de \|deriva N1\| = **1.3 %**.
- **Regla:** no se interpreta ninguna diferencia entre métodos menor que 0.157 m a⁻¹ en `rmse_map` o `rmse_gl`, ni menor que 1.3 % en `integrated_rel_error`.
- Es una cota conservadora: la diferencia entre ventanas incluye la deriva residual, y vale √2 veces el ruido de una sola media. El piso es ~4 % de los errores vistos (~4 m a⁻¹), así que rara vez será el límite.

## Hipótesis (umbrales aprobados por Esteban, 2026-11-10)

- **H3a (pregunta central).** En la prueba final, la media de \|`integrated_rel_error`\| de M7 sobre las 3 corridas es menor que la de M6, **y** M7 es mejor en ≥ 2 de las 3 corridas.
  - Si la diferencia de medias es < 1.3 % (el piso), el resultado es "no concluyente", no PASS.
  - Con 3 corridas, el resultado es un indicio, no una conclusión.
- **H3b.** En `loro`, la media de `rmse_map` de M7 sobre los 25 pliegues principales (no los de `drop_n1`) es ≤ **1.10** × la de M6.
- **H3c (réplica de O26).**
  - **Casos:** las corridas de prueba de `zc_out` (25), `warm_out` (5) y `final` (3), es decir, 33 casos, para M6 y para M7 por separado.
  - **Dispersión** = desviación estándar (ddof = 0) del melt integrado de los 5 miembros / \|melt integrado de la media del ensamble\|.
  - **Pasa** si ningún caso tiene dispersión < **5 %** y, a la vez, \|`integrated_rel_error`\| > **20 %**.
  - El umbral de 10 % de O26 se reporta como tabla de confusión (dispersión > 10 % contra error > 20 %).
- **El orden de los 8 métodos se reporta, no es hipótesis.** Si gana la cuadrática simple (M2), o M3, M4 o M5, **ese es el resultado principal del README**.
- Toda diferencia se clasifica como (a) error nuestro, (b) física faltante, (c) elección del método o (d) error de la referencia. Lo post hoc se etiqueta.

## Divulgación: lo que se vio antes de fijar los umbrales

Todo lo siguiente sale de las corridas de entrenamiento. Ninguna predicción toca Tbot = 2.0 ni `ocean0_4km`.

1. **M5, LORO exploratorio** (lo calculó el revisor de la Task 5 al auditar el GP), sobre las 25 corridas de entrenamiento: `rmse_map` medio **3.94** m a⁻¹, máximo ~**16** m a⁻¹. Es el único número de validación cruzada externa visto antes del PD0.
2. **M5 dentro de la muestra** (el *smoke* de la Task 5, en las 25 corridas): **k = 9** EOF con el **99.2 %** de la varianza; RMSE de reconstrucción 0.97 m a⁻¹; RMSE del GP dentro de la muestra 1.18 de media y 2.24 de máximo.
3. **Ajuste de la parada de M6/M7** (Task 6): `max_epochs`/`patience` pasó de 300/30 a 1000/250 después de mirar las curvas de validación interna en las corridas de entrenamiento del pliegue final. Luego R13 subió el tope a 5000: 0 de 10 miembros lo tocaron, y el más largo paró en la época 2470.
4. **RMSE dentro de la muestra**, con la media del ensamble sobre las 25 corridas: M6 **4.03**, M7 local **3.95**, M2 **5.23** m a⁻¹.
5. **Selección de `tf_kind`** (5 grupos internos, en las 25 corridas), en orden:
   - regla de parada de una sola corrida: ganó nonlocal (4.95 contra 5.72);
   - regla de 3 corridas: ganó local (4.76 contra 5.00);
   - regla actual (R12 + R13): gana **local, 4.51 contra nonlocal 5.18**.

   La elección **se invirtió** al cambiar la regla de parada.

   Las reglas R10–R12 (validación con 3 corridas estratificadas por Tbot, elegidas por miembro) salieron de lo que el revisor de la Task 6 observó en las corridas de entrenamiento:
   - n bimodal;
   - la cola de TF alto subajustada;
   - `ens_T1.5_zc250_s1.0` siempre quedaba fuera como validación.
6. **Exponente ajustado** de M7 local por miembro: n = **2.23–2.47** (mediana 2.39); K entre 16 y 34.
7. **Piso de ruido** (esta tarea): 0.157 m a⁻¹ y 1.3 %.
8. **Del Exp. 2** se conocen el H1 (exponente de 2.4–2.8) y el **melt medio de MITgcm de las 3 corridas a Tbot = 2.0**, que aparece en la Fig. 2a del README del Exp. 2 (círculos abiertos). Es el objetivo de la prueba final, pero ningún método se ha ajustado ni evaluado con él.
9. **Orden de las decisiones:** el *nugget* del GP y el LORO completo para `tf_kind` (Decisiones 3 y 4) se decidieron **después** de ver los puntos 1–5.

## Decisiones de Esteban (2026-11-10)

1. **Umbrales de H3a–H3c: aprobados tal como se propusieron.** H3a: ≥ 2 de 3 corridas, y "no concluyente" por debajo del piso de 1.3 %. H3b: 1.10. H3c: 5 % y 20 %.
2. **Varianza de M7 (R11): aprobada.** σ² = softplus(h₀(x_g))² + (softplus(h₁(x_g))·stopgrad(μ))² + 1e-6, heterocedástica, y la red no ve la temperatura. Se consideró la alternativa de una varianza solo geométrica.
3. **`tf_kind` por LORO interno completo** (`tf_kind_folds = None`, el nuevo valor por defecto de `Hybrid`). La razón es que la elección con 5 grupos fue frágil: se invirtió entre reglas de parada. La alternativa considerada eran 5 grupos (~9 min por ajuste).
4. ***Nugget* del GP de M5: WhiteKernel con cotas [1e-8, 1e-2] y valor inicial 1e-3**, en unidades normalizadas. Reemplaza la cota superior de 1e1; el revisor había propuesto ≤ 1e0.
   - ML21 Ec. 6–7 no tiene *nugget*.
   - Berdahl et al. 2021 (§2.3) usan "a small nugget for numerical conditioning", con la cota inferior de Ranjan et al. 2011.
   - Nuestras corridas tienen una deriva pequeña entre ventanas: el piso de ruido de 0.157 m a⁻¹ equivale a ≈ 1e-3 en unidades normalizadas de los coeficientes. Por eso el *nugget* es pequeño y acotado, y no puede absorber señal.
5. **(Informativa, aceptada.)** El docstring de M2 y el README no atribuyen c² a F19 Ec. 4: F19 usa c a la primera potencia, y la forma con c² es la de B22 Ec. 14. Como K absorbe la constante, el resultado no cambia.
6. **La repetición N1 excluye `final`: aceptado.** La alternativa considerada era repetirla dentro de la única llamada a `--final`.
7. **Regla de parada: 8 h de reloj, con ≤ 5 procesos concurrentes de 1 hilo.**
   - El plan (Task 8 Step 5) decía 6 h. Esteban subió el límite antes de cualquier corrida del protocolo, porque el LORO completo en `tf_kind` pone la estimación en ≈ 6 h o más.
   - La alternativa considerada, leer el límite como tiempo de CPU en serie, se descartó.
8. **Los pisos de 1.3 % (`integrated_rel_error`) y de 0.157 m a⁻¹ (`rmse_map`, `rmse_gl`) valen para todas las comparaciones entre métodos: aceptado.** La alternativa considerada era aplicar el piso solo a los mapas.
9. **H3c sobre 33 casos, por separado para M6 y para M7: aceptado.** La alternativa considerada era evaluarla solo para M6.

## Cómputo estimado de `--cv` y regla de parada

- **Ajustes por método:**
  - 29 pliegues (25 `loro` + 3 `zc_out` + 1 `warm_out`), más `profile_test` para M0 y los métodos por celda: 30;
  - ×2 por la repetición N1: **60 ajustes** (58 para M4 y M5).
- **Tiempos medidos** (Task 6, con hasta 5 procesos concurrentes, 1 hilo cada uno):
  - un ajuste de M7 con `tf_kind` fijo tarda ~50 s (549 s para el automático con 5 grupos, que son 11 ajustes);
  - M6 tarda ~1.1 min; M4 y M5, segundos; M0–M3, un tiempo despreciable.
- **M7 con LORO completo:** con m corridas de entrenamiento son 2m + 1 ajustes, ≈ (2m + 1) × 50 s.
  - `loro` (m = 24): ≈ 41 min;
  - `zc_out` (16–17): ≈ 28–29 min;
  - `warm_out` (20): ≈ 34 min;
  - `profile_test` (25): ≈ 43 min.
  - Pliegues principales ≈ 19.8 h de CPU; repetición N1 (m = 12–20) ≈ 15.9 h.
- **Total:**
  - M7 ≈ 35.7 h de CPU; M6 ≈ 1.1 h;
  - **≈ 37 h de CPU, ≈ 7.4 h de reloj** con 5 procesos. Está por debajo del límite, con poco margen.
- **Regla de parada:** el reloj se mide después de la primera tanda de pliegues. Si la estimación supera **8 h**, se detiene y se reporta, sin cambiar nada del protocolo para que quepa.
- **`--final`:** 8 ajustes, una sola vez; M7 ≈ 43 min.
