# PD2 (fallos R20, R21): convergencia y sensibilidad de (a), plan prerregistrado

**Estado: preparado (casos, configuraciones, builds); NADA corrido.** La campaña de (a) sigue en curso. Revisión 1 (ronda de correcciones): criterios mecánicos, tope de cg3d de 2000 (R21), r_eff y presupuesto condicionado.

## Motivo

En la evaluación provisional (twolayer q0 = 0.01 y 0.03), C2 y C2'' pasan, pero **C2' falla**:

- la velocidad media de la pluma en su ancho, Q/b, es 0.52–0.56 de la teoría;
- la velocidad del núcleo w es 0.88–0.95 de la teoría;
- el ancho es de 20–24 m contra 7.5 m de la teoría, unas 1.5 celdas de dx0 = 5 m.

Atribución propuesta: resolución en la cara más la pared sin deslizamiento (clase c), posiblemente con entrainment efectivo extra (clase b). Se **prueba** esta atribución; la métrica no se redefine.

## Casos (todos `twolayer`, q0 = 0.03 m² s⁻¹, 24 h de modelo, np = 1)

| caso | config / build | dx0 (m) | dz (m) | nx × nz | deltaT (s) | L (m) | cambia respecto a producción |
|---|---|---|---|---|---|---|---|
| `conv_dx10` | `conv_dx10` | 10 | 10 | 170 × 50 | 5.0 | 20 063 | malla más gruesa |
| `conv_dx5` | `conv_dx5` | 5 | 5 | 193 × 100 | 2.5 | 20 140 | solo dz (producción tiene dz = 10) |
| `conv_dx2p5` | `conv_dx2p5` | 2.5 | 2.5 | 216 × 200 | 1.25 | 20 133 | malla más fina |
| `sens_freeslip` | `sens_line2d` | 5 | 10 | 193 × 50 | 3.0 | 20 140 | `no_slip_sides = no_slip_bottom = .FALSE.` |
| `sens_xu2012`: "cierre al estilo Xu (A4 reescalada a la malla)" | `sens_line2d` | 5 | 10 | 193 × 50 | 3.0 | 20 140 | cierre de mezcla (abajo) |
| referencia: `a_line2d_twolayer_q0.03` | `a_line2d` | 5 | 10 | 193 × 50 | 3.0 | 20 140 | producción (48 h; se reevalúa en la misma ventana de 12–24 h) |

- **Lo que se mantiene igual:** perfiles, fuente, EOS, icefront, OBCS, cg3dTargetResidual = 1e-8 y salidas.
- **Diferencias comunes a todos los estudios:**
  - `cg3dMaxIters = 2000` (fallo R21; producción tenía 500);
  - `monitorFreq = deltaT`, es decir, monitor en cada paso. Solo afecta a la salida, para poder hacer la verificación por paso descrita abajo.
- **El campo lejano queda fijo** (r = 1.03, dx_max = 200 m). Por eso el refinamiento no es uniforme: dz se refina exactamente ×2 en cada paso, pero en x el factor depende de la distancia a la cara. Las longitudes L difieren menos de 0.4 % (la suma de dx alcanza L con un exceso menor que la última celda). Se considera despreciable a 20 km de la pluma.

| razón entre mallas | dx medio en 0–25 m de la cara (centros) | dx local a 20 m | dx local a 100 m |
|---|---|---|---|
| dx10 → dx5 | 10.150 / 5.309 = **1.912** | 1.885 | 1.626 |
| dx5 → dx2p5 | 5.309 / 2.822 = **1.881** | 1.777 | 1.445 |

r_eff se define como la razón del dx medio sobre la región de la pluma (centros de celda entre 0 y 25 m de la cara). Se usa en Richardson (r21 = 1.881, r32 = 1.912; ver abajo).

- **dt por CFL:** deltaT = 0.5 s m⁻¹ × min(dx0, dz). La producción dio CFL_w = 0.198 con w_max = 0.66 m/s y dt/dz = 0.3 s/m. Con dt/min(Δ) = 0.5 s/m, CFL ≈ 0.33 si w no cambia, y ≤ 0.45 si el núcleo llega a ~0.9 m/s.
- **Escalado comprobado por prueba** (`tests/test_pd2_cases.py`):
  - `icefrontlength` = 1/dx0 (0.1, 0.2 y 0.4 m⁻¹) y `R_icefront` = H, solo en la columna 0; por tanto `K_icefront` = Nr y el factor parcial vale 1 (`icefront_thermodynamics.F:158-174`);
  - fuente ρ0·q0·dy = 308.1 kg s⁻¹ en la celda del fondo junto a la cara;
  - `SIZE.h` coincide con la malla;
  - `case.json` registra `git_commit` y `git_dirty`.

### Cierre al estilo Xu (A4 reescalada a la malla), caso `sens_xu2012`

Xu et al. (2012), *Ann. Glaciol.* 53(60), 229–234, **p. 231**: νv = 0.1 m² s⁻¹, biarmónica horizontal 300 m⁴ s⁻¹ (con dx = 20 m) y κ = ν (Pr = 1).

- Se pone `viscAr = diffKrT = diffKrS = 0.1` y `viscAh = diffKh = 0`.
- La biarmónica literal es inestable aquí: el límite es L⁴/(32 dt) = 5⁴/96 = 6.5 m⁴ s⁻¹ (`pkg/mom_common/mom_calc_visc.F:48-49`). Se reescala con el mismo tiempo de amortiguamiento de malla dx⁴/A4: A4 = 300·(5/20)⁴ = **1.171875 m⁴ s⁻¹** en `viscA4`, `diffK4T` y `diffK4S`. Para la biarmónica hacen falta OLx, OLy ≥ 3 (`config_check.F:566-575`); se cumple.
- **No es un cambio "solo de cierre".**
  - Con `viscAh = 0`, el arrastre sin deslizamiento sobre w en la cara casi desaparece. El arrastre lateral de w, en `MOM_W_SIDEDRAG` (llamado desde `calc_gw.F:456-465`), es −sideDragFactor·(viscAh·w − viscA4·∇²w)·… (`mom_w_sidedrag.F:73-80`), y con viscAh = 0 solo queda el término biarmónico, pequeño.
  - Además, viscAr sube 100× (1e-3 → 0.1), lo que aumenta mucho el arrastre del fondo sin deslizamiento y la mezcla vertical.
  - Por eso un cambio de R_Q en este caso **no se atribuye al cierre**: se solapa con `sens_freeslip`. Se reporta solo como sensibilidad conjunta "cierre + arrastre de pared".

## Costo (modelo revisado, R21)

Modelo por paso, ajustado a las corridas cortas medidas: con 40 iteraciones de cg3d, 9.44 ms/paso; con ~290, 24.4 ms/paso. Queda

  t_paso = (celdas / 9650) × (7.04 ms + 0.060 ms × N_it)

Comprobación: N_it ≈ 340 da 27.4 ms contra 27.1 ms medidos (dt = 3 s). N_it = número medio de iteraciones de cg3d por paso. En producción es 320–390; con el tope de 2000, una corrida válida tiene N_it ≤ 1999. No se midió el costo extra de monitorear cada paso; se espera < 5 %.

| caso | celdas | pasos (24 h) | N_it = 350 | N_it = 700 | N_it = 1999 (máximo válido) |
|---|---|---|---|---|---|
| `conv_dx10` | 8 500 | 17 280 | 0.12 h | 0.21 h | 0.54 h |
| `conv_dx5` | 19 300 | 34 560 | 0.54 h | 0.94 h | 2.44 h |
| `conv_dx2p5` | 43 200 | 69 120 | 2.41 h | 4.22 h | **10.92 h** |
| `sens_freeslip` | 9 650 | 28 800 | 0.22 h | 0.39 h | 1.02 h |
| `sens_xu2012` | 9 650 | 28 800 | 0.22 h | 0.39 h | 1.02 h |
| **total** | | | **3.5 h** | **6.2 h** | **15.9 h** |

`conv_dx2p5` pasa de 8 h si N_it medio > 1434.

## Presupuesto: sonda y punto de parada

1. Primero, una **sonda de 1 h** de `conv_dx2p5` (2880 pasos; 6–27 min). De ahí salen N_it medio, CFL y la extrapolación a 24 h con el modelo de arriba.
2. Si la extrapolación es ≤ 8 h, se corre completo. **Si es > 8 h, el plan se detiene y se reporta a Esteban** (como en PD1). Partir la corrida en tramos con pickup (C3 verificado) solo se hace **si Esteban lo aprueba** tras ver la sonda.

## Criterios prerregistrados (mecánicos)

### Definiciones

Todas las métricas se calculan con el código de C2' sin cambios: `plume_profile` + `evaluate`, teoría en modo `mitgcm` (fallo R13), en `analysis/fig1_plume_validation.py`.

- **R_Q** = `u_ratio` = media de `w_mean` (= Q/b) en la mitad inferior de la pluma teórica, dividida entre la u de sombrero de copa de la teoría.
- **R_w** = media de `w_max` en la misma máscara, dividida entre la u de la teoría (solo informativa).
- **b** = media de `width` (m) en la misma máscara. Ojo: `width` está cuantizado en caras de celda.

### Ventana

**12–24 h para todas las corridas de estudio y para la referencia de producción** (reevaluada en 12–24 h, no en su ventana de 48 h). El valor de cada caso, X(12–24), sale de `plume_profile(run, (12, 24))`. Las instantáneas horarias usadas son 12, 13, …, 24 h; la de 24 h viene del volcado de estado, porque la instantánea de diagnóstico en endTime no existe.

### Estacionariedad

|R_Q(18–24) − R_Q(12–18)| / R_Q(12–24) < 0.02, con ventanas (12, 18) y (18, 24) (comparten la instantánea de 18 h). Si no se cumple, la corrida es **no estacionaria** y no entra en la clasificación.

### σ por caso

σ_X = desviación estándar (ddof = 1) de las 12 muestras X(t−1, t), t = 13, …, 24 h. Cada muestra es `plume_profile(run, (t−1, t))`: promedia las instantáneas t−1 y t. Las muestras se solapan, así que σ es una medida de variabilidad, no un error estándar. Es conservadora.

### Umbral de diferencia entre dos casos A y B

σ_AB = √(σ_A² + σ_B²). Una diferencia cuenta solo si:

- para R_Q: |R_A − R_B| > **max(2 σ_AB, 0.03)**;
- para b: |b_A − b_B| > **max(2 σ_AB, 2.5 m)** (el dx0 de la malla más fina).

Cada paso de refinamiento (dx10 → dx5, dx5 → dx2p5) se etiqueta **sube**, **baja** o **plano** según ese umbral.

### Validez de cada corrida (se comprueba antes de clasificar)

1. `run_status.json` OK (`require_ok`).
2. **cg3d (R21):** en `output.txt` hay exactamente nTimeSteps líneas `cg3d_iters (last)` (monitor en cada paso), y **todas son < 2000**. Si algún paso llega al tope, la corrida es **inválida**.
3. CFL (`advcfl_uvel_max`, `advcfl_wvel_max`) ≤ 0.5 en todos los pasos.
4. Estacionaria (arriba).

Si alguna corrida `conv_*` es inválida o no estacionaria, la convergencia es **inconclusa** y se dice cuál falló.

### Extrapolación de Richardson

Notación: φ1 = R(dx2p5) (fina), φ2 = R(dx5), φ3 = R(dx10); r21 = 1.881, r32 = 1.912 (r_eff de la tabla); ε21 = φ2 − φ1, ε32 = φ3 − φ2.

- Como r21 ≠ r32, se usa la forma de razón variable (Celik et al. 2008, Procedimiento GCI, Ec. 3–4) `[verificar]`:
  - s = signo(ε32/ε21);
  - p = |ln|ε32/ε21| + q(p)| / ln r21, con q(p) = ln[(r21^p − s)/(r32^p − s)], resuelto por iteración de punto fijo desde q = 0.
  - Si r21 = r32 = r, se reduce a la fórmula del revisor: p = ln[(R10 − R5)/(R5 − R2.5)] / ln r.
- R_ext = (r21^p φ1 − φ2)/(r21^p − 1).
- **La extrapolación es inconclusa, y no se da como número,** si ocurre cualquiera de:
  - s < 0 (secuencia no monótona);
  - ε21 = 0 o ε32 = 0;
  - la iteración no converge;
  - p < 0.5 o p > 3.
- Solo como diagnóstico, no decide nada: p calculado con r = 2 (el refinamiento de dz).

### Clases de resultado para la convergencia (conv_dx10, conv_dx5, conv_dx2p5)

Se evalúan en este orden y vale la primera que aplica.

| # | condición | resultado |
|---|---|---|
| V | alguna corrida inválida o no estacionaria | **Inconclusa (validez)** |
| G | b/dx0 constante: en ambos pasos, \|(b_gruesa/b_fina)/r_eff − 1\| < 0.2 | **Ancho fijado por la malla.** La pluma no está resuelta; la evidencia es **inconclusa para R_Q** sea cual sea su tendencia |
| S | R_Q sube en ambos pasos, b baja (más que el umbral) en ambos pasos, Richardson concluyente y R_ext ∈ [0.75, 1.25] | **Atribución clase (c) a la resolución: APOYADA** |
| P1 | R_Q sube en ambos pasos, b baja en ambos, Richardson concluyente y R_ext < 0.75 | **Apoyo parcial:** se cuantifica la parte numérica (R_ext − R_Q(dx5)); el resto queda sin explicar (candidato clase b) |
| P2 | R_Q sube en ambos pasos, b baja en ambos, Richardson inconcluso o R_ext > 1.25 | **Tendencia a favor sin número:** se reporta la tendencia; no hay valor extrapolado |
| W | R_Q sube en ambos pasos pero b no baja en ambos | **No apoyada:** R_Q cambia sin que la pluma se estreche; se evalúa entrainment o física 2D (clase b) |
| D | R_Q baja en ambos pasos | **Refuta clase (c);** apunta a un error nuestro (a) o a física faltante (b) |
| F | R_Q plano en ambos pasos | **Refuta clase (c):** el déficit no depende de la resolución; evaluar (b) |
| M | cualquier otra combinación (no monótona, o un paso plano y otro con tendencia) | **Inconclusa (no monótona):** sin atribución |

En ningún caso se cambia la métrica. La clase (d) (definición de Q/b frente al sombrero de copa de la teoría) solo se plantea con evidencia independiente.

### Resultados separados, que no entran en la tabla

- **Sensibilidad a dz:** producción (dx0 = 5, dz = 10) contra `conv_dx5` (dx0 = 5, dz = 5), ambas en 12–24 h. ΔR_Q y Δb con el mismo umbral; se reporta sube, baja o plano.
- **Pared:** `sens_freeslip` contra producción. Si R_Q sube más que el umbral, el arrastre sin deslizamiento contribuye (clase c), cuantificado como ΔR_Q. Si no, no hay contribución detectable.
- **Cierre al estilo Xu:** `sens_xu2012` contra producción. ΔR_Q se reporta como sensibilidad conjunta "cierre + arrastre de pared", sin atribuirlo al cierre. No se adopta la mezcla que "mejor acierta" (regla: no ajustar parámetros para que coincida un caso).
- En cada corrida se reporta además si C2 (despegue ±0.10 H) y C2'' (derretimiento, factor 1.5) siguen dentro de tolerancia.

## Cómo lanzar (después de la campaña de (a), uno a la vez)

```
cd ~/Phd_project/portfolio && export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
P=/opt/anaconda3/envs/pyocean/bin/python; S=exp1_fjord_plumes/scripts/run_campaign.py
$P $S --config conv_dx10 --cases conv_dx10 --np 1
$P $S --config conv_dx5 --cases conv_dx5 --np 1
$P $S --config sens_line2d --cases 'sens_*' --np 1
$P $S --config conv_dx2p5 --cases conv_dx2p5 --np 1 --end-time 3600 --run-suffix _probe1h   # sonda; luego PARAR si > 8 h
$P $S --config conv_dx2p5 --cases conv_dx2p5 --np 1                                         # solo si la sonda da <= 8 h
```

## Aclaraciones previas a los resultados (2026-11-05, antes de evaluar)

Al implementar el evaluador mecánico (`analysis/pd2_convergence_eval.py`, Task 5b) aparecieron 11 ambigüedades en los criterios de arriba. Se resuelven aquí (ruling R26) **antes de ver ningún resultado**: `conv_dx10` se usó solo como prueba de humo de los lectores (validez, CFL, dx medio en 0–25 m, estacionariedad); no se vio ni se calculó ninguna clasificación. Las demás corridas no se han evaluado. Los criterios de arriba no cambian; esto solo fija lo que el texto dejaba abierto. El código implementa exactamente lo que sigue.

1. **Monitor de la referencia de producción.** El criterio 2 de validez pide nTimeSteps líneas `cg3d_iters (last)`, pero la producción monitorea cada 200 pasos (`monitorFreq = 600`, 288 líneas). Resolución: para la referencia (y solo para ella) se comprueban las líneas que existen contra **su propio tope** (`cg3dMaxIters` leído del `data` de cada corrida, 500 en producción; 2000 si no se encuentra), con la bandera "muestreado cada 200 pasos, no verificado paso a paso" (`per_step_verified = False`, no invalida). Un tope alcanzado entre líneas muestreadas no es detectable. El CFL se comprueba igual (≤ 0.5) sobre las líneas existentes.
2. **Validez de `sens_*` y de la referencia.** Las mismas comprobaciones de validez y estacionariedad se aplican a **toda** corrida. Una `sens_*` inválida o una referencia inválida deja la sensibilidad correspondiente como "no disponible o inválida" (no entra en la tabla de convergencia).
3. **Corrida `conv_*` ausente** (por ejemplo, `conv_dx2p5` detenida por la regla de 8 h): sin clase en la tabla; el evaluador devuelve `INCOMPLETE` e indica qué falta.
4. **Ancho b.** La observación de que `width` está cuantizado en caras de celda es anterior a R22. Se usa el ancho de R22 (borde interpolado entre centros) con los mismos umbrales (piso 2.5 m).
5. **P1, "R_ext − R_Q(dx5)".** R_Q(dx5) es el de `conv_dx5`, que también tiene dz = 5 (producción tiene dz = 10); la "parte numérica" mezcla refinamiento en x y en z. Se calcula tal como dice el texto (con `conv_dx5`) y se reporta con esa salvedad.
6. **ε = 0 y la iteración de Richardson.** ε = 0 es cero exacto en coma flotante; un ε casi nulo da p muy grande o muy pequeño y cae en la regla p ∉ [0.5, 3]. Tolerancia del punto fijo 1e-10 sobre p, máximo 200 iteraciones; sin convergencia es inconcluso.
7. **Prueba G.** Usa las razones crudas b_gruesa/b_fina (sin el umbral de detección) y los r_eff de la tabla del texto (1.912 para dx10→dx5; 1.881 para dx5→dx2p5), no los calculados de la malla de cada corrida (el dx medio en 0–25 m de cada corrida se guarda en el JSON solo como verificación).
8. **Signo de "sube/baja/plano".** Diferencia = fina − gruesa en los pasos de convergencia, y sensibilidad − producción en las sensibilidades (dz: `conv_dx5` − producción). Desigualdad estricta en el umbral: |Δ| > umbral.
9. **Estacionariedad.** Un solo criterio, sobre R_Q (como está escrito); b no se prueba.
10. **Richardson "concluyente".** Es independiente de las etiquetas de los pasos (puede ser concluyente con un paso "plano"); solo importa dentro de S/P1/P2, que exigen "sube" en ambos pasos. R_ext = 0.75 y 1.25 pertenecen a S (intervalo cerrado).
11. **σ y umbrales.** σ es la desviación estándar (ddof = 1) de las 12 muestras X(t−1, t); en la prueba de humo en `conv_dx10` σ_RQ ≈ 1e-4 y σ_b ≈ 0.016 m (flujo casi estacionario), de modo que en la práctica deciden los pisos (0.03 para R_Q; 2.5 m para b). Se mantiene tal como se prerregistró.
12. **Divergencia monótona (añadida tras la revisión del código del evaluador, R27; sigue siendo antes de ver ningún dato de las corridas de estudio).** La fórmula del orden p usa valor absoluto y aceptaría como concluyente una secuencia monótona **divergente** (por ejemplo R_Q = 0.40, 0.50, 0.70 con incrementos crecientes al refinar daba p = 1.154 y la clase S). Siguiendo el estándar de Celik/Stern, si |ε32/ε21| ≤ 1 (los incrementos no se reducen al refinar) la extrapolación es **inconclusa** con la razón "divergencia monótona", y en `classify` cae en P2 si R_Q sube y b baja en ambos pasos. Esto se suma a las condiciones de inconclusión del texto original (s < 0, ε = 0, sin convergencia, p ∉ [0.5, 3]); el caso |ε32/ε21| = 1 exacto también es divergencia.
13. **Corridas `conv_*` inválidas y anchos nulos (añadida en la misma revisión).** Si alguna corrida `conv_*` es inválida (incluida una que falla `require_ok` y por tanto no tiene métricas), la clase es V sin calcular pasos, Richardson ni la prueba G. Un ancho b cero, negativo o no finito en alguna `conv_*` también da V (razón "zero or non-finite plume width"). El resultado de P1 lleva la salvedad escrita "mezcla refinamiento en x y en z" (ítem 5).


## Resultado (evaluador congelado, 2026-11-05)

**Clase F: refuta la clase (c).** El déficit de R_Q no depende de la resolución, así que lo que sigue es evaluar (b).

| corrida | dx0 cara (m) | R_Q (12–24 h) | b (m) | estacionaria |
|---|---|---|---|---|
| conv_dx10 | 10.15 | 0.538 | 21.0 | sí |
| conv_dx5 | 5.31 | 0.563 | 20.2 | sí |
| conv_dx2p5 | 2.82 | 0.573 | 20.1 | sí |

- **Incrementos de R_Q:** +0.025 y +0.010. Los dos están por debajo del umbral de 0.03, así que cuentan como "plano". b también es plano.
- **Richardson:** concluyente, con p = 1.43 y R_ext = 0.579. Muy por fuera de la tolerancia de C2' (0.75–1.25).
- **G (ancho fijado por la malla):** no aplica. Las razones (b_c/b_f)/r_eff son 0.54 y 0.54, y b ≈ 20 m está resuelto con 7 a 14 celdas.
- **Sensibilidades:**
  - dz 10 → 5 m: ΔR_Q = +0.016, plano;
  - free slip: ΔR_Q = +0.019, plano, sin contribución detectable del arrastre de pared;
  - `sens_xu2012`: inválida (no estacionaria, 2.6 % entre 12–18 h y 18–24 h). Se reporta solo como dato: R_Q ≈ 0.72, b ≈ 9.7 m.
- **C2 y C2'':** pasan en todas las corridas válidas. El flujo de volumen y la altura de desprendimiento están convergidos y concuerdan con la teoría. C2' (velocidad media) **sigue FALLIDO** y no se redefine.

Lectura (hipótesis, todavía sin probar): la pluma de MITgcm es ~2.6× más ancha y ~0.57× más lenta que el sombrero de copa teórico, con un flujo de volumen comparable. Con la mezcla de malla, la pluma 2D resuelta se ensancha más que lo que fija α = 0.1 en la teoría. Candidata a clase (b): mezcla turbulenta no representada, porque el cierre (viscAh/viscAr constantes) no es un LES. La pista que apunta en esa dirección es que `sens_xu2012`, con otro cierre, estrecha la pluma (b ≈ 9.7 m) y acerca R_Q. La clase (d) solo entra con evidencia independiente.

## Extensión de sens_xu2012 (pre-registrada antes de correr, 2026-11-05, aprobada por Esteban)

La corrida de 24 h salió inválida por no ser estacionaria (2.6 % entre 12–18 h y 18–24 h, contra un umbral del 2 %). Se repite **idéntica** salvo `endTime` = 36 h (129600 s), como el caso `sens_xu2012_36h` (`--run-suffix _36h`).

- **Reglas sin cambios:** misma validez del monitor, mismo umbral de estacionariedad del 2 %, mismos umbrales de detección (ΔR_Q 0.03, Δb 2.5 m) y misma comparación contra la producción.
- **Lo único que se desplaza es la ventana, 12 h:** X = 24–36 h, mitades 24–30 y 30–36 h, σ con ventanas horarias de 24 a 36 h (`metric_windows(t0=24)`).
- **Comparación:** contra `a_line2d_twolayer_q0.03` en su ventana original de 12–24 h (estacionaria). Esto supone que la producción ya está en estado estacionario a partir de las 12 h, lo que su estacionariedad respalda.
- **Si sigue sin ser estacionaria a las 36 h:** se reporta como inválida y no se extiende otra vez. Una deriva lenta del cierre de Xu en 2D sería en sí misma un resultado.
- **Si es válida:** sigue siendo una sensibilidad *conjunta* de cierre y arrastre de pared, como ya se explicó arriba. No se atribuye solo al cierre.
