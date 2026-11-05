# PD2 (fallo R20): convergencia y sensibilidad de (a), plan previo a correr

**Estado: preparado (casos, configuraciones, builds); NADA corrido.** La campaña de (a) sigue en curso.

## Motivo

En la evaluación provisional (twolayer q0 = 0.01 y 0.03), C2 y C2'' pasan, pero **C2' falla**:

- la velocidad media de la pluma en su ancho, Q/b, es 0.52–0.56 de la teoría;
- la velocidad del núcleo w es 0.88–0.95 de la teoría;
- el ancho de la pluma es de 20–24 m contra 7.5 m de la teoría, unas 1.5 celdas de dx0 = 5 m.

Atribución propuesta: resolución en la cara más la pared sin deslizamiento (clase c), posiblemente con entrainment efectivo extra (clase b). Se **prueba** esta atribución; la métrica no se redefine.

## Casos (todos `twolayer`, q0 = 0.03 m² s⁻¹, 24 h de modelo, np = 1)

| caso | config / build | dx0 (m) | dz (m) | nx × nz | deltaT (s) | cambia respecto a producción |
|---|---|---|---|---|---|---|
| `conv_dx10` | `conv_dx10` | 10 | 10 | 170 × 50 | 5.0 | malla más gruesa |
| `conv_dx5` | `conv_dx5` | 5 | 5 | 193 × 100 | 2.5 | solo dz (producción tiene dz = 10) |
| `conv_dx2p5` | `conv_dx2p5` | 2.5 | 2.5 | 216 × 200 | 1.25 | malla más fina |
| `sens_freeslip` | `sens_line2d` (malla de producción) | 5 | 10 | 193 × 50 | 3.0 | `no_slip_sides = no_slip_bottom = .FALSE.` |
| `sens_xu2012` | `sens_line2d` | 5 | 10 | 193 × 50 | 3.0 | mezcla de Xu et al. 2012 (abajo) |
| (referencia) `a_line2d_twolayer_q0.03` | `a_line2d` | 5 | 10 | 193 × 50 | 3.0 | producción, 48 h (se compara en la misma ventana de tiempo) |

Lo que se mantiene igual en todos: estiramiento r = 1.03, dx_max = 200 m, L = 20 km, H = 500 m, la columna de tierra al final, perfiles, fuente, EOS, icefront, OBCS, cg3d (500 iteraciones, 1e-8) y salidas. Los builds van en `/Volumes/ESTEBAN_MIT/portfolio_builds/exp1/<config>/`.

- **dt por CFL:** deltaT = 0.5 s m⁻¹ × min(dx0, dz). La producción dio CFL_w = 0.198 con w_max = 0.66 m/s y dt/dz = 0.3 s/m. Con dt/min(Δ) = 0.5 s/m, CFL ≈ 0.33 si w no cambia, y ≤ 0.45 si el núcleo llega a ~0.9 m/s. Guardia: si el monitor da CFL > 0.5 en la primera hora, se para y se reporta.
- **Escalado comprobado por prueba** (`tests/test_pd2_cases.py`):
  - `icefrontlength` = 1/dx0 (0.1, 0.2 y 0.4 m⁻¹) solo en la columna 0;
  - `R_icefront` = H solo en la columna 0, así que `K_icefront` = Nr y el factor parcial del fondo vale 1 (`icefront_thermodynamics.F:158-174`);
  - fuente ρ0·q0·dy = 308.1 kg s⁻¹ en la celda del fondo (k = nz − 1, i = 0);
  - `SIZE.h` coincide con la malla.
- **`sens_xu2012`:** Xu et al. (2012, p. 3, Ann. Glaciol. 53(60)) usan νv = 0.1 m² s⁻¹, biarmónica horizontal 300 m⁴ s⁻¹ (con dx = 20 m) y κ = ν (Pr = 1). Se pone `viscAr = diffKrT = diffKrS = 0.1` y `viscAh = diffKh = 0`. La biarmónica literal es inestable aquí: el límite es L⁴/(32 dt) = 5⁴/96 = 6.5 m⁴ s⁻¹ (`pkg/mom_common/mom_calc_visc.F:48-49`). Por eso se escala a dx0 = 5 m con el mismo tiempo de amortiguamiento de malla dx⁴/A4: A4 = 300·(5/20)⁴ = **1.171875 m⁴ s⁻¹** en `viscA4`, `diffK4T` y `diffK4S`. Para la biarmónica, `SIZE.h` necesita OLx, OLy ≥ 3 (`config_check.F:566-575`); se cumple (3). Pared sin deslizamiento como en producción.

## Costo estimado

Base: 65 s por 2400 pasos en 9 650 celdas → c = 2.81e-6 s por celda·paso. Costo = c × celdas × pasos. Las iteraciones de cg3d (350–390 en producción) deberían crecer al refinar, más o menos con el número de puntos en la dirección refinada. Por eso se da un rango.

| caso | celdas | pasos (24 h) | base | con crecimiento de cg3d |
|---|---|---|---|---|
| `conv_dx10` | 8 500 | 17 280 | 7 min | 3–7 min (×0.5–1) |
| `conv_dx5` | 19 300 | 34 560 | 31 min | 31–62 min (×1–2) |
| `conv_dx2p5` | 43 200 | 69 120 | 2.3 h | 2.3–9.3 h (×1–4) |
| `sens_freeslip` | 9 650 | 28 800 | 13 min | ≈ 13 min |
| `sens_xu2012` | 9 650 | 28 800 | 13 min | ≈ 13 min |
| **total** | | | **3.4 h** | **3.4–10.9 h** |

`conv_dx2p5` podría pasar el límite de 8 h por corrida. Recomendación: cronometrar antes 1 h de modelo de `conv_dx2p5` (6–23 min). Si la extrapolación pasa de 8 h, se corre en dos tramos de 12 h con pickup; esto es legítimo porque C3 (reinicio bit a bit) está verificado.

## Criterios, fijados antes de correr

Métricas: el mismo código y la misma ventana de promedio que la evaluación provisional de PD2. Para cada caso, en la misma ventana de tiempo:

- R_Q = (Q/b)_modelo / (Q/b)_teoría;
- R_w = w_núcleo / w_teoría;
- b en metros y en celdas (b/dx0).

Ruido: la desviación estándar σ de las medias horarias dentro de la ventana. Una diferencia cuenta si es > 2σ.

**La atribución a la resolución (clase c) se apoya si se cumplen las tres:**

1. R_Q crece de forma monótona en dx10 → dx5 → dx2p5, con cada paso > 2σ;
2. b (en m) disminuye de forma monótona al refinar, es decir, b sigue a la malla (b/dx0 casi constante): el ancho es numérico;
3. la extrapolación de Richardson de R_Q (con los tres niveles) queda en [0.75, 1.25], la banda de C2'.

**La contribución de la pared sin deslizamiento (clase c) se apoya** si `sens_freeslip` sube R_Q en más de 2σ respecto a producción en la misma ventana.

**La atribución se refuta si:**

- R_Q no cambia (dentro de 2σ) entre dx10, dx5 y dx2p5;
- o R_Q no es monótona;
- o b no disminuye al refinar.

En ese caso el déficit no es de resolución. Se evalúa entonces entrainment extra o física 2D (clase b), o una discrepancia de la definición de Q/b frente al sombrero de copa de la teoría. Esta última se clasificaría (d), pero solo con evidencia independiente; la métrica no se cambia.

**Apoyo parcial:** R_Q crece de forma monótona, pero el límite extrapolado es < 0.75. Se reporta la parte numérica cuantificada y el resto queda sin explicar (candidato clase b).

**`sens_xu2012`:** si R_Q cambia más de 2σ, el resultado depende del cierre de mezcla (elección numérica, clase c). Se reporta como incertidumbre de modelo. No se adopta la mezcla que "mejor acierta" (regla: no ajustar parámetros para que coincida un caso).

Además, en todos los casos, C2 (despegue ±0.10 H) y C2'' (derretimiento, factor 1.5) deben seguir dentro de tolerancia. Si un refinamiento los rompe, eso también se reporta.

## Cómo lanzar (después de la campaña de (a), uno a la vez)

```
cd ~/Phd_project/portfolio && export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
P=/opt/anaconda3/envs/pyocean/bin/python; S=exp1_fjord_plumes/scripts/run_campaign.py
$P $S --config conv_dx10 --cases conv_dx10 --np 1
$P $S --config conv_dx5 --cases conv_dx5 --np 1
$P $S --config sens_line2d --cases 'sens_*' --np 1
$P $S --config conv_dx2p5 --cases conv_dx2p5 --np 1 --end-time 3600 --run-suffix _probe1h   # cronometrar primero
$P $S --config conv_dx2p5 --cases conv_dx2p5 --np 1
```
