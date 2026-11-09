# PD1a — Cronometraje, paso de tiempo y calibración de ΓT (Exp. 2, Task 4 fase A)

Fecha: 2026-11-08. Criterios vinculantes: `pd0_preregistro.md` (K1, V1, N1–N3). Presupuesto B1: ≤ 32 h de cómputo para todo el Exp. 2 y ≤ 8 h por etapa de corrida. Corridas en `/Volumes/LaCie/portfolio_runs/exp2/`, logs en `…/exp2/logs/`.

## 1. Corridas de cronometraje

Todas con hilos = 1, `pgrep -fl mpirun` vacío antes y después, y ΓT = 0.03 (el valor no afecta el costo). El costo por paso es el "Wall clock time" de la sección `MAIN LOOP` del propio modelo, dividido por el número de pasos.

| Corrida | Malla | np × simultáneas | dt (s) | Duración | s/paso |
|---|---|---|---|---|---|
| `cal2_gT0.03_smoke30d_walls` (Task 2) | 2 km | 5 × 1 | 200 | 30 d | 0.01955 |
| `cal2_gT0.03_tim30d_dt400` | 2 km | 5 × 1 | 400 | 30 d | 0.02040 |
| `tim4m_gT0.03_tim10d` | 4 km | 5 × 1 | 400 | 10 d | 0.00580 |
| `tim4m_gT0.03_tim30d_dt400` | 4 km | 5 × 1 | 400 | 30 d | 0.00660 |
| `tim4m_gT0.03_tim30d_dt800` | 4 km | 5 × 1 | 800 | 30 d | 0.00619 |
| `ens_T*_zc400_s1.0_tim10d` (4 de 5, ver §5) | 4 km | 1 × 4 | 400 | 10 d | 0.02125 c/u |
| `ens_T*_s0.8_tim10d` | 4 km | 1 × 5 | 400 | 10 d | 0.02385 c/u |

- La caja np = 5 a 4 km usa `config/ocean0_4km_mpi` (caso nuevo `tim4m_gT0.03`, misma física que `ocean0_4km`).
- El costo por paso de los primeros 10 días subestima el de 30 días (0.00580 contra 0.00660 a np = 5: factor 1.14, más iteraciones de cg2d). Para extrapolar se aplica ese factor a las corridas serie de 10 días.
- **Rendimiento a 4 km** (dt = 800, por corrida de 2 años, con el factor 1.14): np = 5 en serie, 0.136 h por corrida; np = 1 × 5 simultáneas, 0.594 h por tanda de 5, es decir **0.119 h por corrida**. El arreglo de mayor rendimiento es np = 1 con 5 simultáneas (13 % mejor).

## 2. Paso de tiempo por CFL

**Regla** (escrita después de ver los CFL a dt = 200/400 s y antes de correr los dt dobles): dt = el mayor de la sucesión de duplicaciones desde el valor inicial (200 s a 2 km, 400 s a 4 km) cuyo CFL máximo en 30 días sea ≤ 0.5 en advección (u, v, w con hFac) y en ondas internas. El barotrópico es implícito. Monitor cada hora en las corridas de dt doble y diario en las de Task 2.

| Malla, dt | CFL u | CFL v | CFL w (= W_hf) | CFL ondas internas | f·dt |
|---|---|---|---|---|---|
| 2 km, 200 s (30 d) | 0.084 | 0.057 | 0.134 (día 22) | 0.065 | 0.028 |
| 2 km, **400 s** (30 d) | 0.170 | 0.114 | **0.271** (día 22) | 0.131 | 0.056 |
| 4 km, 400 s (30 d) | 0.094 | 0.053 | 0.159 (día 26) | 0.065 | 0.056 |
| 4 km, **800 s** (30 d) | 0.186 | 0.106 | **0.317** (día 30) | 0.131 | 0.112 |

- **CFL de ondas internas:** c₁·dt/Δx, con c₁ = velocidad del primer modo baroclínico, calculada del perfil discreto de 36 niveles con la EOS lineal (problema de autovalores N²w = c²·(−d²w/dz²)). Ocean0: c₁ = 0.654 m s⁻¹ (H = 720 m, igual a NH/π). Ensamble (rampa, d = 300 m): el máximo es c₁ = 0.91 m s⁻¹ (Tbot = −0.5, z_c = 400), que da un CFL de 0.091 (4 km, 400 s) y 0.18 (4 km, 800 s).
- **Duplicar otra vez** (800 s a 2 km, 1600 s a 4 km) daría un CFL w ≈ 0.54 y 0.63, por encima de 0.5. Ese dt queda descartado sin correrlo, porque el CFL escala linealmente con dt.
- **Elección: dt = 400 s a 2 km y 800 s a 4 km.** dt/Δx es el mismo en las dos mallas. Sustituye los valores provisionales de la Task 2 (200 y 400 s) en `cases.yaml`.
- **Sensibilidad al dt** (30 días, ΓT = 0.03, `pd1_smoke_check.py`):

  | | dt | melt plataforma (m a⁻¹) | melt profundo (m a⁻¹) | S1 / S2 / S4 / S5 | S3 |
  |---|---|---|---|---|---|
  | 2 km | 200 | 19.827 | 57.781 | PASS (S2 2.31e-6) | FAIL |
  | 2 km | 400 | 19.836 (+0.05 %) | 57.811 (+0.05 %) | PASS (S2 2.30e-6) | FAIL (T mín −2.507) |
  | 4 km | 400 | 20.498 | 61.327 | PASS (S2 1.77e-6) | FAIL |
  | 4 km | 800 | 20.488 (−0.05 %) | 61.265 (−0.10 %) | PASS (S2 1.75e-6) | FAIL (T mín −2.449) |

  S3 falla igual que en la Task 2 (R10: el criterio no se redefine). Las demás condiciones del smoke se mantienen con el dt nuevo.
- **Limitación:** el CFL w a 4 km y 800 s alcanza su máximo el día 30, el último de la corrida, así que puede crecer más adelante. Se vigila en las corridas de 2 años (§4).

## 3. Costo del Exp. 2 contra B1

Horas de pared con los s/paso medidos (2 km, 0.02040 s/paso a 400 s: 0.894 h por corrida de 2 años; 4 km, ver §1):

| Bloque | Corridas | Arreglo | Horas |
|---|---|---|---|
| Calibración 2 km | 4 × 2 años (3 + verificación) | np = 5, una a la vez | 3.58 |
| Calibración 4 km | 3 × 2 años simultáneas + verificación | np = 1 (×3, luego ×1) | ≤ 1.06 |
| Ocean1 2 km | 20 años = 4 etapas × 5 años (2.24 h c/u ≤ 8 h) + prueba de reinicio | np = 5 | 8.94 + ≤ 0.25 |
| Ensamble 4 km | 28 × 2 años, tandas 5,5,5,5,5,3 | np = 1 × 5 | 3.56 |
| Ya gastado | smokes de la Task 2 y cronometraje | — | ≈ 0.3 |
| **Total** | | | **≈ 17.7 h ≤ 32 h (B1 se cumple)** |

- Con una segunda iteración de K1 en las dos mallas (+0.89 + 0.53 h), el total sería ≈ 19.1 h.
- **Para comparar:** con los dt provisionales (200 y 400 s), el mismo plan costaría ≈ 33 h (6.85 + 2.1 + 17.1 + 7.0), por encima de B1.
- La cota de "calibración 4 km" usa el costo con 4 simultáneas (0.02125 s/paso × 1.14), que es una cota superior para 3 simultáneas y para una corrida serie sola.

**Arreglo de la calibración a 4 km:** np = 1, 3 simultáneas, y la verificación `ocean0_4km` también serie. np = 1 es el arreglo de mayor rendimiento medido (§1). Con solo 3 corridas, np = 5 una tras otra terminaría ≈ 0.1 h antes (3 × 0.136 h = 0.41 h contra ≈ 0.53 h), y la verificación sola ≈ 0.3 h antes. Se prefiere np = 1 porque así la calibración, la verificación y el ensamble usan el mismo binario. La verificación `ocean0_4km` es además la referencia del ensamble. La diferencia (≤ 0.5 h) no afecta a B1.

## 4. Calibración K1 (pendiente: cola en marcha)

- Cola: `exp2_isomip_cavity/scripts/queue_calibration.sh` (ver `task-4-report`).
- Procedimiento: `analysis/calibrate.py` aplica pd0 K1 de forma mecánica:
  - lee cada corrida con `metrics.load_run`, que aplica `require_ok` (N3);
  - calcula el melt profundo (< −300 m) de los registros 19–24;
  - interpola ln⟨m⟩ contra ln ΓT, lineal por tramos, hacia 30 m a⁻¹, sin extrapolar;
  - escribe `figures/calibration_<res>.json`.
- Resultados, K1, V1, N1–N3 y estadística de sobreenfriamiento (R15): se rellenan al terminar la cola.

## 5. Hallazgo aparte (no afecta a la calibración)

`ens_T0.0_zc400_s1.0_tim10d` no arrancó: `INI_THETA: found 10142 wet grid-pts with theta=0 identically` (STOP en `ini_theta.F`, `checkIniTemp`).
- Con Tbot = 0.0, el perfil rampa da T = 0 exactamente bajo z_c + d/2, y MITgcm lo toma como un campo sin inicializar.
- Afecta a los 3 casos del ensamble con Tbot = 0.0 (z_c = 250, 400, 550).
- No afecta a Ocean0 ni a Ocean1.
- Solución probable, que queda para decisión antes del ensamble (Task 5): `checkIniTemp=.FALSE.` en PARM05. Es solo una comprobación, sin efecto numérico.
- Esta corrida fallida no entra en el cronometraje: la tanda de 5 se repitió con los casos s = 0.8.
