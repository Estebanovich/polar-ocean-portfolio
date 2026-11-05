# PD3 — ¿`addMass < 0` retira agua local? (Task 7, Step 1)

**Veredicto: NO, si `temp_addMass`/`salt_addMass` están definidos.** Con un valor
escalar definido, la celda con `addMass < 0` retira agua **a `temp_addMass` y
`salt_addMass`**, no a la T y S locales: su T y S cambian (y se alejan de la
salida de la pluma). Solo si ambos quedan **sin definir** (`UNSET_RL`) el término
desaparece y la masa entra **y** sale con la T y S locales, en todas las celdas.

## 1. Fuente (MITgcm en `/Volumes/ESTEBAN_MIT/MITgcm`, sin modificar)

- `model/src/apply_forcing.F:504-530` (T) y `:873-898` (S): si
  `selectAddFluid ≠ 0` y `temp_addMass ≠ UNSET_RL`, con superficie libre lineal y
  `convertFW2Salt = -1` (lo fija `ini_parms.F:651` cuando `selectAddFluid ≥ 1`):

      gT += addMass·mass2rUnit·(temp_addMass − θ) / (rA·drF·hFacC)

  La fórmula no distingue el signo de `addMass`.
- La advección es en forma advectiva: `pkg/generic_advdiff/gad_calc_rhs.F:764-779`
  resta `θ_local·div(U)`, así que la divergencia que genera `addMass`
  (`calc_div_ghat.F:133`, `integr_continuity.F:130`) no cambia θ por sí sola.
  Entonces, en la celda, V dθ/dt = (m/ρ0)(temp_addMass − θ). Con m < 0, esto es
  exactamente retirar agua de temperatura `temp_addMass`. Si `temp_addMass` < θ, la
  celda se calienta, y viceversa: es antidifusivo y saca a θ del rango inicial
  (viola C5').
- Sin definir: `set_defaults.F:264-267` deja `temp_addMass`, `salt_addMass` y
  `temp_EvPrRn` en `UNSET_RL`, pero **`salt_EvPrRn = 0`** (`:265`), y
  `ini_parms.F:1579-1580` copia `*_EvPrRn` en `*_addMass` si están sin definir.
  Para que S también sea local hay que poner `salt_EvPrRn = 1.234567E+05`
  (= `UNSET_RL`, `eesupp/inc/EEPARAMS.h:81`). `config_summary.F:416-421`
  documenta "UNSET = use local T/S".
- (`external_forcing.F:383-401` tiene la misma fórmula, pero solo se usa con
  `USE_OLD_EXTERNAL_FORCING`, `apply_forcing.F:448-461`.)

## 2. Spike (caja 10 × 10 × 10, hidrostática, np = 1)

- Config: `exp1_fjord_plumes/config/spike_addmass/` (`code/`, `input/`,
  `gendata.py`, `analyse.py`). Build:
  `/Volumes/ESTEBAN_MIT/portfolio_builds/exp1/spike_addmass/`. Corridas:
  `/Volumes/LaCie/portfolio_runs/exp1/spike_addmass/{tset,local}` (OK, 300
  pasos de 60 s, 2 s de reloj cada una).
- dx = dy = 200 m, dz = 10 m, EOS lineal como en (a). T = 2 °C y S = 34 psu
  uniformes, difusión cero. Fuente de +2 m³/s en (k, j, i) = (2, 4, 7) y sumidero
  de −2 m³/s en (7, 4, 2) (`addMass` = ±ρ0·Q). Con campos uniformes, cualquier
  cambio de T o S en el sumidero solo puede venir del término `addMass`.
- Variante `tset`: `temp_addMass = 10`, `salt_addMass = 30`.
- Variante `local`: `temp_addMass` sin definir, `salt_EvPrRn = UNSET_RL`. La T y
  S de la salida se imponen con RBCS **solo en la celda fuente**, con
  máscara = Q·τ/V = 1 y τ = 2e5 s. Así, la tendencia RBCS −(máscara/τ)(θ − 10)
  (`pkg/rbcs/rbcs_add_tendency.F:88-111`) es idéntica a m(10 − θ)/(ρ0 V).

| t (h) | 0 | 0.5 | 1 | 2 | 3 | 5 |
|---|---|---|---|---|---|---|
| T sumidero `tset` (°C) | 2.000 | 1.9405 | 1.9526 | 1.9649 | 1.9540 | 1.9548 |
| S sumidero `tset` (psu) | 34.000 | 34.0297 | 34.0237 | 34.0175 | 34.0230 | 34.0226 |
| T sumidero `local` | 2.000 | 2.000 | 2.000 | 2.000 | 2.000 | 2.000 |
| S sumidero `local` | 34.000 | 34.000 | 34.000 | 34.000 | 34.000 | 34.000 |

- `tset`: el sumidero se enfría y sala (Δθ = −0.06 °C en 0.5 h). Sin la
  renovación por advección, la predicción de "retiro a `temp_addMass`" daría
  1.928 °C. La celda sale del rango inicial: `dynstat_theta_min` = 1.9548.
- `local`: |ΔT| = |ΔS| = **0 exacto** en el sumidero (`dynstat_theta_min`
  = 2.0000000000000).
- La celda fuente es igual en las dos variantes (máx. |ΔT| = 7.2e-5 °C,
  |ΔS| = 3.6e-5 psu en 5 h). Esto confirma que RBCS con máscara Q·τ/V
  reproduce la fuente a `temp_addMass`.
- Volumen: el flujo neto por las caras del sumidero es −2.000000 m³/s (= −Q),
  y Σ η·A se mantiene en |·| < 1e-10 m³ (la suma neta de addMass es 0). El
  volumen se conserva.

## 3. Alternativas para (c) (decisión de Esteban)

- **(i) Solo fuente en la salida + relajación RBCS en la cara** (alternativa
  del plan): `addMass = +Q0` en la salida, a `temp_addMass`/`salt_addMass`, sin
  sumideros, y RBCS hacia el perfil inicial en la columna de la cara. Es
  sencillo, pero el arrastre (≈ 3600 m³/s para Q0 = 100) no sale del
  fiordo como volumen. RBCS solo restaura T y S, no genera el flujo de entrada
  profundo. La circulación forzada por la pluma queda ≈ 37 veces más débil.
- **(ii) Copia local de `apply_forcing.F` en `code/`**: cambiar `addMass` por
  `MAX(addMass, 0)` en las dos fórmulas (T y S). Es un cambio de dos líneas,
  pero hay que mantener una copia de unas 1000 líneas de la fuente, y la T/S de
  salida sigue siendo un escalar.
- **(iii) Recomendada: `temp_addMass`/`salt_addMass` sin definir + RBCS con
  máscara solo en la(s) celda(s) de salida.** Los sumideros retiran agua local
  (probado: ΔT = ΔS = 0). La salida recibe la tendencia exacta
  m(T_out − θ)/(ρ0 V) con máscara = q·τ/V, `relaxTFile`/`relaxSFile` =
  T_out/S_out de `plume_source_sink`. No toca la fuente ni copia código, admite
  T/S por celda y es la misma ecuación que (ii), verificada en el spike (7e-5 °C).
  Requisitos: `useRBCS`, `pkg/rbcs` en `packages.conf`, `salt_EvPrRn=1.234567E+05`
  y τ ≤ min(V/q) para que la máscara quede ≤ 1 (no es obligatorio, pero es más
  claro). Ojo: si (c) usa RBCS también para otra cosa (p. ej. una esponja),
  comparte `tauRelaxT`, y la máscara se escala en consecuencia.

**Esperando decisión de Esteban antes de Step 4.**
