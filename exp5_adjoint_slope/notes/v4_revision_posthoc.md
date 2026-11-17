# V4 y GM: revisión post hoc (registrada antes de correr)

Fecha: 2026-11-16. Decidido por Esteban tras la puerta del día 1 y la falla de V4. **Todo lo de esta nota es post hoc**: se decidió después de ver fallar V4 en G0 y el adjunto en GS. El README lo declara así.

## 1. GM en forma skew-flux (decisión B)
- **Motivo:** en forma advectiva (`GM_AdvForm`), el adjunto de Tapenade crece ×10¹³ por paso (×300 con `GM_Small_Number = 1e-8`) y llega a NaN en GS. En forma skew-flux es estable en 10 días con el valor por defecto. Es la forma que usan todos los `verification/*/code_tap` de MITgcm (ver `puerta_dia1.md`).
- **Cambio en J:** J = −ρ₀c_p [Σ v·hFacS·Δz·(θ̄ − θ_f)·Δx + Σ F^GM_y]·10⁻¹² TW, donde F^GM_y es el flujo meridional de θ del tensor GM/Redi integrado en la cara (`gmredi_ytransport.F`):
  - F^GM_y = −hFacS·Δz·Δx·[K_vy·(θ_j − θ_{j−1})/Δy_C + K_vz·∂θ/∂z|_v].
  - En forma skew, ese flujo contiene el bolo y la difusión isopicnal de Redi. La J anterior solo tenía el bolo, y por eso el balance de calor de la plataforma no cerraba.
  - Con la nueva J, la comparación con el fino (v′θ′ resuelto) es más justa: en los dos casos J es el transporte total de calor por remolinos más la media.
- `GM_Small_Number` vuelve al valor por defecto (1e-20).

## 2. Revisión del forzamiento superficial (V4)
**Diagnóstico de G0:**
- el CDW inundó la plataforma (0.76 °C en el fondo para y = 5–85 km);
- la SST de la plataforma está en −1.7 °C contra un objetivo de −1.85;
- el objetivo de salinidad superficial sobre la plataforma (34.3) es **más salado** que el de mar abierto (34.2), al revés de una plataforma fresca (Thompson et al., 2018).

**Mecanismos de la literatura ya revisada** (`rev_literatura.md`):
- más viento del este sobre el talud hunde las isopicnas costeras y enfría la plataforma (Stewart & Thompson, 2012; Morrison et al., 2023);
- la plataforma fresca necesita agua superficial dulce (Thompson et al., 2018).

**Variantes de cribado** (10 años desde las condiciones iniciales, GM skew, todo lo demás igual a G0):

| Id | τ_E (y < 150 km) | Objetivo de S superficial en la plataforma |
|---|---|---|
| S0 | −0.05 N m⁻² (sin cambio) | sin cambio (34.3) |
| SW | **−0.10 N m⁻²** | sin cambio |
| SF | −0.05 N m⁻² | **33.8**, con la misma mezcla tanh hacia 34.2 en mar abierto |
| SWF | −0.10 N m⁻² | 33.8 |

Las condiciones iniciales no cambian. Solo cambia el forzamiento: el viento y el objetivo de relajación de S.

**Criterio de selección (fijado antes de correr)**, sobre la media del año 10:
1. θ media en volumen de la plataforma (y < 80 km) < −1.0 °C;
2. J (nueva definición) > 0.

Se elige la **primera** variante que cumpla ambas condiciones, en este orden: S0, SW, SF, SWF. Es decir, primero el menor cambio, y el viento antes que el agua dulce porque el viento es el control que estudia el experimento. Si ninguna cumple, G0 se gira con SWF y **V4 se declara fallida**. Se reportan las cuatro variantes, no solo la elegida.

Después se gira G0 completo (30 años) con la variante elegida, y V4 se evalúa sobre G0 con el criterio congelado del PD0.

## 3. Hallazgo (2026-11-16, después de registrar y antes de elegir): el dominio era periódico en y
- **El primer S0 (10 años) no cumplió:** plataforma a −0.80 °C (0.92 °C en el fondo) y J = −2.0 TW. Pero el balance de calor al sur de la línea no cerraba: J llegaba a −9 TW mientras dH/dt era de ±1 TW por año.
- **Causa:** cruzaban 0.168 Sv netos hacia el norte por **todas** las filas. La batimetría no tenía tierra en la primera ni en la última fila, y MITgcm es periódico en y salvo que la tierra lo cierre. El océano de 3000 m del norte se conectaba con la plataforma de 500 m del sur: el CDW entraba a la plataforma **a través de la frontera**, no a través del talud.
- **Clasificación:** (a), error nuestro (gendata). Invalida G0, el diagnóstico de V4 de G0 (plataforma "cálida") y el primer S0. Afecta también al fino `chan2p5`, que usa la misma batimetría.
- **Arreglo:** tierra en las filas 0 y ny−1 (`fields.bathymetry`, con prueba); el canal sigue periódico en x.
- **Consecuencia para la revisión:** la falla de V4 que motivó revisar el forzamiento se debía, al menos en parte, a este error. El orden registrado ya lo contempla: primero se repite **S0 sin cambios de forzamiento**. Si cumple, no se cambia el forzamiento. Los resultados del primer S0 se reportan como invalidados por el error, no como evidencia física.

## 4. La forma skew-flux es inestable en el directo: se pasa a la opción A (Esteban, 2026-11-16)
- **Hallazgo:** con las paredes puestas, S0 en forma skew se volvió numéricamente inestable desde el año 2. θ llegó a ir de −15 a +12 °C, con J oscilando ±5 TW entre años. El mismo arranque, en corridas de 2 años:

  | Variante | θ máx / mín |
  |---|---|
  | V1: forma advectiva, DST3 sin limitador | 0.807 / −1.79 (acotada) |
  | V2: forma skew, DST3 con limitador (esquema 33) | 2.38 / −7.7 (inestable) |
  | V3: forma advectiva con `GM_Small_Number = 1e-6` | 0.801 / −1.79 (acotada) |

  La inestabilidad viene de la forma skew, no de la advección.
- **Decisión (opción A):** forma advectiva con `GM_Small_Number = 1e-6`. Es la única combinación con directo estable (V3) y adjunto estable: 10 días sin crecimiento, en la prueba del 15 nov sin paredes, que se repite con paredes en GS. El piso equivale a N² ~ 1e-8 s⁻² y solo toca celdas casi neutras. Se declara como regularización, igual en directo y adjunto.
- **J (definitiva):** J = −ρ₀c_p [Σ (v + v*)·hFacS·Δz·(θ̄ − θ_f)·Δx + Σ F^Redi_y]·10⁻¹² TW, es decir, resuelto + bolo (`GM_PsiY`) + flujo de Redi (K_vy y K_vz; `GM_ExtraDiag` activo en forma advectiva).
  - En las medias mensuales, las tres partes salen exactas de `VVELTH`/`VVEL`, `GM_vbT` y `DFyE_TH`.
  - El ancla del adjunto compara la J de Fortran con la de Python, rearmada desde `GM_PsiY`, `GM_Kvy` y `GM_Kvz`, y con los flujos del propio modelo.
- La sección 1 (forma skew) queda **reemplazada**. El primer S0 con paredes (en forma skew) también se invalida.

## 5. Resultado de V4 sobre G0 con la configuración Dettling (2026-11-16)
G0 giró 30 años (np 5), estable: θ entre −1.66 y 0.83 °C en el año 30. Media de los 12 meses del año 30 (`G0/screen_v4.json`):

| Magnitud | Valor |
|---|---|
| θ medio de la plataforma (y < 80 km) | −1.48 °C (fondo −1.47 °C); S 34.44 |
| J total | **−0.030 TW** (hacia mar afuera) |
| resuelto / bolo / Redi | −0.035 / +0.002 / +0.003 TW |

**V4 no se cumple.**
- "Frente de plataforma fresca": falla por construcción (régimen de plataforma densa, `pd1b_config_dettling.md`).
- J > 0: también falla. J es negativa desde el año ~17.

**Balance de calor al sur de la línea J (y < 170 km), por años:**

| Año | J (TW) | dH/dt (TW) | θ plataforma |
|---|---|---|---|
| 8 | +0.134 | +0.150 | −1.52 |
| 14 | +0.059 | +0.054 | −1.41 |
| 20 | −0.034 | −0.034 | −1.41 |
| 30 | −0.030 | −0.028 | −1.47 |

- **El signo de J está bien.** dH/dt ≈ J: el contenido de calor al sur de la línea cambia en lo que J dice, con un flujo neto de volumen < 10⁻⁶ Sv.
- **Lo que falta.** El flujo superficial al sur de la línea es casi nulo (|dH/dt − J| ≲ 0.006 TW). La relajación a θ_f = −1.9 °C no actúa como sumidero, porque la capa superficial ya está cerca de −1.8 °C. Sin sumidero en la plataforma, el estado de equilibrio solo admite J → 0. La plataforma pierde hoy el calor que ganó en los años 5–15 (J decae de −0.036 a −0.030 TW entre los años 23 y 30).
- **Clasificación (provisional):** (b) física ausente. Falta el hielo marino de S&T y Dettling, que extrae calor en la plataforma: el flujo de dos ecuaciones contra el punto de congelación in situ, más la convección de polinia que lleva agua cálida a la superficie. Además, el estado no ha llegado al equilibrio.
- **Del orden de lo publicado:** no se evalúa, porque el signo ya falla.
- **No se cambia el criterio ni el forzamiento.** El adjunto (∂J/∂u) sigue definido con J < 0, así que GS, GA y GB siguen adelante. El README reporta V4 como fallo, con este balance.

## 6. Corrección de §5: el fallo de V4 era nuestro error (a), no física ausente (2026-11-16)
- **Hallazgo.** El flujo superficial casi nulo de §5 no venía de la física. G0 corrió 30 años **sin relajación superficial**:
  - MITgcm activa `doThetaClimRelax` y `doSaltClimRelax` solo si `tauThetaClimRelax` y `tauSaltClimRelax` son > 0 (`model/src/set_parms.F:250-253`). `FORCING_SURF_RELAX` solo se llama con esas banderas (`external_forcing_surf.F:173`).
  - Al pasar a Dettling se quitaron los τ y quedaron solo `lambdaThetaFile` y `lambdaSaltFile`, que se leen pero no se aplican.
  - La prueba `test_dettling_physics` exigía que los τ no estuvieran, así que fijaba el error.
- **Corrección.**
  - Se añaden `tauThetaClimRelax = 164000` y `tauSaltClimRelax = 2592000`; las tasas efectivas siguen saliendo de los archivos λ. Esto vale en chan10 y en chan2p5.
  - Nueva prueba: `test_surface_relaxation_switched_on`.
  - Una corrida de 30 días confirma `doThetaClimRelax = T` y `doSaltClimRelax = T`, y que el modelo abre los 14 archivos de entrada.
- **Auditoría del eco de parámetros:** viscosidades, Leith, arrastre, κ_z, EOS, f₀, GM (κ, mapa, taper, piso), rbcs y forzamientos coinciden con `pd1b_config_dettling.md`.
- **Consecuencias.**
  - G0 sin relajación queda como `G0_norelax`; su V4 (§5) se invalida.
  - GS se detuvo, porque usaba el mismo estado.
  - G0 se gira de nuevo 30 años y V4 se evalúa sobre él. El criterio no cambia.
- **Lo que no cambia:** el ancla GJ. Compara Fortran con Python sobre las mismas corridas, así que no depende de la relajación. Se repite igual con el G0 nuevo.

## 7. V4 sobre G0 con relajación (2026-11-16)
G0 giró de nuevo 30 años, estable: θ entre −1.94 y 0.84 °C. La plataforma quedó en el punto de congelación: θ medio −1.897 °C, fondo −1.898 °C, S 34.38.

**Año 30** (`G0/screen_v4.json`):

| Magnitud | Valor |
|---|---|
| J total | **−0.0002 TW** |
| resuelto (media euleriana) | −0.0022 TW, hacia mar afuera |
| bolo | +0.0017 TW, hacia la costa |
| Redi | +0.0002 TW, hacia la costa |
| J por años (14, 20, 26, 30) | −0.007, −0.003, −0.001, −0.0002 TW, convergiendo a 0 |
| dH/dt al sur de la línea | −0.0035 TW (casi en equilibrio) |

- **V4 no se cumple como está escrito.** J no es > 0: en el año 30 es ≈ 0 y su signo no es robusto.
- **Lo que publican Dettling et al. (2024, p. 1679), con esta misma configuración:** "at both resolutions, the eddy heat flux is balanced by a mean offshore heat flux over the shelf and slope". El G0 reproduce esa estructura:
  - el flujo de los remolinos parametrizados (bolo + Redi = +0.0019 TW) va hacia la costa;
  - la media euleriana (−0.0022 TW) lo compensa hacia mar afuera;
  - la J neta tiende a 0 a medida que la plataforma se acerca al punto de congelación.
- **Clasificación: (a) nuestro error, en la expectativa, no en el modelo.** "J > 0" se fijó para el régimen de plataforma fresca (Thompson et al., 2018). Al pasar a la plataforma densa (`pd1b`) no se revisó si seguía aplicando, y en este régimen la J neta en equilibrio es ≈ 0. El criterio no se redefine: V4 se reporta como no cumplido, con esta explicación y la partición.
- **Consecuencia para el adjunto.** J sigue siendo un funcional bien definido: es el calor neto que entra a la plataforma, igual en equilibrio a la pérdida superficial ahí. Pero es un residuo pequeño de dos términos opuestos de ±0.002 TW, así que las sensibilidades conviene leerlas también por componente. La interpretación de H5 se discute con esa partición.
- **Ancla GJ** (desde el G0 nuevo):
  - fc = 4.93293e-5 TW contra Python 4.93153e-5 TW. Diferencia absoluta 1.4e-8 TW, igual que antes de la corrección.
  - rel_diff = 2.8e-4 y 4.1e-5 con el último registro extrapolado. La relación respecto de J está mal condicionada porque J ≈ 0.
  - Respecto del flujo bruto (|resuelto| + |bolo| + |Redi| = 0.0035 TW), la diferencia es 3.9e-6, y 5.8e-7 con el último registro extrapolado.
  - Python y los flujos del propio modelo coinciden a 1e-15.
