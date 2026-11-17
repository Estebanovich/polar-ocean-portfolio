# PD0 — preregistro del Exp. 5 (criterios congelados antes de GA, GB, GV, FC y FP)

Fecha: 2026-11-15. Fuente: el spec de diseño del Exp. 5 (2026-11-14), §6–7. Los umbrales se copian **literalmente** del spec. Lo marcado *(propuesta)* lo fijé yo dentro de lo que el spec deja abierto; Esteban lo revisa, pero ningún número se evalúa antes de este commit.

## Cambios de implementación ya registrados (antes de evaluar nada)
- **Ventanas del adjunto:** GA = **360 días** (no 365) y GB = 180 días; J es la media de los últimos 30 días. Sin `pkg/cal`, los controles mensuales son periódicos sobre `externForcingCycle`, que tiene que ser múltiplo de 30 días. El tiempo reinicia en 0 desde el pickup final de G0 (`pickupSuff`).
- **Línea de J:** y = 90 km, en la isóbata de ~920 m (cara v común a las dos mallas).
- **G0:** viscAh = 1000 m² s⁻¹ y Δt = 1200 s, porque con 200 m² s⁻¹ divergió en el año 2.

## Validación

| Id | Prueba | Criterio (congelado) |
|---|---|---|
| V1 | `grdchk` adjunto contra diferencias finitas centradas y contra el tangente lineal, en ≥ 10 puntos: 4 de τx (registros 1, 6 y 12 de GA y uno sobre el quiebre), 3 de T₀ y 3 de κᵥ, donde la sensibilidad de GA es grande | error relativo < 1 % en **todos** los puntos |
| V2 | linealidad: corridas directas con +Δu y −Δu del mismo patrón, con la amplitud de V3 | \|Δ⁺ + Δ⁻\| / \|Δ⁺ − Δ⁻\| < 0.1 |
| V3 | perturbación finita guiada por el adjunto, tres patrones: el máximo de ∂J/∂τx, el de ∂J/∂T₀ y el de ∂J/∂κᵥ | \|ΔJ_pred − ΔJ_dir\| / \|ΔJ_dir\| < 0.30 |
| V4 | estado medio del grueso | cualitativo y declarado: frente de plataforma fresca presente y J > 0 en la media de G0 |
| V5 | fracción de remolinos de J en el fino (v′θ′ respecto de la media zonal y temporal) | > 50 %, o se declara que no |

- **eps de V1:** el que dé menor error en la prueba de humo GS, barriendo tres valores por control; se anota en `notes/puerta_dia1.md` antes de V1.
- **Amplitudes Δu (V2, V3)** *(propuesta)*: el patrón se normaliza a su máximo y se multiplica por:
  - 0.01 N m⁻² para τx (constante en el tiempo, como en H5b);
  - 0.1 °C para T₀;
  - 1×10⁻⁵ m² s⁻¹ para κᵥ.

## Hipótesis

- **H5a (rotacional contra magnitud).** Se comparan dos patrones de τx(y), constantes en el tiempo y con la misma norma RMS sobre y < 150 km:
  - (i) aumento uniforme del viento del este en y < 150 km;
  - (ii) dipolo centrado en el quiebre (y₀ = 100 km), derivada de gaussiana de semiancho 25 km.

  Se cumple si |ΔJ_pred(ii)| > |ΔJ_pred(i)|, con ΔJ_pred = Σ ∂J/∂τx · δτx sobre todos los registros de GA. Es **"robusta"** solo si (ii) gana también con semianchos de 12.5 y 50 km. La confirmación directa (V3 sobre (i) y (ii) de 25 km) se reporta, pero no cambia el veredicto.
- **H5b (prueba cruzada).** En el fino, el patrón de mayor sensibilidad de GB, constante en el tiempo y con amplitud 0.01 N m⁻². Se cumple si se dan las tres condiciones:
  1. **detectable:** |ΔJ| > 2σ_Δ, con σ_Δ = √(s_C²/N + s_P²/N);
  2. **signo** igual al que predice GB;
  3. **magnitud** dentro de un factor 2 de la predicha.

  Si no es detectable, se reporta "no detectable con N = 3", que no es confirmación. Con N = 2 se reporta solo como indicativa.
- **H5c (anticipación).** |∂J/∂τx| integrado en y, como función del tiempo de anticipación de cada registro mensual de GA. Se cumple si cae a la mitad de su máximo entre 30 y 180 días.

## Escalado de sensibilidades *(propuesta)*
σ típicas para comparar controles (Mazloff 2011; Jones et al. 2019):
- τx: 0.02 N m⁻² por registro mensual;
- T₀: 0.1 °C;
- S₀: 0.02;
- κᵥ: 1×10⁻⁵ m² s⁻¹.

## Enmienda PD1b (2026-11-16, antes de correr GA, GB, GV, FC y FP)
Cambian la geometría y el forzamiento (configuración de Dettling, `pd1b_config_dettling.md`). **Los umbrales de V1–V5 y H5a–H5c no cambian.** Cambia solo lo que depende de la geometría:
- **Línea de J:** y = 170 km, en la isóbata de ≈ 920 m (fila 17 del grueso y 68 del fino). Reemplaza la línea de 90 km.
- **H5a.** El quiebre está en y₀ = 125 km. Las escalas se toman de la batimetría nueva:
  - el semiancho del dipolo es la escala del tanh del talud, 37.5 km (antes ~25 km), y la robustez se prueba con la mitad y el doble: 18.75 y 75 km;
  - la norma RMS se toma sobre la plataforma y el talud, y < 275 km (antes y < 150 km);
  - el patrón (i) es el aumento uniforme del viento del este en y < 275 km.
  - La motivación es el bombeo de Ekman sobre el talud (S&T, p. 437), no el mecanismo de Rodriguez et al. (2016).
- **V4.** El resultado ya está reportado (`v4_revision_posthoc.md` §7): no se cumple, con la explicación. No se reevalúa.
- **eps de V1** (de `puerta_dia1.md`):
  - τx: ≤ 1e-5 N m⁻². Se toma el menor error entre 1e-5 y 1e-6, porque con 1e-4 la respuesta ya no es lineal en la fila de J.
  - T₀: 1e-3 °C.
  - κᵥ: 1e-5 m² s⁻¹.
- **J secundaria (post hoc, opcional).** El flujo hacia la costa de los remolinos parametrizados (bolo + Redi), la cantidad que valida Dettling et al. (2024). Si el tiempo alcanza, se corre como segundo funcional, con su propio binario, y se reporta como **post hoc**. J total sigue siendo el funcional principal.

## Enmienda PD2 (2026-11-16, 22:40; decisión de Esteban; **post hoc**, escrita antes de ver EA, GV y los ensambles directos)

**Hallazgo.** El grueso de 10 km es caótico (permite remolinos), al contrario de lo que suponía el spec (§4: "el grueso con GM no debería tener ese límite").
- Ensamble directo E0–E4 (360 días desde G0, ruido en θ₀ ≤ 10⁻⁶ °C): la desviación del J mensual es 3×10⁻⁸ TW en el mes 1, 6×10⁻⁵ en el mes 2 y 3×10⁻⁴ en el mes 3, es decir, del orden de J. El error crece con una escala de e-plegamiento de ~4 días.
- E0z (controles nulos por la vía `xx_from`) reproduce E0 bit a bit, así que la divergencia viene del ruido.
- Los adjuntos de GA (360 d) y GB (180 d) crecen exponencialmente hacia atrás: max |∂J/∂θ₀| = 8×10⁻³ a 30 d, 2×10¹¹ a 180 d y 9×10⁴³ a 360 d. El de 90 días (D90) es estable (2×10⁻³).
- Leith no es la causa (D90_noleith). En D90, ±10⁻³ °C en una celda profunda cambia el J de 90 días en −6.7×10⁻⁴ TW con los dos signos.
- Clasificación: **(c) elección numérica** (la configuración de Dettling a 10 km con Leith no es laminar), con un **error nuestro (a)** en la premisa del spec. El acortamiento previsto (180 días) no basta.

**Lo congelado que dependía de GA/GB se reporta tal cual:**
- el adjunto determinista de ventana larga **falla**;
- V1, V2 y V3 sobre la ventana de GA no se evalúan con GA;
- H5a, H5c y la predicción de H5b no se evalúan con GA/GB.

Los umbrales no cambian. Todo lo siguiente es **post hoc** y se etiqueta así en el README.

**A. Determinista corto: ventana de 30 días desde G0 año 30 (la de GM30, donde el sistema todavía es lineal: dispersión de 3×10⁻⁸ TW).**
- **V1d:** `grdchk` en 10 puntos elegidos con `validate.top_points` sobre GM30:
  - 4 de τx: los 3 mayores separados, más uno en la fila del quiebre (y ≈ 125 km);
  - 3 de T₀ y 3 de κᵥ.

  Se compara contra diferencias finitas y contra el TLM, con los eps de PD1b. El criterio es < 1 % en todos los puntos.
- **V2d y V3d:** los tres patrones de GM30 (τx, T₀ y κᵥ), normalizados a su máximo y multiplicados por las amplitudes del PD0, con corridas directas ±Δu de 30 días. Los criterios son los congelados (0.1 y 30 %).

**B. Adjunto de ensamble (Lea et al., 2000): EA21–EA30, ventanas de 90 días desde los pickups anuales de G0 de los años 21 a 30.**
- **Sensibilidad:** la media de los 10 miembros, con su error estándar (desviación / √10). Los tiempos de anticipación de los registros 1, 2 y 3 son 60, 30 y 0 días (`adjsens.lead_time_days`).
- **Advertencia:** sin `pkg/cal`, los registros son periódicos sobre la ventana, así que el registro 3 también recibe la sensibilidad de los primeros 15 días.
- **V2e y V3e:**
  - el patrón de τx es la media de EA sumada sobre los registros, constante en el tiempo, con máximo 0.01 N m⁻²;
  - corridas directas +Δu y −Δu desde los 10 pickups, de 90 días;
  - Δ± = J± − J₀, donde J₀ es el fc del miembro de EA (la misma trayectoria);
  - V2e = |⟨Δ⁺ + Δ⁻⟩| / |⟨Δ⁺ − Δ⁻⟩| < 0.1;
  - V3e: |ΔJ_pred − ⟨Δ⁺ − Δ⁻⟩/2| / |⟨Δ⁺ − Δ⁻⟩/2| < 0.30, y además la respuesta directa debe superar dos veces su error estándar. Si no lo supera, se reporta "no detectable".
- **H5a:** el criterio congelado, con la sensibilidad media de EA sumada sobre los registros, y las escalas de PD1b (37.5 km, con robustez a 18.75 y 75 km).
  - Se añade, sin cambiar el veredicto, si la diferencia |ΔJ(ii)| − |ΔJ(i)| supera dos veces su error estándar, emparejando por miembro.
  - La confirmación directa son las corridas ±Δu de (i) y (ii) a 37.5 km con 5 miembros (años 26 a 30). Se reporta, pero no cambia el veredicto.
- **H5c:** se usa |Σ_y ∂J/∂τx| medio de EA en los tiempos de anticipación 0, 30 y 60 días. El tiempo de caída a la mitad, t½, se obtiene interpolando linealmente.
  - Cumple si 30 ≤ t½ ≤ 180 días.
  - Si no ha caído a la mitad a los 60 días, el resultado es **"no evaluable más allá de 60 días"**, que no cuenta ni como fallo ni como confirmación.
- **H5b:** la predicción del grueso es un ensamble directo del grueso.
  - Diseño: el patrón de FP, constante, con ±Δu, 180 días, J en los últimos 30, desde los 10 pickups. La predicción es ⟨Δ⁺ − Δ⁻⟩/2.
  - Ya no viene de GB.
  - El patrón de FP es el de mayor sensibilidad de EA (el mismo de V3e), con máximo 0.01 N m⁻² y constante.
  - Los criterios del fino (detectable, signo, factor 2) no cambian.

### Corrección del 2026-11-16, 23:40 (error nuestro, (a); los criterios no cambian)

**Qué pasaba.** Sin `periodicExternalForcing`, `fu` se lee una sola vez, pero `CTRL_MAP_FORCING` le suma `xx_fu` en cada paso (`pkg/ctrl/ctrl_map_forcing.F:119`, `model/src/forward_step.F:572`). Solo el forzamiento periódico repone `fu` desde el archivo (`model/src/external_fields_load.F:268`). El control de viento se acumulaba paso a paso.

**Cómo se detectó.** La primera corrida V3E (0.01 N m⁻²) explotó al día 7.

**Qué afecta:**
- Todas las sensibilidades a τx anteriores (GS, GSU, GM30, GA, GB y la primera tanda de EA) eran respecto de un incremento por paso. En GM30 eran ~690 veces mayores.
- No se ven afectados θ₀, κᵥ ni el diagnóstico de caos: E0–E4 corrieron con controles nulos, y E0z coincide con el valor anterior a 10⁻¹⁰.

**Arreglo** (`scripts/run_adjoint.py`): forzamiento periódico con un registro por ciclo (`externForcingPeriod` = `externForcingCycle` = ventana), recargado en cada paso.

**Resultados con el arreglo:**
- GSU repetido: τx pasa en las tres filas (j18: 4.9×10⁻⁶; j10: 7.6×10⁻⁶; j30: 6.7×10⁻⁴).
- El fallo de τx en GS (3.1 %) era, probablemente, la no linealidad de la acumulación.

**Archivo.** Las corridas viejas se conservan como `*_accum`. EA, GM30, V1d, V2d/V3d, V2e/V3e, H5a, H5b y H5c se rehacen con el arreglo.

### Corrección del 2026-11-17, 08:00 (error nuestro, (a); los criterios no cambian)

**Qué pasaba.** Ni el spin-up grueso G0 ni el fino F0/F0b recibieron una perturbación que rompiera la simetría zonal: batimetría, estado inicial y forzamiento son uniformes en x.
- G0 y todos los estados de EA, V3E, H5A y H5BC son exactamente 2-D: max|θ − media zonal| = 2.2×10⁻¹⁶.
- F0b también lo es (3.9×10⁻¹¹): el spin-up "con remolinos" de 12 meses no tuvo remolinos. El proxy de energía cinética total no podía detectarlo, porque lo domina el chorro medio.
- La revisión final de la rama lo encontró.

**Consecuencias.**
- **Grueso.** Se declara como modelo GM 2-D *deliberado*: su adjunto es el de la parametrización, no el de los remolinos.
  - El "caos" de E0–E4 y el crecimiento de los adjuntos GA y GB son, sobre todo, el crecimiento de modos no simétricos sembrados por el ruido de θ₀ o por el redondeo (pliegue e ≈ 5 d), no caos de un equilibrio 3-D.
  - EA, V1–V3, H5a, H5c y la predicción de H5b no cambian. Se reinterpretan como resultados de un estado zonalmente simétrico.
  - Los 10 puntos de V1d son 4 gradientes distintos: las copias en x tienen el mismo gradiente.
- **Fino.** FC1–3 y FP1–3 muestrearon el inicio de los remolinos (J mensual de −0.8 a 25 GW).
  - Su H5b y su V5 no son válidos. Se conservan como `F[CP][1-3]_transient` y no se usan para ningún veredicto.

**Rediseño del fino (fijado antes de correr).**
- **F1.** 365 días desde el pickup de F0b, con ruido de θ ≤ 10⁻³ °C (miembro 0, misma semilla que FC), en una sola etapa: `pChkptFreq` = 182.5 d no admite etapas de 90 días (ajuste técnico de las 07:25, antes de cualquier resultado).
- **EKE.** EKE(mes) = ½ ⟨UVELSQ − ū² + VVELSQ − v̄²⟩, con ū y v̄ las medias zonales de UVEL y VVEL del mes. ⟨·⟩ es la media de volumen sobre las celdas mojadas.
- **Equilibrio.** Se acepta si la EKE media de los meses 10–12 difiere menos de 10 % de la de los meses 7–9.
  - Si no se cumple, F1 se extiende una vez 180 días (F1b) y se repite la prueba con los meses −6…−4 y −3…−1.
  - Si tampoco se cumple, se sigue de todos modos y se declara.
- **FC1–3 y FP1–3.** Mismo diseño que antes: 180 días, J en los últimos 30, ruido por miembro y el mismo patrón de FP, pero desde el final de F1 (o de F1b).
- **Criterios.** H5b (detectable, signo, factor 2) y V5 (> 50 %) no cambian.

**Resultado de F1 (09:13) y extensión F1b.**
- EKE mensual de F1 (m² s⁻²): 1.1×10⁻², 3.5×10⁻², 7.8×10⁻³, 2.5×10⁻³, 1.2×10⁻³, 9.8×10⁻⁴, 7.7×10⁻⁴, 7.4×10⁻⁴, 7.8×10⁻⁴, 7.0×10⁻⁴, 6.8×10⁻⁴, 6.0×10⁻⁴.
- Meses 10–12 contra 7–9: −13.4 %. **No cumple** la regla del 10 %.
- Se aplica la extensión prevista. F1b dura **182.5 días** y no 180, sin ruido nuevo, para que su final caiga en un múltiplo de `pChkptFreq` y escriba el pickup. Es un ajuste técnico, fijado antes de correr F1b.
- La prueba de F1b compara sus meses 1–3 con sus meses 4–6. FC/FP parten del final de F1b (iteración 328500).

**Resultado de F1b (10:08).**
- EKE mensual de F1b (m² s⁻²): 6.9×10⁻⁴, 5.7×10⁻⁴, 5.2×10⁻⁴, 4.7×10⁻⁴, 4.3×10⁻⁴, 4.9×10⁻⁴.
- Meses 4–6 contra 1–3: −21.6 %. **Tampoco cumple.** Según la regla, se sigue de todos modos: FC/FP parten del final de F1b y la deriva se declara en el README.
- Lectura: tras la ráfaga inicial de inestabilidad (mes 2 de F1), la EKE decae lentamente y sigue decayendo 18 meses después; el modelo fino no ha llegado a un estado estadísticamente estacionario de remolinos. ΔJ de H5b se mide entre miembros pareados que comparten esa deriva, que se cancela en la diferencia en la medida en que es común a los dos brazos.
