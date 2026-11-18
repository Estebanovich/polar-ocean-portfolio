# Exp. 6: puerta del día 1 (D2), cronometría y decisiones (17–18 nov 2026)

Todas las pruebas son de 30 días, 10 días o 1 año, corridas en la laptop (núcleos de rendimiento, hilos en 1). Las salidas están en `/Volumes/LaCie/portfolio_runs/exp6/_timing/`.

## Mediciones

| Prueba | Config | np | Paso | Duración | Pared | Resultado |
|---|---|---|---|---|---|---|
| T10_30d | g50 | 1 | 1000 s | 30 d | 68 s | OK; 75 % del tiempo en termodinámica (6 trazadores con advección 77 y GM-Redi) |
| g50_np5 | g50 | 5 | 1000 s | 30 d | 35 s | OK; solo ×1.9 con 5 rangos (mosaicos de 8 filas con halo de 4) |
| dt2000 | g50 | 1 | 2000 s | 1 año | 414 s | OK; CFL máximos: u 0.021, v 0.017, w 0.081 |
| dt3000 | g50 | 1 | 3000 s | 1 año | — | **NaN** (inestable) |
| f10_dt1000 | f10 | 5 | 1000 s | 10 d | 115 s | OK, desde el reposo (todavía sin remolinos) |
| f10_dt500 | f10 | 5 | 500 s | 10 d | 222 s | OK |
| g0acc | g50p5 | 5 | 2000 s mom., 8000 s trazadores | 8 años de trazadores | 315 s | OK; DIC y PO4 positivos; CFL vertical de trazadores 0.71 |

## Decisiones (antes de seguir; clasificación (c), elección numérica)

1. **`g50` con `deltaT` = 2000 s,** no los 1000 s del tutorial.
   - Los CFL advectivos quedan por debajo de 0.1. A 3000 s el modelo diverge, así que el margen es de 1.5×.
   - Ahorra la mitad del costo de todas las corridas gruesas.
2. **El spin-up G0 está acelerado** (Bryan 1984 [verificar]): `deltaTtracer` = `deltaTClock` = 8000 s y `deltaTMom` = 2000 s, con 5 rangos (`g50p5`).
   - Los trazadores avanzan 4 años por cada año de dinámica. Es válido para llegar al equilibrio, no para el ciclo estacional, así que después vienen **48 años síncronos (G0s)** antes de GC y GT.
   - V6a se evalúa sobre los años de trazador de G0 y se reporta también la deriva en G0s.
   - G0 usa múltiplos de 4 años para que la fase estacional no salte al pasar a G0s.
   - **Regla registrada:** si G0 produce NaN o PO4 < 0, se repite con factor 2 (`deltaTtracer` = 4000 s).
3. **`f10` con `deltaT` = 1000 s** (provisional).
   - **Regla registrada:** si F0 diverge con remolinos desarrollados, se baja a 500 s y se aplica la regla de recorte.
4. **El grueso corre con 1 rango,** con hasta 5 casos a la vez. G0, G0s y GTa corren con 5 rangos, porque en esos casos no hay paralelismo entre casos.

## Proyección (horas de pared)

| Bloque | Años de modelo | Costo | Horas |
|---|---|---|---|
| G0 (acelerado, np 5) | 1000 de trazador | 0.66 min/año | 11 |
| G0s (np 5) | 48 | 3.6 min/año | 3 |
| GC y GT en paralelo, con las 8 ramas en los núcleos libres | 251 + 8 × 25 | 6.9 min/año por núcleo | 29 |
| GTa + 2 ramas (np 5, después de GT y GC) | 251 + 50 | 3.6 min/año | 18 |
| **Grueso total** | | | **≈ 61** |
| F0_1900, F0_2060 (≤ 8 años c/u) | 16 | 70 min/año | 19 |
| FT: 2 años de rama × 3 miembros × (c, w) × 15 años | 180 | 70 min/año | 210 |
| **Fino total** | | | **≈ 229** |

- **La regla de recorte no se activa:** el fino suma unas 229 h, por debajo de 250 h, así que se corren 3 miembros.
- Si el fino tuviera que bajar a 500 s, serían 458 h. Con 2 miembros serían 313 h, que también pasan del tope. Quedaría entonces una sola pareja c/w por año: 76 años, 175 h, como ilustración sin hipótesis.
- **Calendario**, si D1 llega el 18 nov:
  - grueso del 18 al 21 nov;
  - fino del 21 nov al ~2 dic;
  - README del 3 al 8 dic;
  - correo el 9–10 dic.

## Evento del 18 nov: G0 a 8000 s diverge; se aplica la regla registrada (factor 2)

- **Primera etapa de G0** (`deltaTtracer` = 8000 s, 100 años de trazador):
  - del año de trazador 1 al 8 es estable: u_max ≈ 0.54 m/s y CFL vertical < 0.75;
  - desde el año 9, el CFL vertical de los trazadores pasa de 1 (`advcfl_W_hf_max` llega a 2.4);
  - en el año 19 el DIC da NaN y después el PO4 se vuelve negativo;
  - en el año 31 la dinámica también da NaN.
- **Causa:** w crece con el spin-up de la circulación profunda.
  - En T10, síncrono, w_max sube de 0.3 a 1.3 cm/s en 10 años, cerca de la pared sur a unos 1000 m.
  - Con celdas de ~100 m y 8000 s, el CFL queda por encima de 1, y el esquema 7 (OS7MP) no es estable así.
  - El namelist es idéntico al del tutorial salvo `deltaT`.
- **Clasificación:** (c), elección numérica, decidida el día 1 a partir del CFL de un solo año.
- **Acción:** la regla registrada de antemano.
  - `deltaTtracer` = `deltaTClock` = 4000 s, con 7884 pasos por año de trazador.
  - G0 se reinicia desde cero. G0 dura el doble: unas 22 h en lugar de 11.
- **Si 4000 s también diverge,** G0 pasa a ser síncrono (2000 s). Se rotula *post hoc*. (Superado por el evento siguiente: 2000 s resultó inestable.)

## Evento del 18 nov (2): G0 a 4000 s también diverge; la causa es el paso de momento de 2000 s

- **G0 a 4000 s** (`deltaTMom` = 2000 s): eta da NaN hacia el año de trazador 10 (año ~5 de dinámica), con ke_max = 4.7. Archivado en `_failed_G0_dt4000`.
- **Diagnóstico** (np 1, salidas en `/Volumes/LaCie/portfolio_runs/exp6/_diag/`):

| Prueba | Configuración | Resultado |
|---|---|---|
| D1 | T10 continuado del año 10 al 20, síncrono a 2000 s | ke_max ya vale 1.14 en el año 10 (a 1000 s es 0.166) y sube a 4.1 en el año 11; luego NaN. **2000 s es inestable aun en síncrono** |
| D2 | síncrono a 1000 s, 10 años desde el reposo | estable: ke_max 0.166 en el año 10; CFL máx. 0.011; PO4 mín. 3.9×10⁻⁴ |
| D3 + D3b | momento 1000 s, trazadores y reloj 4000 s, 60 años de trazador (15 de dinámica) | estable: ke_max 0.169 al final; CFL máx. 0.044; PO4 mín. 4.7×10⁻⁴; 210 s de pared por año de trazador con np 1 |

- **Causa:** la decisión 1 (`g50` a 2000 s) se tomó con el CFL de un solo año. La inestabilidad aparece al desarrollarse la circulación (picos de u de 1–2 m/s y w de ~1.3 cm/s junto a la pared sur, a ~1000 m), después del año 10. Con la regla del factor 2 se bajó el paso equivocado (el de trazadores).
- **Clasificación:** (c), elección numérica. También cuenta como error nuestro de procedimiento (a): la prueba del día 1 era demasiado corta para el paso de momento.
- **Acción** (*post hoc*, antes de congelar el PD0):
  - `g50` y `g50p5` a `deltaT` = 1000 s, el valor del tutorial. Son 31 536 pasos por año.
  - G0 con `deltaTMom` = `deltaTFreeSurf` = 1000 s y `deltaTtracer` = `deltaTClock` = 4000 s (×4), con 7884 pasos por año de trazador. Los namelists de G0 son idénticos a los de D3, salvo `nTimeSteps`.
  - GC, GT y sus ramas no se habían corrido; no hay que repetir nada más.

## Proyección revisada (18 nov, con 1000 s)

Costos medidos: 13.9 min por año síncrono con np 1 (D2, corriendo junto a D3b); 3.5 min por año de trazador acelerado con np 1 (D3). Para np 5 se usa el ×1.9 del día 1.

| Bloque | Años de modelo | Costo | Horas |
|---|---|---|---|
| G0 (acelerado, np 5) | 1000 de trazador | ~1.85 min/año | ~31 |
| G0s (np 5) | 48 | ~7.3 min/año | ~6 |
| GC y GT en paralelo, con las 8 ramas en los núcleos libres | 251 + 8 × 25 | 13.9 min/año por núcleo | ~58 |
| GTa (np 5) y luego sus 2 ramas (np 1, en paralelo) | 251 + 2 × 25 | 7.3 y 13.9 min/año | ~36 |
| **Grueso total** | | | **≈ 131** (antes 61) |
| **Fino total** (sin cambio; `f10` ya estaba en 1000 s) | | | **≈ 229** |

- La regla de recorte del fino no se activa (229 h < 250 h).
- Grueso y fino suman ~360 h, unos 15 días de máquina. Con G0 arrancando el 18 nov, el grueso termina hacia el 24 nov y el fino hacia el 4 dic. Quedan 4–5 días para el README antes del correo del 9–10 dic, sin margen para fallos.

## Decisión del 18 nov: se quita GTa

- **Velocidad real de G0:** la primera etapa (100 años de trazador) tardó 6965 s, es decir 1.16 min por año de trazador con np 5, contra los ~1.85 estimados. G0 completo baja de ~31 h a ~20 h.
- **Recorte:** Esteban aprobó quitar GTa y sus dos ramas (la sensibilidad a la esponja con C_ant; V6e era descriptiva y no entra en ninguna hipótesis). Ahorra ~36 h.
- **Proyección resultante:** grueso ≈ 84 h (G0 ~20, G0s ~6, GC y GT con sus ramas ~58) y fino ≈ 229 h; ~313 h en total. El grueso termina hacia el 22 nov y el fino hacia el 2 dic.
- **Limitación que se reporta:** sin GTa, la esponja norte no lleva carbono antropogénico, así que el C_ant del canal entra solo por la superficie.
