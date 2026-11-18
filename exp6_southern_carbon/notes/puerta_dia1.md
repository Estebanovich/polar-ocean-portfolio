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
