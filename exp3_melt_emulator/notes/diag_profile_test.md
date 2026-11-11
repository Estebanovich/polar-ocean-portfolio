# Diagnóstico del pliegue `profile_test` (M7 en `ocean0_4km`) — POST HOC

Fecha: 2026-11-11. **Análisis post hoc**, no preregistrado. No cambia ningún resultado del protocolo: `cv_results.json`, `final_test.json`, `protocol_config.json`, `hypotheses.json` y el código de los métodos quedan intactos. No se corrió `run_protocol.py`.

- Script: `analysis/profile_test_diag.py`.
- Números: `figures/profile_test_diag.json`.
- Arreglos grandes: `/Volumes/LaCie/portfolio_runs/exp3/diag/profile_test_M7_decomposition.npz` (no versionado).

## El hallazgo que hay que explicar

En el pliegue `profile_test` todos los métodos se entrenan con las 25 corridas y se prueban en `ocean0_4km`, el Ocean0 de ISOMIP+: perfil lineal de superficie a fondo, d = 720 m, nunca visto en el entrenamiento. Resultados de la variante principal (unidades `main__profile_test__*.json`):

| Método | `rmse_map` (m a⁻¹) | `integrated_rel_error` | dispersión |
|---|---|---|---|
| M0 | 5.18 | −33 % | — |
| M1 | 9.96 | +62 % | — |
| M2 | 6.76 | **−4.8 %** | — |
| M3 | 9.79 | +59 % | — |
| M6 | 8.47 | +53 % | 9.3 % |
| M7 (local) | **53.5** | **+413 %** | 25 % |

El −4.8 % de M2 tiene signo negativo; el "+5 %" del encargo era su valor absoluto. En la repetición `drop_n1`, M7 da +438 % con el mismo mecanismo (n = 2.16–2.26).

## Clasificación

**No es error nuestro (a).** Es un comportamiento genuino fuera de distribución de la forma elegida para M7: **(c) elección del método**, con una componente de **(b) física que falta** en el factor local. Está sostenido por las pruebas 1–5 de abajo.

- La forma separable m = K·TF₊ⁿ·g_θ(x_g) aprende g_θ condicionada a la relación TF ↔ (posición respecto de z_c) de la familia de entrenamiento, que tiene d = 300 m fijo.
- En `ocean0_4km` esa relación cambia: con d = 720 m, la rampa empieza en la superficie. El producto multiplica entonces un TF entre 3 y 9 veces mayor por un g grande, que se aprendió para compensar un TF casi nulo.
- No hay (d): la referencia es la corrida MITgcm del Exp. 2, con su perfil verificado (prueba 1).

## Evidencia

### 1. Las entradas de `ocean0_4km` se construyen igual que las del ensamble (descarta (a))

- **Perfiles:** el T/S del dataset coincide exactamente con `profiles.ocean0` (AD16 Ec. 21–22, Tabla 6) y con `profiles.ramp(Tbot = 1, z_c = 360, d = 720)`. La diferencia máxima es 0.0 en ambos casos.
- **z_c = 360 tiene sentido.** Es el centro de la rampa con d = 720 m, que es la definición del generador. Por eso `draft_minus_zc` = −draft − 360 está bien calculado. Lo que cambia es su significado físico: el gradiente de T por metro de `draft_minus_zc` es 300/720 del de las corridas de entrenamiento.
- **Malla:** `draft` y `shelf_mask` son idénticos a los de las corridas de entrenamiento con s = 1.0 (696 celdas).
- **TF:** coincide exactamente con un recálculo independiente, T − T_f(S, p) (AD16 Tabla 4, R2); la diferencia máxima es 0.0.
- **Interpolación de T:** difiere en ≤ 0.018 °C de la evaluación analítica. La diferencia sale solo del recorte en las celdas con calado de menos de 10 m, y se aplica igual a todas las corridas.
- **Unidades:** el melt medio es 10.19 m a⁻¹, del mismo orden que el ensamble.
- **Reproducción:** M2, recalculado en forma cerrada, reproduce el −4.8 % de la unidad.
- **Reajuste de M7:** se reajustó una vez con los valores congelados (`tf_kind = "local"`, el que eligió la unidad; 79 s). Reproduce **exactamente** el melt integrado de los 5 miembros (diferencia relativa 0.0) y los mismos K y n. El producto K·TF₊ⁿ·g reconstruye μ con un error de 6e-14.
- **No hay error de escalado ni de unidades.** En las celdas profundas (`draft_minus_zc` ≥ −50 m), el g aprendido coincide con el g necesario (ver la prueba 4).

### 2. Cada variable, por separado, está dentro del rango de entrenamiento

Para las 10 `FEATURES`, la fracción de celdas de `ocean0_4km` fuera del rango [mín, máx] del entrenamiento es **0 %**:

- TF va de 0.00 a 2.70 °C, frente a −0.04 a 3.83 °C en el entrenamiento, así que **n no se extrapola**;
- ⟨TF⟩ = 1.10, frente a 0.13–1.63;
- `draft_minus_zc` va de −354 a 205 m, frente a −544 a 315 m.

### 3. La combinación (posición, TF) no se vio nunca en el entrenamiento

Por bandas de `draft_minus_zc`, el TF de `ocean0_4km` se compara con el máximo del entrenamiento en la misma banda:

| banda (m) | celdas ref. | TF ref. (mediana) | TF entrenamiento (máx.) | celdas ref. por encima del máx. |
|---|---|---|---|---|
| [−400, −250) | 45 | 0.20 | 0.18 | 56 % |
| [−250, −150) | 369 | 0.85 | 0.29 | **100 %** |
| [−150, −50) | 135 | 1.12 | 1.41 | 8 % |
| [−50, 250) | 147 | 1.7–2.6 | 2.6–3.8 | 0 % |

**405 de 696 celdas (58 %)** tienen un TF mayor que cualquier TF visto a la misma distancia del centro de la termoclina.

- En el entrenamiento (d = 300 m), las celdas con `draft_minus_zc` < −150 m están por encima de la rampa: T = T₀ y TF ≈ 0.1 °C.
- En Ocean0 (d = 720 m), esas mismas celdas ya están dentro del gradiente, con TF de 0.5 a 1.0 °C.

### 4. Descomposición de la predicción: el exceso viene de g_θ × TF₊ⁿ en la mitad superior

| banda (m) | Σ melt real | Σ M7 | exceso | g aprendido (mediana) | g necesario = melt/(K TF₊ⁿ) (mediana) |
|---|---|---|---|---|---|
| [−400, −250) | 140 | 362 | +221 | 8.8 | 3.5 |
| [−250, −150) | 1164 | **25 316** | **+24 152** | **4.2** | **0.14** |
| [−150, −50) | 1255 | 6226 | +4971 | 1.4 | 0.21 |
| [−50, 50) | 1657 | 1672 | +15 | 0.30 | 0.30 |
| [50, 150) | 1833 | 1812 | −21 | 0.25 | 0.22 |
| [150, 250) | 1039 | 976 | −63 | 0.17 | 0.19 |

(Σ en m a⁻¹ × celda.)

- **Exceso total:** 29 275. El **84 %** cae en las celdas de la prueba 3 que superan el TF máximo de su banda. **Todo** el exceso está por encima de −50 m: +29 344 ahí, y −69 en las bandas más profundas.
- **Por debajo del centro de la termoclina, M7 acierta:** el error integrado está entre −6 % y +1 % por banda.
- **Sobre el rango de g:** g_θ en `ocean0_4km` (mediana 2.9, máximo 13.1) está **dentro** de su rango en las celdas de entrenamiento (máximo 15.4). Esto era de esperar, porque x_g está en rango.
- **Por qué g es grande arriba:** en el entrenamiento, el melt de las celdas someras se sostiene con un TF local casi nulo. El calor llega por la pluma desde abajo, un efecto no local que TF₊ⁿ no ve. Para reproducirlo, g_θ aprendió un g grande en función de la posición respecto de z_c.
- **Por qué falla en Ocean0:** ahí el TF local de esas celdas es real y grande. El factor g no sabe que ya no hace falta compensar, y el producto explota: con n ≈ 2.4, (0.85/0.1)^2.4 ≈ 170.

### 5. K y g no están identificados por separado

- Con g = 1, K·TF₊ⁿ por sí solo da +323 % (media de los miembros), porque K vale de 16 a 34, frente a K₀ ≈ 6.3 de M2.
- En profundidad, un g de 0.2 a 0.3 compensa ese K. Solo el producto está identificado, como ya decía el docstring de `nets.py`.
- Por eso no tiene sentido culpar a un solo factor. El fallo está en la **interacción** entre el TF fuera de distribución condicional y el g aprendido.
- **Por qué la dispersión es de 25 %:** los miembros difieren en K, n y g justo en la zona sin datos (integrado con g = 1: de 20 500 a 43 300). Aun así, la dispersión es muy inferior al error de +413 %. Es el patrón que H3c busca, aunque `profile_test` no entra en los 33 casos de H3c.

## Por qué M2 y M6 no fallan igual

- **M2** (K c² |TF| TF) no tiene ningún término condicionado a la posición. Sobreestima en la banda [−250, −150) (1618 frente a 1164, +39 %) y subestima en [−50, 150) (−24 % y −14 %). Los errores se compensan en parte en el integrado (−4.8 %), aunque su `rmse_map` es 6.8. M2 no tiene ningún factor aprendido que dependa de la posición. En esa banda su error es +39 %, frente a +2070 % de M7 (25 316 frente a 1164): el g aprendido de M7 es unas 30 veces el necesario (4.2 frente a 0.14).
- **M6** recibe TF, T, S y ⟨TF⟩ como entradas aditivas de la red. Un TF grande no se amplifica de forma multiplicativa con exponente 2.4, así que su error se queda en +53 %.

## Conclusiones

1. **No se encontró ningún error.** Se descartan, con las pruebas 1–4:
   - las entradas, el TF y el liquidus;
   - las máscaras y las unidades;
   - el escalado y la reproducción.
2. **El fallo de M7 es (c) elección del método, con algo de (b):**
   - La forma separable K·TF₊ⁿ·g(x_g) con TF local no representa el aporte no local de calor de la pluma.
   - g_θ lo absorbe como función de la posición respecto de z_c, válida solo para la forma de perfil de entrenamiento (d = 300 m).
   - Ningún diagnóstico marginal de rango lo detecta: es un fuera de distribución **condicional**, porque cada variable está en rango por separado.
3. **Matiz para el README:** el preregistro excluyó M4 y M5 de `profile_test` porque "no pueden representar d = 720 m". El diagnóstico muestra que M7 tampoco lo representa, aunque formalmente reciba las variables por celda. Su g_θ hereda implícitamente d = 300 m a través de `draft_minus_zc`. M6 también recibe `draft_minus_zc`, pero lo atenúa con las entradas de temperatura.
4. **Nada de esto cambia los resultados ni las hipótesis preregistradas.** Es interpretación post hoc. Una variante de M7 con un forzamiento no local en el término físico, o sin `draft_minus_zc` en x_g, sería un experimento nuevo y post hoc. No se probó aquí.
