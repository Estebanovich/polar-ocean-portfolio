# PD2 — Ensamble a 4 km: H1, H2 y N1–N3

Fecha: 2026-11-09. Cola de producción (PID 23111): el ensamble terminó a las 22:03 CST (`END ens rc=0`). Ocean1 sigue corriendo y no entra aquí.

Definiciones: `pd0_preregistro.md`, incluida la aclaración R22, que se fijó a las 16:41 CST (commit 968410f), antes de evaluar cualquier salida del ensamble. Todo se calcula sobre los últimos 6 meses (registros 19–24) de corridas de 2 años a 4 km, con ΓT_cal,4km = 0.029356.

Código y salidas:
- código: `analysis/ensemble.py`, `fig2_melt_sensitivity.py` y `fig3_meltwater_export.py`;
- tabla: `figures/ensemble_table.json`;
- veredictos: `figures/fig2_h1.json` y `figures/fig3_h2.json`.

## N1–N3 (29 corridas: referencia y 28 `ens_*`)

- **N3:** pasan las 29; no hay huecos.
- **N2:** pasan las 29 (máx. |residuo| = 3.1e-7 < 1 %).
- **N1:** fallan 5, con deriva > 3 %:

  | Corrida | N1 |
  |---|---|
  | `ens_T-0.5_zc250_s1.0` | 3.7 % |
  | `ens_T-0.5_zc550_s1.0` | 6.5 % |
  | `ens_T1.5_zc550_s1.0` | 3.4 % |
  | `ens_T-0.42_zc402_s0.8` | 8.8 % |
  | `ens_T-0.39_zc503_s1.2` | 4.9 % |

  Son, sobre todo, casos fríos o de termoclina profunda, con melt pequeño (0.3–4.5 m a⁻¹), en los que la misma deriva absoluta pesa más en términos relativos. Según la aclaración R22, siguen en el ajuste y se señalan. Ninguna corrida se alargó.
- **Verificación cruzada:** `ensemble.collect` reproduce exactamente los N1 y N2 de `ensemble_n_check.py`, el paso `ncheck` de la cola (diferencia máxima 0.0).

## H1 — forma cuadrática (criterio: ⟨m⟩ de toda la plataforma, R22)

T_f,GL = −2.3341 °C (Sbot = 34.7; base mojada más profunda 565.2 m).

| z_c (m) | R² cuadrática | n (ley de potencia) | k | Veredicto | Corridas con N1 fallido |
|---|---|---|---|---|---|
| 250 | 0.9998 ✔ | 2.43 ✘ | 1.04 | **FALLA** | T −0.5 |
| 400 | 0.9995 ✔ | 2.78 ✘ | 0.29 | **FALLA** | — |
| 550 | 1.0000 ✔ | 2.70 ✘ | 0.069 | **FALLA** | T −0.5, T 1.5 |

**H1 falla en los tres z_c.** La forma cuadrática ajusta muy bien, pero el exponente queda por encima del rango preregistrado [1.3, 2.0]: el melt crece más rápido que (Tbot − T_f,GL)².

**Secundario, no criterio (R22):** con el melt profundo (base < −300 m) sale n = 1.99 (z_c 250, dentro del rango), 2.49 (400) y 2.62 (550), con R² ≥ 0.9998. No cambia el veredicto.

**Robustez (*post hoc*):** se calcularon tres variantes, todas sobre el melt de toda la plataforma.

| z_c (m) | n sin un punto (5 ajustes) | n sin las corridas con N1 fallido | pendientes locales d ln m / d ln ΔT |
|---|---|---|---|
| 250 | 2.40–2.48 | 2.40 | 2.19–2.55 |
| 400 | 2.75–2.82 | 2.78 | 2.68–2.99 |
| 550 | 2.60–2.76 | 2.67 | 2.43–2.88 |

En ninguna variante n baja de 2.19. El fallo no depende de un punto aislado ni de las corridas que fallan N1.

**Lectura (hipótesis, sin verificar):** en la familia "rampa", la base del hielo poco profunda está en la termoclina o por encima de ella. Al subir Tbot, crecen a la vez la temperatura del agua que llega a la base y el área de base expuesta a agua cálida. Una variable explicativa que es Tbot, y no la temperatura en la base, puede dar un exponente efectivo mayor que 2. Clasificación provisional: (c), por elección del diseño (la variable Tbot de la rampa), no (a). Para probarlo haría falta un ajuste contra la temperatura media en la base, que sería *post hoc*.

## H2 — profundidad de salida del agua de deshielo

| z_c (m) | D_out (m) para Tbot = −0.5, 0.0, 0.5, 1.0, 1.5 | Spearman ρ | Veredicto |
|---|---|---|---|
| 250 | 250.6, 199.3, 141.8, 127.9, 117.6 | −1.0 | **PASA** |
| 400 | 374.2, 322.4, 181.2, 133.4, 102.6 | −1.0 | **PASA** |
| 550 | 460.4, 334.5, 113.7, 97.4, 90.7 | −1.0 | **PASA** |

**H2 pasa en los tres z_c:** D_out baja de forma monótona con Tbot. Las corridas de prueba (Tbot = 2.0: 113, 99 y 97 m) siguen la tendencia, pero no entran en H2. Ocean0: D_out = 119 m.

## Figuras

- `figures/fig2_melt_sensitivity.{png,pdf}`, con tres paneles:
  - (a) ⟨m⟩ contra Tbot, con el ajuste cuadrático, el punto de Ocean0 (10.19 m a⁻¹) y Tbot = 2.0 como prueba;
  - (b) mapa de calor;
  - (c) efecto de s, solo descriptivo.
- `figures/fig3_meltwater_export.{png,pdf}`: (a) secciones en x_front para Tbot −0.5 y 1.5 con z_c = 400; (b) D_out contra Tbot; (c) mapa de ∫C dz para Tbot 1.5 y z_c 400.
- `figures/fig3_anim_ocean0_4km_r21.gif`: 0.71 MB, 24 meses.

ψ_bt no se usa en ningún resultado: a 4 km su máximo en la cavidad es 2.2 veces el de 2 km.
