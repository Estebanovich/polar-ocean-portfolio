# PD1 — Puente con la tesis: `bay_thesis`, V1 y R1k (Task 3)

(borrador; se completa con los resultados)

## R1k frente a R1

<!-- R1k-R1 -->
| archivo | grupo | clave | R1 | R1k | motivo |
|---|---|---|---|---|---|
| `data.pkg` | PACKAGES | useKPP | — | `.TRUE.` | activa KPP |
| `data` | PARM01 | viscAr | `1.E-5` | `1.E-4` | viscosidad de fondo de KPP (PD0 §6, L94 ec. 29) |
| `data` | PARM01 | ivdc_kappa | `1.` | — | KPP no admite otra convección (`pkg/kpp/kpp_check.F:128-133`) |
| `data` | PARM01 | implicitViscosity | — | `.TRUE.` | KPP la exige (`pkg/kpp/kpp_check.F:141-145`) |
| `data.free` | PARM01 | viscAr | `1.E-5` | `1.E-4` | ídem, etapa libre |
| `data.free` | PARM01 | ivdc_kappa | `1.` | — | ídem |
| `data.free` | PARM01 | implicitViscosity | — | `.TRUE.` | ídem |
| `data.free` | PARM03 | nTimeSteps | `7200` | `21600` | R1k corre 240 h (60 h forzadas + 180 h libres) |
| `data.kpp` | KPP_PARM01 | * | — | present | grupo vacío: todos los parámetros de KPP por omisión (PD0 §6) |
| `data.diagnostics` | DIAGNOSTICS_LIST | fields(1,4) | — | `'KPPhbl'` | salida: capa límite de KPP |
| `data.diagnostics` | DIAGNOSTICS_LIST | fileName(4) | — | `'kpp2d'` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | frequency(4) | — | `-900.` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | timePhase(4) | — | `0.` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | fields(1:2,5) | — | `'KPPviscA','KPPdiffT'` | salida: viscosidad y difusividad de KPP |
| `data.diagnostics` | DIAGNOSTICS_LIST | fileName(5) | — | `'kpp3d'` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | frequency(5) | — | `-10800.` | salida |
| `data.diagnostics` | DIAGNOSTICS_LIST | timePhase(5) | — | `0.` | salida |
<!-- /R1k-R1 -->
