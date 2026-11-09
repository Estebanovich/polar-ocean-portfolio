# External data for Exp. 2

Files in this directory are not tracked by git. Download them as described below.

| File | Source | Downloaded | SHA-256 |
|---|---|---|---|
| `Ocean1_input_geom_v1.01.nc` | Cornford, S. L. and Asay-Davis, X. S. (2016): ISOMIP+ ice-shelf topography, GFZ Data Services, https://doi.org/10.5880/PIK.2016.002 — file from https://datapub.gfz.de/download/10.5880.PIK.2016.002/ | 2026-11-08 | `765735a954e8b89e8845a2a11fc205e67514def3ab9f6319ab5be79a20e12be9` |

`Ocean1_input_geom_v1.01.nc` is the geometry for both Ocean0 and Ocean1 (Asay-Davis et al. 2016, §3.1.1). It holds a 1 km grid of 480 × 80 cells with centres at x = 320.5–799.5 km and y = 0.5–79.5 km. The variables, all dimensioned (y, x), are `upperSurface`, `lowerSurface`, `bedrockTopography`, `floatingMask`, `groundedMask` and `openOceanMask`.

```bash
curl -L -o exp2_isomip_cavity/data_external/Ocean1_input_geom_v1.01.nc \
  https://datapub.gfz.de/download/10.5880.PIK.2016.002/Ocean1_input_geom_v1.01.nc
shasum -a 256 exp2_isomip_cavity/data_external/Ocean1_input_geom_v1.01.nc
```
