# External data for Exp. 2

Files in this directory are not tracked by git. Download them as described below.

| File | Source | Downloaded | SHA-256 |
|---|---|---|---|
| `Ocean1_input_geom_v1.01.nc` | Cornford, S. L. and Asay-Davis, X. S. (2016): ISOMIP+ ice-shelf topography, GFZ Data Services, https://doi.org/10.5880/PIK.2016.002 — file from https://datapub.gfz.de/download/10.5880.PIK.2016.002/ | 2026-11-08 | `765735a954e8b89e8845a2a11fc205e67514def3ab9f6319ab5be79a20e12be9` |
| `osf_jpl/MITgcm_JPL_Summary.pdf` | Nakayama, Y. and Menemenlis, D.: ISOMIP+ Ocean0–2 COM MITgcm-JPL results, summary, OSF project ISOMIP+/MISOMIP1 model output, https://doi.org/10.17605/OSF.IO/QNA89 | 2026-11-09 | `908834b1c723c60f606072c294d198a1bf376e0713e4d94974640f0cb2c3b836` |
| `osf_jpl/Ocean0_yosi_COM_MITgcm_JPL_new.nc` | same OSF project (MITgcm-JPL Ocean0 COM output, 12 monthly records) | 2026-11-09 | `cd9e50bd93533fa29dc34546886abbc4522a4327807b9581de11c05e33477e9a` |
| `osf_jpl/Ocean1_yosi_COM_MITgcm_JPL_new.nc` | same OSF project (MITgcm-JPL Ocean1 COM output) | 2026-11-09 | `8a0427ff231d823389d7db4f4b035863350db3c7b489745b2bb5d708d523a58d` |

`Ocean1_input_geom_v1.01.nc` is the geometry for both Ocean0 and Ocean1 (Asay-Davis et al. 2016, §3.1.1). It holds a 1 km grid of 480 × 80 cells with centres at x = 320.5–799.5 km and y = 0.5–79.5 km. The variables, all dimensioned (y, x), are `upperSurface`, `lowerSurface`, `bedrockTopography`, `floatingMask`, `groundedMask` and `openOceanMask`.

```bash
curl -L -o exp2_isomip_cavity/data_external/Ocean1_input_geom_v1.01.nc \
  https://datapub.gfz.de/download/10.5880.PIK.2016.002/Ocean1_input_geom_v1.01.nc
shasum -a 256 exp2_isomip_cavity/data_external/Ocean1_input_geom_v1.01.nc
```

The `osf_jpl/` files are the MITgcm-JPL (Nakayama, Menemenlis) ISOMIP+ submission, used only for the post hoc comparison in `notes/pd1_comparacion_jpl.md` (`analysis/jpl_compare.py`). They are git-ignored like every file here except this one; download them from the OSF project (DOI above) into `exp2_isomip_cavity/data_external/osf_jpl/` and check them with `shasum -a 256`.
