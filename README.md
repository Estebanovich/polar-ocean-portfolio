# Polar ocean modelling portfolio — Esteban Cruz Isidro

Idealized, reproducible MITgcm experiments on polar ocean processes (ice–ocean interaction, buoyant plumes, fjord and shelf dynamics). Each experiment pairs a model configuration with theory, validation tests and analysis scripts. Its success criteria are fixed before the runs, and failed criteria are reported as failed. This repository supports my PhD applications.

## Experiments

1. **[Subglacial discharge plumes at a tidewater glacier front](exp1_fjord_plumes/README.md)**. Three configurations:
   - a 2-D line plume validated against plume theory with the same melt law;
   - two interacting 3-D point plumes in a stratified ambient, compared with laboratory experiments;
   - a rotating fjord forced separately by a plume, wind and the shelf.

How to run the model, the tests and the analysis: [RUNNING.md](RUNNING.md).

## Previous work

- **[pyOCEAN](https://github.com/Estebanovich/pyOCEAN)**: a regional and coastal ocean model I am writing from scratch in Python (NumPy + Numba). It has z\*, σ and HYCOM-style hybrid vertical coordinates and ports of the MOM6, CROCO and CVMix closures (GLS, KPP, Jackson shear mixing, bottom-boundary-layer mixing). It is checked against MOM6, MITgcm and CROCO in standard benchmarks, and against tide gauges, satellite SST and gliders in a regional run of the Mid-Atlantic Bight. Below is the DOME dense overflow (Legg et al. 2006), the kind of flow that carries dense shelf water into the abyss: without shear mixing (left) the plume stays too dense, with it (middle) it entrains like MITgcm (right).

  <img src="https://raw.githubusercontent.com/Estebanovich/pyOCEAN/52c035c/docs/media/dome_overflow.gif" width="90%" alt="DOME overflow bottom dye in pyOCEAN without and with shear mixing, against MITgcm">

- [canyon_mx](https://github.com/Estebanovich/canyon_mx)
- [BayIW_Rectan](https://github.com/Estebanovich/BayIW_Rectan)
- [IWBayMITgcm-postprocessing](https://github.com/Estebanovich/-IWBayMITgcm-postprocessing)


## License

Code is released under the [MIT License](LICENSE). Text and figures are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Satellite imagery in Figure 0 keeps its original terms (NASA Worldview/GIBS; modified Copernicus Sentinel data), as credited in its caption.
