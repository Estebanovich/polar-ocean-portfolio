# Polar ocean modelling portfolio — Esteban Cruz Isidro

Idealized, reproducible MITgcm experiments on polar ocean processes (ice–ocean interaction, buoyant plumes, fjord and shelf dynamics). Each experiment pairs a model configuration with theory, validation tests and analysis scripts. Its success criteria are fixed before the runs, and failed criteria are reported as failed. This repository supports my PhD applications.

## Experiments

1. **[Subglacial discharge plumes at a tidewater glacier front](exp1_fjord_plumes/README.md)**. Three configurations:
   - a 2-D line plume validated against plume theory with the same melt law;
   - two interacting 3-D point plumes in a stratified ambient, compared with laboratory experiments;
   - a rotating fjord forced separately by a plume, wind and the shelf.

How to run the model, the tests and the analysis: [RUNNING.md](RUNNING.md).

## Previous work

- [pyOCEAN](https://github.com/Estebanovich/pyOCEAN)
- [canyon_mx](https://github.com/Estebanovich/canyon_mx)
- [BayIW_Rectan](https://github.com/Estebanovich/BayIW_Rectan)
- [IWBayMITgcm-postprocessing](https://github.com/Estebanovich/-IWBayMITgcm-postprocessing)
- [HYCOM-IOA](https://github.com/Estebanovich/HYCOM-IOA)

## License

Code is released under the [MIT License](LICENSE). Text and figures are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Satellite imagery in Figure 0 keeps its original terms (NASA Worldview/GIBS; modified Copernicus Sentinel data), as credited in its caption.
