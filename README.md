# Polar ocean modelling portfolio — Esteban Cruz Isidro

Idealized, reproducible MITgcm experiments on polar ocean processes (ice–ocean interaction, buoyant plumes, fjord and shelf dynamics). Each experiment pairs a model configuration with theory, validation tests and analysis scripts. Its success criteria are fixed before the runs, and failed criteria are reported as failed. This repository supports my PhD applications.

## Experiments

1. **[Subglacial discharge plumes at a tidewater glacier front](exp1_fjord_plumes/README.md)**. Three configurations:
   - a 2-D line plume validated against plume theory with the same melt law;
   - two interacting 3-D point plumes in a stratified ambient, compared with laboratory experiments;
   - a rotating fjord forced separately by a plume, wind and the shelf.
2. **[An ISOMIP+ ice-shelf cavity: melt sensitivity to warm deep water and the destination of meltwater](exp2_isomip_cavity/README.md)**. Three parts:
   - the ISOMIP+ protocol, calibrated as prescribed and checked against the published results and the MITgcm-JPL output: melt and barotropic circulation in range, cavity overturning below the pre-registered range;
   - a 28-run ensemble at 4 km varying the deep-water temperature, the thermocline depth and the ice draft: the melt is a near-exact quadratic in temperature but steeper than the pre-registered power law (that hypothesis fails), and the meltwater leaves the cavity at shallower depths as the deep water warms (that hypothesis holds);
   - the ensemble exported as the training and test data of a melt emulator (Experiment 3).
3. **[Emulating ice-shelf basal melt from a small cavity ensemble](exp3_melt_emulator/README.md)**. Eight methods, from a one-parameter quadratic melt law to a deep ensemble of neural networks and a physics hybrid, trained on the Experiment 2 runs and scored only on runs they had not seen, under a pre-registered protocol with a final test at warmer deep water evaluated once:
   - with 25 training runs, most differences among the better methods are not resolved: on the melt map the deep ensemble beats the one-parameter melt parameterisations (except when the thermocline depth is held out) but not a per-cell regression or a Gaussian process on the three forcing parameters, and none of these is resolved from the regression inside the training range. When the thermocline depth is held out, the one-parameter local quadratic law is resolved best;
   - the deep ensemble has a resolved lead on the integrated melt only inside the training range (not when the thermocline depth is held out), and its spread flagged its large errors on the pre-registered out-of-range cases, replicating a published result (DeepMelt-GL) in a controlled case, but not on a profile shape it had not seen;
   - the physics hybrid did not extrapolate better than the deep ensemble (that hypothesis fails) and, in a post hoc diagnosis, failed on a profile shape it had not seen because its geometry network had learned the training profile shape.

How to run the model, the tests and the analysis: [RUNNING.md](RUNNING.md).

## Previous work

- **[pyOCEAN](https://github.com/Estebanovich/pyOCEAN)** ([model manual](https://estebanovich.github.io/pyOCEAN/)): a regional and coastal ocean model I am writing from scratch in Python (NumPy + Numba). It has z\*, σ and HYCOM-style hybrid vertical coordinates and ports of the MOM6, CROCO and CVMix closures (GLS, KPP, Jackson shear mixing, bottom-boundary-layer mixing). It is checked against MOM6, MITgcm and CROCO in standard benchmarks, and against tide gauges, satellite SST and gliders in a regional run of the Mid-Atlantic Bight. Below is the DOME dense overflow (Legg et al. 2006), the kind of flow that carries dense shelf water into the abyss: without shear mixing (left) the plume stays too dense, with it (middle) it entrains like MITgcm (right).

  <img src="https://raw.githubusercontent.com/Estebanovich/pyOCEAN/52c035c/docs/media/dome_overflow.gif" width="90%" alt="DOME overflow bottom dye in pyOCEAN without and with shear mixing, against MITgcm">

- [canyon_mx](https://github.com/Estebanovich/canyon_mx)
- [BayIW_Rectan](https://github.com/Estebanovich/BayIW_Rectan)
- [IWBayMITgcm-postprocessing](https://github.com/Estebanovich/-IWBayMITgcm-postprocessing)


## License

Code is released under the [MIT License](LICENSE). Text and figures are released under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Satellite imagery in Figure 0 keeps its original terms (NASA Worldview/GIBS; modified Copernicus Sentinel data), as credited in its caption.
