# Running the experiments

A portfolio of idealized, reproducible MITgcm polar-ocean experiments. Each experiment pairs a model configuration with theory, validation tests and analysis scripts. Experiment 1 (`exp1_fjord_plumes/`) studies subglacial plumes in a fjord; Experiment 2 (`exp2_isomip_cavity/`) an ISOMIP+ ice-shelf cavity.

## Commands

All Python runs in one conda environment (numpy, scipy, matplotlib, pyyaml, pytest):

```bash
PY=/path/to/env/bin/python
$PY -m pytest                                                                     # all tests
$PY -m pytest exp1_fjord_plumes/tests/test_plume.py::<test_name>                  # single test

# Build MITgcm with <config>/code as -mods (the MITgcm source is never modified)
scripts/build_mitgcm.sh exp1_fjord_plumes/config/<cfg> <build_dir>/<cfg> [--mpi]
# Run one case: copies binary + input, writes run_status.json (OK | FAILED | RUNNING)
scripts/run_mitgcm.sh <build_dir> <input_dir> <run_dir> <np>

# Campaign: generate inputs, build once, run every matching case
$PY exp1_fjord_plumes/scripts/run_campaign.py --config a_line2d --cases '<glob>' --np N \
    [--end-time S] [--run-suffix _short2h] [--set KEY=VALUE] [--skip-ok] [--dry-run]

# Everything (figures, metrics and README) from existing run output
bash exp1_fjord_plumes/reproduce.sh --quick
bash exp2_isomip_cavity/reproduce.sh --quick      # also writes the Exp. 3 dataset (outside the repo)
```

Experiment 2 uses the same scripts with its own configurations (`ocean0_2km`, MPI np 5; `ocean0_4km`, serial):

```bash
scripts/build_mitgcm.sh exp2_isomip_cavity/config/ocean0_4km /Volumes/ESTEBAN_MIT/portfolio_builds/exp2/ocean0_4km
$PY exp2_isomip_cavity/scripts/run_campaign.py --config ocean0_4km --cases 'ens_*' --np 1 \
    --max-concurrent 5 --run-suffix _r21 --skip-ok          # Ocean1: --stages 4 (pickup restarts)
# The campaigns as actually run (each wrapped in caffeinate so an idle Mac does not pause them):
caffeinate -is bash exp2_isomip_cavity/scripts/queue_calibration.sh   # K1 calibration + verification
caffeinate -is bash exp2_isomip_cavity/scripts/queue_production.sh    # 28-run ensemble, N1-N3, Ocean1
caffeinate -is bash exp2_isomip_cavity/scripts/queue_r24.sh           # post hoc R24 runs
bash exp2_isomip_cavity/reproduce.sh --quick [--gifs]
```

## Pipeline (Experiment 1)

The pipeline runs `cases.yaml` → `gendata` → MITgcm run → `analysis`; `theory/` provides the reference values for validation.

- **`cases.yaml`** is the case registry. `configs:` holds the base parameters of each MITgcm configuration, and `cases:` override them. The case-name prefix selects the configuration:
  - `a_line2d`: 2-D, non-hydrostatic, serial;
  - `b_twoplume3d`: 3-D with two plumes;
  - `c_fjord3d`: a 3-D fjord with wind and shelf forcing, run in two stages from a pickup.
- **`gendata/make_case.py`** writes the input binaries and renders the namelist templates in `config/<cfg>/input_templates/`. It also writes `case.json` with the *actual* parameters, e.g. a source separation snapped to the grid. Analysis uses these actual values, never the nominal ones.
- **`scripts/run_campaign.py`** builds the model and runs the cases. It refuses to start while another MITgcm run is active. For two-stage (c) cases it copies the stage-1 pickups, including the open-boundary pickup needed for an exact restart.
- **`pfcommon/`** holds code shared across experiments:
  - binary and MDS I/O;
  - `runctl.require_ok()`, which makes the analysis refuse runs that are not OK or did not reach `nTimeSteps`;
  - figure styling.
- **`theory/`** holds the plume models (Jenkins 2011 line plume, Slater et al. 2016 half-cone plume, Morton–Taylor–Turner cone) and the three-equation melt law. The `mitgcm` melt mode reproduces `pkg/icefront/icefront_thermodynamics.F` line by line.
- **`notes/`** (in Spanish) is the decision record: pre-registered criteria, results and every follow-up.

## Pipeline (Experiment 2)

The same chain, `cases.yaml` → `gendata` → MITgcm run → `analysis`, with these differences:

- **Geometry** (`gendata/geometry.py`) comes from the ISOMIP+ file `Ocean1_input_geom_v1.01.nc`, which is not tracked: download it into `exp2_isomip_cavity/data_external/` as described, with its DOI and checksum, in `data_external/SOURCES.md` (the MITgcm-JPL output for the post hoc comparison is there too).
- **Run suffix `_r21`.** Every valid run directory ends in `_r21`; runs made before the correction of the friction-velocity law (ruling R21) are void. `case.json` carries `melt_law`, and `analysis/metrics.load_run` refuses any run without the current one, as well as any run that is not OK.
- **Analysis** (`analysis/`): `calibrate.py` and `pd1_calibration_eval.py` (K1, V1), `fig1_ocean0_validation.py` (Ocean1, V2, V3), `ensemble.py`, `fig2_melt_sensitivity.py` and `fig3_meltwater_export.py` (H1, H2), `h1_posthoc_tf.py` (post hoc R23, R24), `jpl_compare.py` (post hoc), `export_exp3.py` (dataset for Experiment 3, written to `$EXP2_RUNS/exp3_dataset.npz`).
- **README:** rendered by `scripts/render_readme.py` from `README.template.md` and `figures/*.json` (`--check` verifies that it is up to date). `figures/void_pre_r21/` holds void records that nothing reads.

## Locations

Run output, builds and the MITgcm source default to the paths of the original machine. Override them with environment variables:
- `EXP1_RUNS`, `EXP2_RUNS`: run-output roots (defaults `/Volumes/LaCie/portfolio_runs/exp1` and `.../exp2`);
- `EXP1_BUILDS`, `EXP2_BUILDS`: build roots (defaults `/Volumes/ESTEBAN_MIT/portfolio_builds/exp1` and `.../exp2`);
- `MITGCM_ROOT`: MITgcm checkout.

They are defined in `pfcommon/paths.py`. `reproduce.sh` also accepts `RUNS`.

## Conventions

- **MITgcm I/O:** 64-bit binaries with single-CPU I/O. Output is MDS diagnostics only.
- **Compute limits:** up to 5 MPI ranks; thread variables (`OPENBLAS/OMP/VECLIB_NUM_THREADS`) set to 1. Each run is timed briefly before its campaign.
- **Rigor:**
  - Success criteria are fixed before running, and a failed criterion is never redefined.
  - Metrics defined after seeing a result are labelled post hoc.
  - Every model–reference difference is classified as (a) our error, (b) missing physics, (c) a numerical choice or (d) an error in the reference.
  - Parameters are never tuned to match a single case.
- **Citations in code:** each equation cites its reference (author, year, equation) and units. Behaviour taken from MITgcm cites `pkg/…/file.F:line`.
