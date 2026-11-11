# Running the experiments

A portfolio of idealized, reproducible MITgcm polar-ocean experiments. Each experiment pairs a model configuration with theory, validation tests and analysis scripts. Experiment 1 (`exp1_fjord_plumes/`) studies subglacial plumes in a fjord; Experiment 2 (`exp2_isomip_cavity/`) an ISOMIP+ ice-shelf cavity; Experiment 3 (`exp3_melt_emulator/`) emulators of the Experiment 2 basal melt.

## Commands

Experiments 1 and 2, and the analysis of Experiment 3, run in one conda environment (numpy, scipy, matplotlib, pyyaml, pytest); the Experiment 3 protocol needs a second one (below):

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
bash exp3_melt_emulator/reproduce.sh --quick      # hypotheses, figures and README from the JSON
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

## Pipeline (Experiment 3)

The emulators are trained on the dataset that `exp2_isomip_cavity/analysis/export_exp3.py` writes outside the repository (`$EXP2_RUNS/exp3_dataset.npz`; its sha256 is checked before every use, and it is read with `allow_pickle=False`). Experiment 3 reads it at the fixed path `/Volumes/LaCie/portfolio_runs/exp2/exp3_dataset.npz` and keeps its run products (cross-validation units, logs, final-test maps) in `/Volumes/LaCie/portfolio_runs/exp3/`; `EXP2_RUNS` does not move them.

- **Two environments.** The protocol (`scripts/run_protocol.py`, `methods/`) needs PyTorch and scikit-learn, which are in a second conda environment, `anemoi`; the analysis and figure scripts (`analysis/`) run under the main environment (`pyocean`). Nothing is installed into either.
- **Tests under anemoi.** That environment has no pytest. `exp3_melt_emulator/scripts/pytest_anemoi.sh [pytest args]` symlinks the pure-Python pytest stack of the main environment into a temporary directory and runs it under anemoi; with no path among the arguments it tests `exp3_melt_emulator`. Under the main environment the modules that need PyTorch or scikit-learn skip themselves, so the whole suite still runs there.

```bash
export ANEMOI=/path/to/anemoi/bin/python                                      # pytest_anemoi.sh and reproduce.sh read ANEMOI (and PYOCEAN / PY)
bash exp3_melt_emulator/scripts/pytest_anemoi.sh -q                       # Exp. 3 tests under anemoi
$PY -m pytest -q                                                          # whole suite (Exp. 3 torch/sklearn modules skip)
$ANEMOI exp3_melt_emulator/scripts/run_protocol.py --freeze-config       # once, at the pre-registration (refuses to overwrite)
$ANEMOI exp3_melt_emulator/scripts/run_protocol.py --cv [--wall-limit-hours H --amendment "<text>"]   # resumable
$ANEMOI exp3_melt_emulator/scripts/run_protocol.py --final               # the final test: runs once, then refuses
bash exp3_melt_emulator/reproduce.sh [--quick]                            # full mode re-runs --cv (compared with, then restored to, the committed record), never --final
```

- **`run_protocol.py`** has three modes. `--freeze-config` writes `figures/protocol_config.json` (dataset sha256, folds, every method's hyperparameters, thresholds, noise floor) and refuses to overwrite it. `--cv` verifies that frozen configuration, the dataset and the folds, fits every method on every non-final fold in up to 5 single-threaded processes, writes one file per unit to `/Volumes/LaCie/portfolio_runs/exp3/cv_units/` (a restart skips finished units) and assembles `figures/cv_results.json`; its frozen wall-clock stop rule can be replaced only together with the text of a documented amendment. `--final` refuses if `figures/final_test.json` or its start marker exists, or if the cross-validation is incomplete; once an attempt has completed there is no override (`--resume-after-crash` repeats, from scratch, only an attempt recorded as crashed).
- **Analysis** (`analysis/`): `noise_floor.py` (noise floor from the Exp. 2 runs), `hypotheses.py` (H3a–H3c and the ranking), `figures.py` (Figures 1, 2, 4), `physics_diag.py` (Figure 3), `profile_test_diag.py` (post hoc, under anemoi), `bootstrap.py` (bootstrap intervals over runs), `compare_records.py` (re-run against committed record, used by `reproduce.sh`).
- **README:** rendered by `scripts/render_readme.py` from `README.template.md` and `figures/*.json` (`--check` verifies that it is up to date); a test fails if the template contains a number that does not come from a placeholder.

## Locations

Run output, builds and the MITgcm source default to the paths of the original machine. Override them with environment variables:
- `EXP1_RUNS`, `EXP2_RUNS`: run-output roots (defaults `/Volumes/LaCie/portfolio_runs/exp1` and `.../exp2`);
- `EXP1_BUILDS`, `EXP2_BUILDS`: build roots (defaults `/Volumes/ESTEBAN_MIT/portfolio_builds/exp1` and `.../exp2`);
- `MITGCM_ROOT`: MITgcm checkout.

They are defined in `pfcommon/paths.py`. `reproduce.sh` also accepts `RUNS`. For Experiment 2, the analysis scripts and `exp2_isomip_cavity/reproduce.sh --quick` honour `EXP2_RUNS` / `RUNS`, but the queue scripts (`exp2_isomip_cavity/scripts/queue_*.sh`, and therefore the full mode of `reproduce.sh`) hardcode `/Volumes/LaCie/portfolio_runs/exp2` and refuse to start unless that disk is mounted.

## Conventions

- **MITgcm I/O:** 64-bit binaries with single-CPU I/O. Output is MDS diagnostics only.
- **Compute limits:** up to 5 MPI ranks; thread variables (`OPENBLAS/OMP/VECLIB_NUM_THREADS`) set to 1. Each run is timed briefly before its campaign.
- **Rigor:**
  - Success criteria are fixed before running, and a failed criterion is never redefined.
  - Metrics defined after seeing a result are labelled post hoc.
  - Every model–reference difference is classified as (a) our error, (b) missing physics, (c) a numerical choice or (d) an error in the reference.
  - Parameters are never tuned to match a single case.
- **Citations in code:** each equation cites its reference (author, year, equation) and units. Behaviour taken from MITgcm cites `pkg/…/file.F:line`.
