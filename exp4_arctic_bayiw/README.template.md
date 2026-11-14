<!-- Rendered by scripts/render_readme.py: edit README.template.md, never README.md. Numbers come from figures/*.json and cases.yaml. -->
# Experiment 4 — Wind-driven internal waves in an Arctic fjord with and without landfast ice: Van Mijenfjorden (Svalbard), and a bridge from an idealised bay

**Esteban Cruz** · MITgcm, realistic fjord geometry and an idealised bay · PhD application portfolio (2027)

## Abstract

A short pulse of wind over a stratified fjord sets the water column oscillating: near-inertial currents, internal Kelvin waves that run around the basin along its coasts, and seiches. In my MSc thesis, a rectangular bay at mid latitude radiated most of that energy out through its mouth, at near-inertial frequencies. Here we ask what happens in a real Arctic fjord, Van Mijenfjorden in Svalbard ({{lat}}°N), whose inner basin is closed by an island and connected to the shelf only through two straits narrower than the internal Rossby radius, and what landfast ice does to the response. We ran MITgcm on a {{box_km}} domain at {{dx_m}} m with the measured bathymetry and representative summer and late-winter stratification, forced by one wind pulse along the fjord axis, and reran the thesis bay (copied bit for bit) as a bridge, at its own latitude and at {{lat}}°N. Criteria and hypotheses were fixed before any run (`notes/pd0_preregistro.md`); the flux that decides the hypotheses was amended after the first runs, before any hypothesis ratio was computed (§2.6). The copy of the thesis reproduces it exactly (V1, maximum surface-elevation difference {{V1_deta}} m), and the model reproduces the tidally forced internal Kelvin wave of Støylen & Fer (2014) in this fjord (V3: phase speed {{V3_c}} m s⁻¹, trapping scale {{V3_decay}} km). {{n_failed_Word}} validation criteria failed and are reported as failed: the near-inertial peak on the shelf is blue-shifted to {{V2a_peak}} f (V2a), the pre-registered Kelvin-front metric does not measure a propagating wave under a spatially uniform wind (V2b), and the net baroclinic energy flux through the straits changes sign between {{dx400_m}} m and {{dx_m}} m (V4). Of the {{n_hyp}} hypotheses, the band partition of the radiated flux (H4a) is inconclusive, because the net flux through the straits is dominated by periods longer than the record resolves ({{H4a_low_share}} of it); landfast ice reduces the internal-wave energy of the inner fjord below half that of open water (H4b(i), ratio {{H4b_ratio}}) but the energy that remains does not enter mainly as a baroclinic flux (H4b(ii)); summer and winter open water are not distinguishable (H4c); and replacing the thesis mixing by KPP raises the radiated fraction of the wind work in the bay (H4d, ratio {{H4d_R_ratio}}; with the flux of Amendment 1 it would be {{a1_R_ratio}}, post hoc, §2.6) but gives a negative fraction in the fjord ({{H4d_F1_eta}} %), so H4d is refuted there by the frozen arithmetic, a sign effect of the low-frequency exchange and not KPP radiating less. The bay radiates outward and mostly near-inertially at both latitudes; the fjord's net is low-frequency and inward. That contrast is descriptive. It is consistent with the geometry rather than f, but it cannot isolate the geometry: the bay runs also differ from the fjord runs in the spatial pattern of the wind (a local patch over the open bay against a uniform pulse along the narrow fjord, which projects differently onto near-inertial motion), the horizontal grid (about {{bay_dx}} km in the bay against {{dx_m}} m), the time step ({{bay_dt}} s against {{deltaT}} s) and the outer domain. Part of the fjord's low-frequency term is an artefact of the fixed density reference of the energy budget (a fraction {{mr_drift_frac}} of it in a post hoc check), and none of these fluxes is resolution-converged. Every number below is read from the metrics files in `figures/` and the case registry `cases.yaml`.

![Figure 0](figures/fig0_schematic.jpg)

*Figure 0.* The setting and the processes; no model result is drawn. (a) Van Mijenfjorden in a Sentinel-2 image of {{fig0_date}}, with landfast ice covering the inner fjord and open water in Bellsund, and the island Akseløya across the mouth; the inset locates it in Svalbard (MODIS). The yellow dashed frame is the model domain ({{box_km}}, open boundaries to the west, north and south, with a sponge of {{n_sponge}} cells, lightly shaded); the red lines are the strait sections of Akselsundet ({{aks_w}} km) and Mariasundet (two passages, the wider {{mar_w}} km), which close the inner-fjord box of the energy budget; the late-winter runs put landfast ice on the sea cells inside the dashed blue polygon; the dots are the analysis stations (yellow: the shelf station of V2a; white: the south-coast stations of V2b). (b) A plan view, not to scale, of the processes the experiment asks about: (1) a pulse of wind along the fjord axis toward the mouth; (2) the near-inertial response, currents turning clockwise at the inertial period; (3) the internal Kelvin wave, trapped within about one Rossby radius of the coast and travelling cyclonically, with the coast on its right; (4) radiation of internal-wave energy through the straits to the shelf; (5) landfast ice, which keeps the wind from the water and puts an ice–ocean drag on its surface. Imagery: Contains modified Copernicus Sentinel data ({{fig0_year}}), processed by ESA; inset: MODIS Terra true colour, NASA Worldview/GIBS. Drawn with VectorCraft by `figures/src/fig0_schematic.py`.

## 1. Introduction

Wind blowing over a stratified fjord moves the surface layer, tilts the pycnocline against the coasts, and leaves behind, when it stops, a field of internal waves: near-inertial oscillations, internal Kelvin waves trapped along the coasts and travelling with the coast on their right in the northern hemisphere, and Poincaré waves and seiches (Gill, 1982). How much of that energy leaves the fjord, and at which frequencies, sets how much is left to mix the fjord and how much feeds the shelf. In a Swedish sill fjord, Arneborg & Liljebladh (2001) found the net baroclinic energy flux at the mouth to be inward; in Svalbard, Støylen & Fer (2014) modelled tidally forced internal Kelvin waves in Van Mijenfjorden.

My MSc thesis (BayIW_Rectan) studied a rectangular bay ({{bay_L}} × {{bay_W}} km, {{bay_H}} m) forced by a wind patch at f = {{f_thesis}} s⁻¹: most of the baroclinic energy flux left through the mouth, at near-inertial frequencies, but only a small fraction of the wind work left the bay, and its discussion proposed that a boundary-layer closure might change that fraction. Experiment 4 takes that question to a real Arctic fjord. At {{lat}}°N the inertial period is {{Ti78}} h instead of {{Ti_thesis}} h, the internal Rossby radius is only {{R1_summer}} km in summer and {{R1_winter}} km in late winter, and the inner basin of Van Mijenfjorden ({{fj_area}} km², {{fj_L}} km long, {{fj_W}} km wide, at most {{fj_Hmax}} m deep) opens to Bellsund only through Akselsundet and Mariasundet, both narrower than the Rossby radius. In winter it is covered by landfast ice.

{{n_hyp_Word}} hypotheses were fixed before any run (`notes/pd0_preregistro.md` §5; the flux that decides them was amended after the first runs, §2.6): **H4a**, the radiated baroclinic flux is mostly near-inertial, as in the thesis bay, or, the alternative, carried by the internal Kelvin wave; **H4b**, landfast ice cuts the internal-wave energy of the inner fjord by at least half, and what remains enters through the straits; **H4c**, summer stratification radiates more than late winter in open water; **H4d**, a boundary-layer closure (KPP; Large et al., 1994) radiates a larger fraction of the wind work than the constant mixing of the thesis. A bridge of runs connects the thesis bay to the fjord (§2.3).

## 2. Methods

### 2.1 Model, domain and bathymetry

MITgcm (checkpoint69j), hydrostatic, on a Cartesian grid in UTM 33N: {{box_km}} at {{dx_m}} m, {{nz}} levels from {{dz_top}} m at the surface to the bottom of the grid ({{H_grid}} m), partial cells. The bathymetry is IBCAO v5.2 (Jakobsson et al., 2024) on the coastline of the Norwegian Polar Institute (S100), with the two straits checked against Støylen & Fer (2014) (`notes/pd1_batimetria.md`, `figures/bathy_report.json`). Akselsundet is {{aks_w}} km wide with a sill at {{aks_sill}} m on the grid; Mariasundet has two passages, the wider {{mar_w}} km, which the grid opens as {{mar_grid_w}} km of water, with a sill at {{mar_sill}} m. In the inner fjord, {{tid_measured}} % of the model cells rest on measured soundings and {{tid_interp}} % on interpolated ones. The western, northern and southern boundaries are open (Orlanski radiation) with a sponge of {{n_sponge}} cells.

### 2.2 Stratification and equation of state

The representative profiles are medians of the UNIS Hydrographic Database casts inside the inner fjord ({{n_summer}} summer casts, August–September; {{n_winter}} late-winter casts, March–April, from only {{n_winter_years}} years), on the model levels and made stable with the least change. The equation of state is linear, fitted to TEOS-10 over the range of the profiles (tα = {{eos_ta}} K⁻¹, sβ = {{eos_sb}}; maximum error {{eos_err}} kg m⁻³). The first-mode phase speed c₁ (flat bottom, rigid lid; `theory/modes.py`) is {{c1_summer}} m s⁻¹ in summer and {{c1_winter}} m s⁻¹ in late winter (mean to maximum depth), so R₁ = c₁/f is {{R1_summer}} km and {{R1_winter}} km and the grid resolves R₁ with {{R1dx_summer}} cells in summer but only {{R1dx_winter}} in winter. A Kelvin wave circles the inner fjord ({{fj_P}} km) in {{TK_summer}} h in summer. Because the fjord is narrow, the first cross-fjord Poincaré mode of the first vertical mode lies at {{poincare1}} f (Gill, 1982, §10.3, Eq. 10.3.5), well above the near-inertial band; in the thesis bay it lies at {{poincare_bay}} f, inside it.

### 2.3 Forcing and runs

The fjord runs are forced by one pulse of wind stress, a half Gaussian in time of peak {{tau_max}} N m⁻², width σ = {{sigma_h}} h and centre at {{center_h}} h, spatially uniform, blowing along the fjord axis toward the mouth (from {{wind_from}}° grid; down-fjord winds are the strong, channelled winds of Svalbard fjords, Frank et al., 2023), except {{cross_axis_runs}}, where the same pulse blows across the axis. They are forced for {{forced_h}} h and then free, {{t_end_h}} h in all, with Δt = {{deltaT}} s (staggered); V3 has no wind and the M2 tide at its open boundaries instead. The bay runs keep the thesis set-up: the same pulse in time (peak {{bay_tau}} N m⁻², σ = {{bay_sigma}} h, centre at {{bay_center}} h), applied over a circular patch of radius {{bay_patch_r}} grid cells centred on the bay; the thesis grid (about {{bay_dx}} km in the bay, the open sea limited to the depth of the thesis grid) and Δt = {{bay_dt}} s; forced for {{bay_forced}} h, and {{R1_t_end}} h in all for R1, {{R1k_t_end}} h for R1k and {{R2_t_end}} h for R2 and R3.

| Run | Domain | Stratification | Surface | Vertical mixing | Role |
|---|---|---|---|---|---|
| F1 | fjord, {{dx_m}} m | summer | open water | KPP | reference |
| F2 | fjord | late winter | landfast ice | KPP | H4b |
| F3 | fjord | late winter | open water | KPP | H4b, H4c |
| F4 | fjord | late winter | drifting ice | KPP | descriptive |
| F5 | fjord | summer | open water, cross-axis wind | KPP | descriptive |
| F1c | fjord | summer | open water | constant, as the thesis | H4d |
| F1g | fjord | summer | open water | GGL90 | descriptive |
| V3 | fjord | summer | no wind, M2 tide at the open boundaries | KPP | V3 |
| V4 | fjord, {{dx400_m}} m | summer | open water | KPP | V4 |
| R1 | thesis bay, thesis f | thesis | open water, thesis wind patch | thesis | V1, H4d |
| R1k | thesis bay, thesis f | thesis | open water, thesis wind patch | KPP | H4d |
| R2 | thesis bay at {{lat}}°N | fjord summer | open water, thesis wind patch | KPP, fjord numerics | bridge |
| R3 | as R2 without the bay | fjord summer | open water | KPP, fjord numerics | bridge |

R1 is the thesis run reproduced with the thesis code and namelists; only the output format changed. R1k adds KPP and nothing else that a namelist comparison does not list (`notes/pd1_puente.md`); R1k runs {{R1k_t_end}} h and is compared with R1 over the {{R1_window}} h of the thesis. R2 changes f, profile, equation of state and the fjord mixing at once; it is the thesis bay in the fjord's ocean.

### 2.4 Landfast ice

Landfast ice is MITgcm sea ice ({{ice_h}} m thick, concentration {{ice_A}}) over the inner fjord with the ice dynamics off, so the ice does not move, without thermodynamics, and the ocean under it receives the ice–ocean drag τ = ρ₀ C_w |u_ice − u| (u_ice − u) with the ice at rest and the default C_w of MITgcm, in the range measured under Arctic ice (Fer et al., 2022). MITgcm does not compute that drag when no rheology solver runs, so MITgcm files are replaced by modified copies (`config/<cfg>/code/` of every fjord configuration, {{mod_configs}}, each with the sha256 of the original and every change between marked lines): `seaice_dynsolver.F` computes the drag coefficient also without a rheology solver, so that the ocean under fixed ice feels it; `obcs_check.F` skips two stops that forbid Orlanski open boundaries and sea-ice open-boundary conditions together, which is safe here because the ice stays at least {{ice_min_cells}} cells from the open boundaries in every record of F2 and F4 (`figures/ice_ob_check.json`); `obcs_add_tides.F` ramps the tide of V3 up over one M2 cycle. The F2 and F4 restarts are bit-identical (`figures/pd1_restart.json`).

### 2.5 Mixing and numerics

KPP with MITgcm defaults and the background viscosity and diffusivity of Large et al. (1994); Smagorinsky horizontal viscosity with a small biharmonic part and quadratic bottom drag, the configuration Støylen & Fer (2014) used for this fjord; a second-order flux-limited advection of temperature and salinity, whose numerical mixing is declared (Klingbeil et al., 2014). The time step is staggered (`notes/pd1_config.md`).

### 2.6 Energy fluxes, budget and bands

The fluxes and the energy budget follow the hydrostatic decomposition of Kang & Fringer (2012, Eqs. 13–17), as the thesis does: barotropic and baroclinic kinetic energy, the exact available potential energy relative to a fixed background profile ρ_b (the initial one), the conversion C, the wind work P, and on each open face the pressure work F and the energy advection A; the residual is D. The box is the inner fjord closed by stepped sections across Akselsundet and Mariasundet, or the bay closed at its mouth. On R1, the implementation reproduces the thesis scripts with a maximum relative difference of {{kf_rel}} for the budget and {{nash_rel}} for the Nash et al. (2005) flux of the thesis (tolerance {{contrast_tol}} %). The flux that decides the hypotheses is the baroclinic pressure work F_bc (Amendments 1–3 of the pre-registration; advection is reported as descriptive). No hypothesis ratio had been computed when they were written, but two flux results were already known, as the pre-registration records. Amendment 1 had put the advection into the H4d numerator (F_bc + A_bc); Amendment 2 took it out again knowing that on R1 the pressure work carries {{r1_F_out}} GJ out of the bay and the advection {{r1_A_in}} GJ into it. Amendment 3 made the pressure work decide H4a, H4b(ii) and H4c(ii) when the band split of the F1 pressure work through each strait (the H4a input) was already in the V4 diagnostics (`figures/pd2_v4_diagnostics.json`). With the Amendment 1 numerator, the bay ratio of H4d would be {{a1_R_ratio}} (η = {{a1_R1k_eta}} % against {{a1_R1_eta}} %), which the frozen rule would read as "{{a1_R_outcome}}": the one supported verdict of this experiment depends on that amendment (post hoc, descriptive). A flux is split into frequency bands by the exact cospectral partition of its time mean, with bands fixed before the runs: low (periods longer than {{low_period}} h, not resolved by the record), Kelvin/subinertial ({{k_period}} h), near-inertial ({{ni_lo}}–{{ni_hi}} f, {{ni_period}} h, {{ni_bins}} frequency bins), superinertial ({{si_period}} h) and seiche (shorter than {{seiche_period}} h).

### 2.7 Pre-registered criteria

Every model–reference difference is classified as (a) our error, (b) missing physics, (c) a numerical choice or (d) an error in the reference, and no failed criterion is redefined.

- **V1:** R1 equals the saved thesis run (surface elevation within {{V1_tol}} m).
- **V2a:** in F1, the near-inertial peak of surface kinetic energy at a shelf station lies within one frequency bin of f. **V2b:** the internal Kelvin wave speed, from the arrival of the pycnocline-displacement front at successive coastal stations, lies in {{v2b_range}} m s⁻¹.
- **V3:** with the M2 tide alone, an internal Kelvin wave trapped at the coast (decay scale {{V3_decay_range}} km), travelling cyclonically, at {{V3_c_range}} m s⁻¹ (within a factor of two of Støylen & Fer, 2014).
- **V4:** V4 ({{dx400_m}} m) against F1: the net outward baroclinic flux within {{V4a_tol}} % (V4a) and the share of each resolved band within {{V4b_tol}} (V4b).
- **Hypotheses:** a ratio must exceed the resolution uncertainty, so the minimum factor is {{factor}} (a ratio between {{factor_inv}} and {{factor}} is "not distinguishable"). H4a: if the low band carries more than {{H4a_low_max}} of the net flux, inconclusive; otherwise supported if the near-inertial share exceeds {{H4a_phi}}, the alternative if the Kelvin share does. H4b: (i) baroclinic energy of the inner fjord without the low band, F2/F3 ≤ {{H4b_max}}; (ii) the baroclinic inflow through the straits at least as large as the barotropic one. H4c: (i) near-inertial wind work and (ii) radiated flux, F1/F3 ≥ {{factor}}. H4d: η = F_bc,out / P, the fraction of the wind work radiated, with KPP over the thesis closure ≥ {{factor}}, in the bay (R1k/R1) and in the fjord (F1/F1c).

## 3. Results

### 3.1 The bridge to the thesis (V1) and the method contrast

R1 is bit-identical to the saved thesis run: the surface elevation differs by {{V1_deta}} m in all {{V1_n}} common records, and the restart files at the end of both stages are identical byte for byte (V1 **{{V1_verdict}}**). The thesis therefore enters this experiment as it was, including two features found on the way and classified (d), errors of the reference: the thesis ran on a β-plane, not an f-plane, and its spectral scripts placed f at a value that the β-plane reaches only near the southern boundary of the domain, far from the bay (`notes/pd1_puente.md` §7).

### 3.2 Validation: {{n_failed}} criteria fail

![Figure 1](figures/fig1_validation.png)

*Figure 1.* Validation against the frozen criteria. (a) V2a: surface kinetic energy spectrum at the shelf station of F1, the near-inertial band shaded, f dashed, the peak marked. (b) V2b: arrival of the half-maximum pycnocline displacement along the south coast, with the pre-registered speed range. (c, d) V3: decay of the M2 pycnocline amplitude away from the coast and its phase along the south coast. (e) V4: share of the outward flux in each resolved band at {{dx_m}} m and {{dx400_m}} m.

| Criterion | Threshold | Result | Verdict | Class |
|---|---|---|---|---|
| V2a: near-inertial peak on the shelf | within one bin of f | {{V2a_peak}} f, {{V2a_off}} bins | **{{V2a_verdict}}** | (c) |
| V2b: Kelvin front speed | {{v2b_range}} m s⁻¹ | {{V2b_c}} m s⁻¹ ({{V2b_n}} stations) | **{{V2b_verdict}}** | (a) |
| V3: trapping, sense, speed | {{V3_decay_range}} km; cyclonic; {{V3_c_range}} m s⁻¹ | {{V3_decay}} km; cyclonic; {{V3_c}} m s⁻¹ | **{{V3_verdict}}** | — |
| V4a: net outward flux, {{dx400_m}} m against {{dx_m}} m | {{V4a_tol}} % | {{V4_F1}} GJ against {{V4_V4}} GJ ({{V4a_rel}} %) | **{{V4_verdict}}** | (c) |
| V4b: band shares | {{V4b_tol}} | differences {{V4b_min}}–{{V4b_max}} | **{{V4_verdict}}** | (c) |

**V2a** fails: the near-inertial peak on the shelf is at {{V2a_peak}} f, {{V2a_off}} bins from f. The blue shift is domain-wide (post hoc: median {{v2a_med}} f over the wet columns, central range {{v2a_p5}}–{{v2a_p95}} f; only {{v2a_frac}} % of the columns peak within a bin of f) and the motion is clockwise. A post hoc column test, fixed before it was run, put the F1 configuration in a doubly periodic slab: the peak came at {{col_peak}} f ({{col_off}} bins), which excludes an error in the Coriolis term or the time stepping and, by the rule fixed beforehand, makes the class (c). The test excludes only that: the cause of the blue shift is not identified. The candidates are the bounded geometry and its lateral modes, the open boundaries and the sponge, and background vorticity or mean flow; the {{box_km}} domain width alone gives too small a shift (Gill, 1982, §10.3).

**V2b** fails, class (a): our metric. With a spatially uniform wind every coastal station responds locally during the pulse ({{V2b_early}} of {{V2b_n}} stations reach the half-maximum before the pulse is over, {{pulse_end_h}} h), so the fitted {{V2b_c}} m s⁻¹ is not a phase speed (Gill, 1982, §10.10, Eq. 10.10.21; Yao, 1986). A post hoc, descriptive free-wave speed (crests of the band-passed displacement after the pulse; Hodges et al., 2000; Beletsky et al., 1997) is {{fw_c}} ± {{fw_se}} m s⁻¹, which does not separate a free wave from the forced response. V3 measures the model's internal Kelvin wave directly.

**V3** passes: the tidally forced internal Kelvin wave is trapped within {{V3_decay}} km of the coast, travels eastward along the south coast, and its phase speed is {{V3_c}} m s⁻¹, as in Støylen & Fer (2014).

**V4** fails, class (c): at {{dx400_m}} m (R₁/Δx = {{R1dx_summer_400}}) the net baroclinic flux through the straits is outward ({{V4_V4}} GJ) while at {{dx_m}} m it is inward ({{V4_F1}} GJ). The gross flux falls to {{v4_aks_ratio}} of its {{dx_m}} m value at Akselsundet, which carries {{aks_gross_share}} % of it, and to {{v4_mar_ratio}} at Mariasundet, which the coarse grid barely opens. By the pre-registered rule the limitation is declared. A consequence follows: the resolution bound of {{V4a_tol}} % that justified the factor {{factor}} of the hypotheses is not demonstrated.

### 3.3 Hypotheses

| Hypothesis | Result | Outcome |
|---|---|---|
| H4a: partition of the radiated flux (F1) | low band {{H4a_low_share}} of the net | **{{H4a_outcome}}** |
| H4b(i): energy under landfast ice, F2/F3 | {{H4b_EF2}} / {{H4b_EF3}} GJ = {{H4b_ratio}} | **{{H4b_i}}** |
| H4b(ii): inflow through the straits (F2) | baroclinic {{H4b_Fbc}} GJ against barotropic {{H4b_Fbt}} GJ | **{{H4b_ii}}** |
| H4c: summer against winter (F1/F3) | (i) {{H4c_i}}; (ii) {{H4c_ii}}, both fluxes inward | **{{H4c_outcome}}** |
| H4d: KPP against the thesis closure, bay (R1k/R1) | η = {{H4d_R1k_eta}} % against {{H4d_R1_eta}} %, ratio {{H4d_R_ratio}} (F + A, descriptive: {{a1_R_ratio}}; §2.6) | **{{H4d_R_outcome}}** |
| H4d: KPP against the thesis closure, fjord (F1/F1c) | η = {{H4d_F1_eta}} % against {{H4d_F1c_eta}} %, ratio {{H4d_F_ratio}} | **{{H4d_F_outcome}}** (by sign) |

![Figure 2](figures/fig2_spectra_partition.png)

*Figure 2.* The fjord against the thesis bay and the bridge. (a) Surface kinetic energy spectra: F1 (south-coast mean and shelf station), R1k (thesis bay), R2 (the bay at {{lat}}°N) and R3 (no bay), with f of each latitude. (b) Share of the surface kinetic energy in each band at those stations. (c) Signed share F_b / Σ₅|F_b| of the net outward baroclinic pressure work in each band, through the straits of the fjord runs and the mouth of the bay; positive is outward.

**H4a is inconclusive.** The net baroclinic pressure work through the straits of F1 over the record is {{H4a_low}} GJ in the low band (inward), against {{H4a_K}} (Kelvin), {{H4a_NI}} (near-inertial), {{H4a_SI}} (superinertial) and {{H4a_S}} GJ (seiche): the low band carries {{H4a_low_share}} of the total, above the {{H4a_low_max}} limit, so by the rule evaluated first H4a is "{{H4a_outcome}}". The resolved bands sum to {{H4a_res}} GJ, an inflow, so their shares (Kelvin {{H4a_phi_K}}, near-inertial {{H4a_phi_NI}}) fall outside the unit interval and mean nothing. With advection included the low share is {{H4a_FA_share}} (descriptive).

**H4b: landfast ice removes most of the internal-wave energy, but what remains does not enter as a baroclinic flux.** (i) Without the low band, the time-mean baroclinic energy of the inner fjord is {{H4b_EF2}} GJ under landfast ice and {{H4b_EF3}} GJ in open water, a ratio of {{H4b_ratio}}: met. Most of the F2 value is available potential energy from the drift of the background profile (§3.4); the post hoc ratio of the baroclinic kinetic energy alone is {{H4b_ke_ratio}}, so the outcome holds with either number. (ii) The baroclinic pressure work entering through the straits ({{H4b_Fbc}} GJ) is smaller than the barotropic one ({{H4b_Fbt}} GJ): not met.

![Figure 4](figures/fig4_ice_maps.png)

*Figure 4.* Depth-integrated, time-mean baroclinic kinetic energy in each resolved band in the inner fjord, under landfast ice (F2, top) and in open water with the same late-winter profile (F3, bottom); one logarithmic colour scale per band.

**H4c is not distinguishable.** The near-inertial wind work in summer is {{H4c_WF1}} GJ and in late winter {{H4c_WF3}} GJ, a ratio of {{H4c_i}}, in neither tail. The radiated flux ratio, {{H4c_ii}}, is the quotient of two inward fluxes ({{H4c_F1}} and {{H4c_F3}} GJ) and is flagged; it cannot change the outcome. The winter first mode is poorly resolved at {{dx_m}} m (R₁/Δx = {{R1dx_winter}}).

**H4d is supported in the bay and refuted, by sign, in the fjord.** In the bay over the first {{R1_window}} h, KPP radiates η = {{H4d_R1k_eta}} % of the wind work against {{H4d_R1_eta}} % with the thesis closure (the thesis value), a ratio of {{H4d_R_ratio}}: supported. Most of the change is in the denominator: with KPP the wind does {{H4d_W_ratio}} times less work ({{H4d_WR1k}} against {{H4d_WR1}} GJ), while the radiated flux grows only by a factor of {{H4d_num_ratio}}. In the fjord, the net baroclinic pressure work through the straits in F1 is inward ({{H4d_F1_num}} GJ), so η(F1) = {{H4d_F1_eta}} % against {{H4d_F1c_eta}} % in F1c and the frozen rule gives "refuted" with a flag on the negative numerator. The refutation comes from the sign of the low-frequency net exchange, not from KPP radiating less. Post hoc, with the {{n_resolved}} resolved bands only, η is {{ph_eta_F1}} % (F1) and {{ph_eta_F1c}} % (F1c); with the closure of the thesis the Kelvin band carries the net outward ({{F1c_K}} GJ).

![Figure 3](figures/fig3_strait_flux_budget.png)

*Figure 3.* (a) Baroclinic pressure work into the inner fjord through each strait in F1. (b) Through both straits in F1, F1c, F3 and F2. (c, d) Cumulative baroclinic budget of the inner fjord in F1 and F2: wind work, conversion, inflow (pressure work and advection), change of energy and the residual D; the growth of the energy in F2, without wind over the fjord, is the drift of the fixed reference (§3.4).

### 3.4 Post hoc checks

These were written after the results they examine and change no verdict.

- **Drift of the fixed reference.** The background profile of the inner fjord is smoothed during every long run (largest change {{drift_F1c}} kg m⁻³ in F1c, {{drift_F1}} in F1, {{drift_F2}} in F2, {{drift_bay}} in the thesis bay). Measured against the fixed initial ρ_b that the budget prescribes, that appears as available potential energy that grows even without wind: {{ape_F2_end}} GJ at the end of F2 (its baroclinic kinetic energy ends at {{ke_F2_end}} GJ), {{ape_F1_end}} GJ in F1 and {{ape_R2_end}} GJ in R2. The residual D of these budgets is therefore not dissipation, class (c). The same drift enters the pressure perturbation of the strait fluxes.
- **Moving reference.** With the density reference taken as the box-mean profile at each time (method fixed before computing), the low band of F1 is still inward and dominant ({{mr_F1_low}} GJ instead of {{H4a_low}} GJ; low share {{mr_F1_share}}), and η(F1) is still negative ({{mr_F1_eta}} %). The drift explains a fraction {{mr_drift_frac}} of the low band, not its sign. But the resolved sum of F1 changes sign ({{mr_F1_res}} GJ, mostly the Kelvin band, {{mr_F1_K}} GJ), and in F1c the low band turns outward ({{mr_F1c_low}} GJ). The band shares of the fjord are not robust to the reference. In the bay the reference hardly matters (the H4d ratio would be {{mr_ratio_R}}).
- **Arneborg & Liljebladh band.** Band-passed as in Gullmar fjord (Arneborg & Liljebladh, 2001), the baroclinic pressure work enters the fjord in F1 ({{al_F1}} GJ) and leaves it in F1c ({{al_F1c}} GJ) and in the bay (R1 {{al_R1}}, R2 {{al_R2}} GJ, negative is outward). A sign comparison with another fjord, without resolution backing.
- **Narrow straits.** Mariasundet is {{mar_wR1}} R₁ wide in summer and Akselsundet {{aks_wR1}} R₁; Kelvin waves are transmitted poorly through straits much narrower than the Rossby radius (Johnson & Garrett, 2006, as cited by Jackson et al., 2018).

### 3.5 The bridge (descriptive)

Of the net outward baroclinic pressure work, the near-inertial band carries {{br_R1k_ni}} and the low band {{br_R1k_low}} in the thesis bay (R1k), and {{br_R2_ni}} and {{br_R2_low}} in the same bay at {{lat}}°N with the fjord's stratification and mixing (R2): the bay radiates outward and mostly near-inertially at both latitudes. In the fjord the signed low share is {{br_F1_low}} (inward) and the near-inertial share {{br_F1_ni}}. The surface kinetic energy points the same way: {{br_R2_station_ni}} near-inertial inside R2, but on the south coast of F1 {{br_F1_coast_k}} in the Kelvin band and {{br_F1_coast_ni}} near-inertial. Since R2 keeps the bay's near-inertial character at the fjord's f, f alone does not explain the change. It is consistent with the geometry (narrow basin, narrow straits), as the pre-registration predicted from the Poincaré frequencies (§2.2), but R2 and F1 also differ in the spatial pattern of the wind (a local patch over the open bay against a uniform pulse along the narrow fjord, which projects differently onto near-inertial motion), the horizontal grid (about {{bay_dx}} km in the bay against {{dx_m}} m), the time step ({{bay_dt}} s against {{deltaT}} s) and the outer domain (§2.3), so the bridge cannot attribute the change to the geometry alone. The fjord's low-frequency term carries the drift caveat of §3.4, and none of these fluxes is resolution-converged (V4).

### 3.6 Animations (descriptive)

Four animations show the response to the wind pulse in plan view over the first four days. They are descriptive: no number in this README comes from them, and no criterion or hypothesis uses them (`analysis/animate.py`). The field a₁ is the first-mode isopycnal displacement: the isopycnal displacement relative to the fixed initial reference profile of the energy budget, projected with the N² weight on the local first vertical mode (scaled to one at its maximum), positive upward. The filtered field a₁′ is a₁ minus its centred running mean over exactly one inertial period of each run ({{Ti78}} h in the fjord, {{Ti_thesis}} h in the thesis bay); frames within half a period of the ends of the record have no full window and are left out. Each animation has its own fixed colour scale.

![Animation 1](figures/anim_fjord_ice_raw.gif)

*Animation 1.* a₁ in late winter, open water (F3, left) and landfast ice (F2, right), same profile and wind pulse, one shared scale; the purple line is the initial ice edge and the black lines are the strait sections. Features that keep the same sign from frame to frame (coastal fringes, a strip on the shelf, the open-boundary sponges) may be wind mixing or the drift of the reference profile rather than internal waves, and under the ice of F2 there is no wind stress, so the two panels differ in mixing as well as in waves.

![Animation 2](figures/anim_fjord_ice_filt.gif)

*Animation 2.* The same runs, a₁′. In open water the inner fjord shows bands of opposite sign along its two coasts that change from frame to frame; under the ice a₁′ of the inner fjord stays close to zero on this scale (no wind stress under the ice in F2), while outside the ice the two runs look alike. The running-mean filter attenuates sub-inertial signal (Kelvin waves, slow seiches) and has ripple above f, so a₁′ shows mainly near- and super-inertial oscillations.

![Animation 3](figures/anim_bay_vs_fjord_raw.gif)

*Animation 3.* a₁ in the thesis bay (R1, left, cropped to the bay and the shelf in front of it, scale set by the bay) and in Van Mijenfjorden in summer (F1, right). The scales are separate: the two runs differ in profile, depth, f and the spatial pattern of the wind, so the amplitudes are not comparable. The caveat of Animation 1 on same-sign features applies, and the bridge is not resolution-converged (§3.5).

![Animation 4](figures/anim_bay_vs_fjord_filt.gif)

*Animation 4.* The same runs, a₁′, each filtered with its own inertial period. In the bay the two side walls carry displacements of opposite sign; in the fjord the signal is spread over the basin and the shelf. As in Animation 2, the filter attenuates sub-inertial signal (Kelvin waves, slow seiches) and has ripple above f: a₁′ shows mainly near- and super-inertial oscillations.

Two further animations draw one isopycnal in three dimensions over the model bathymetry, with the land from IBCAO, seen from the southwest over Bellsund looking up the fjord (`analysis/blender_export.py`, `analysis/blender_render.py`; Blender). They are descriptive in the same sense. The surface is the depth of the isopycnal in each density snapshot; its colour is the displacement of the isopycnal from its initial depth, positive upward, minus its centred running mean over one inertial period, with the frames that lack a full window left out. The running-mean filter attenuates sub-inertial signal (Kelvin waves, slow seiches), so the colours show mainly near- and super-inertial oscillations, not the Kelvin waves of the fjord. The vertical exaggeration is {{bl_ve_sea}}× below sea level and {{bl_ve_land}}× on land, and each animation has one fixed colour scale.

[![Animation 5](figures/blender_shotA_F1.gif)](figures/blender_shotA_F1.mp4)

*Animation 5.* Summer, open water (F1): the isopycnal {{blA_sigma}} kg m⁻³, at {{blA_z}} m at the start, the level of strongest stratification below the surface mixed layer; scale ±{{blA_lim}} m. Click for the MP4.

[![Animation 6](figures/blender_shotB_F3F2.gif)](figures/blender_shotB_F3F2.mp4)

*Animation 6.* Late winter, open water (F3, left) and landfast ice (F2, right; the white slab marks the ice cells), same profile and wind pulse, one shared scale of ±{{blB_lim}} m. The winter profile is a staircase whose surface step the pulse mixes away within hours, so the isopycnal is {{blB_sigma}} kg m⁻³, at {{blB_z}} m: of the steps above the sill of Akselsundet, the one whose isopycnal stays in the largest share of the water columns over the whole record in both runs; it sits on a density step of only {{blB_step}} kg m⁻³, so metre-scale values along the coasts may be mixing rather than internal waves. Under the ice of F2 there is no wind stress, so the two panels differ in forcing and mixing as well as in waves. The open boundaries and their {{blB_sponge}}-cell sponges are masked in grey. Click for the MP4.

## 4. Discussion

### 4.1 What the experiment shows

In a real Arctic fjord closed by an island, a wind pulse does not radiate a clean near-inertial flux through the mouth, as the idealised bay did. The baroclinic energy flux through the narrow straits is a small net residual of a larger gross exchange, dominated over a {{t_end_h}} h record by periods too long to resolve, and inward in the reference run. Under landfast ice, with no wind work over the fjord, the internal-wave kinetic energy of the inner fjord nearly vanishes. In the bay, KPP raised the radiated fraction of the wind work mainly by lowering the wind work itself.

### 4.2 What we cannot claim

- That the fjord's net flux is physically inward: part of it is a reference artefact (§3.4), and it changes sign at {{dx400_m}} m.
- Any partition of the radiated flux by frequency in the fjord (H4a inconclusive; the shares change with the reference and the resolution).
- That KPP radiates less in the fjord: the refutation of H4d in the fjord is a sign effect.
- Anything about other fjords, other wind events or tidal conditions: one fjord, one event.

### 4.3 Limitations

- **V2a FAIL, class (c)** by the rule fixed before the column test, which excludes only an error in the inertial frequency; the cause of the blue shift (median {{v2a_med}} f) is not identified (candidates: bounded geometry, open boundaries and sponge, background vorticity).
- **V2b FAIL, class (a):** the front-arrival metric is ill-posed under a spatially uniform wind (Gill, 1982, §10.10, Eq. 10.10.21; Yao, 1986). The post hoc free-wave speed is {{fw_c}} ± {{fw_se}} m s⁻¹.
- **V4 FAIL, class (c): resolution.** The band partition is not robust, and the {{V4a_tol}} % resolution bound behind the factor {{factor}} of the hypotheses is not demonstrated.
- **The low band** (periods longer than {{low_period}} h, unresolved with T_fin = {{t_end_h}} h) dominates the net strait flux; H4a is inconclusive.
- **Fixed-reference drift:** the growth of the exact available potential energy and the residual D are not dissipation. The drift partly contaminates the low-band strait pressure work; the post hoc moving-reference check removes a fraction {{mr_drift_frac}} of it and the sign stays. The F1 band fractions are not robust.
- **H4d in the fjord** is "refuted" literally, from a negative η (net inward low-frequency exchange); η from the resolved bands is {{ph_eta_F1}} % against {{ph_eta_F1c}} % (post hoc).
- **Winter resolution:** the first mode is poorly resolved at {{dx_m}} m (R₁/Δx = {{R1dx_winter}}); this weighs on H4b and H4c.
- **Bathymetry:** {{tid_interp}} % of the inner-fjord cells rest on interpolated depths and {{tid_land}} % on cells that IBCAO codes as land, so only {{tid_measured}} % rest on measured soundings (IBCAO type identifiers, `figures/bathy_report.json`).
- **Modified MITgcm files:** `obcs_check.F`, `seaice_dynsolver.F` and `obcs_add_tides.F` are modified copies (§2.4), declared as our own code.
- **Mariasundet** is {{mar_wR1}} R₁ wide: reduced Kelvin transmission is expected (Johnson & Garrett, 2006, as cited by Jackson et al., 2018), and the grid opens it wider than the passages.
- **One fjord and one wind event; no tides in the study physics** (only V3 has the tide); **ice without thermodynamics.**
- **Late-winter stratification** rests on {{n_winter}} casts from {{n_winter_years}} years.
- **Two analysis details that affect descriptive numbers only:** the bay records hold one sample more than T/Δt, so their low band takes three frequency bins instead of two (the bridge shares and the R1/R1k post hoc numbers); and the near-inertial wind work of H4c(i) is band-passed with an unwindowed circular FFT, so its absolute values carry wrap-around from the end of the record (F1 and F3 are treated alike).

## 5. Conclusions

- The thesis bay was reproduced bit for bit, and the model reproduces the tidally forced internal Kelvin wave of Van Mijenfjorden (V3). {{n_failed_Word}} validation criteria failed and are declared: a blue-shifted near-inertial peak of unidentified cause, an ill-posed Kelvin-front metric, and a strait flux that is not resolution-converged.
- In the fjord, the net baroclinic flux through the straits is dominated by unresolved low frequencies and is inward in the reference run, so the frequency partition of the radiated energy (H4a) is inconclusive.
- Landfast ice nearly removes the internal-wave energy of the inner fjord (H4b(i) met), and what remains does not enter mainly as a baroclinic flux (H4b(ii) not met). Summer and winter open water are not distinguishable (H4c).
- KPP raises the radiated fraction of the wind work in the bay, mostly by reducing the wind work (H4d supported); in the fjord the frozen rule gives "refuted" because the net flux is inward.
- The bay radiates outward and near-inertially at both latitudes; the fjord does not. Descriptively this is consistent with the geometry rather than f, but the bridge runs also differ in the wind pattern, the grid, the time step and the outer domain, and the fjord's flux carries the drift caveat. To measure the fjord's radiation, the next step is higher resolution in the straits, longer records, and an energy budget with an evolving reference profile.

## Reproducing

```bash
bash exp4_arctic_bayiw/reproduce.sh --quick     # every JSON, figure and this README from existing run output
bash exp4_arctic_bayiw/reproduce.sh --dry-run   # the full plan: backup, every campaign, --quick, comparison
bash exp4_arctic_bayiw/reproduce.sh --quick --animations   # also redraws the four animations of §3.6 (about 2 h more)
bash exp4_arctic_bayiw/reproduce.sh             # full: inputs, builds and every run, then --quick and the comparison
```

The quick mode reads the runs (`/Volumes/LaCie/portfolio_runs/exp4`), the saved thesis run and the external data listed in `data_external/SOURCES.md`, and recomputes every cached analysis, so it takes hours. The full mode needs MITgcm, MPI and gfortran and writes to a new run root, so the original runs are never overwritten; it backs up the committed JSONs and compares them with the regenerated ones. Both refuse to start while a model run is active. See [RUNNING.md](../RUNNING.md). Notes in Spanish, the decision record: `notes/pd0_preregistro.md` (pre-registration and its amendments), `pd1_*.md`, `pd2_validacion.md`, `pd3_hipotesis.md`.

## References

- Arneborg, L., & Liljebladh, B. (2001). The internal seiches in Gullmar Fjord. Part I: Dynamics. *J. Phys. Oceanogr.*, 31, 2549.
- Beletsky, D., O'Connor, W. P., Schwab, D. J., & Dietrich, D. E. (1997). Numerical simulation of internal Kelvin waves and coastal upwelling fronts. *J. Phys. Oceanogr.*, 27, 1197.
- Fer, I., et al. (2022). Upper-ocean turbulence structure and ocean-ice drag coefficient estimates during the MOSAiC drift. *J. Geophys. Res. Oceans*. https://doi.org/10.1029/2022JC018751
- Frank, L., Jonassen, M. O., Skogseth, R., & Vihma, T. (2023). Atmospheric climatologies over Isfjorden, Svalbard. *J. Geophys. Res. Atmos.* https://doi.org/10.1029/2022JD038011
- Gill, A. E. (1982). *Atmosphere–Ocean Dynamics*. Academic Press, International Geophysics Series 30 (read in an Academic Press reprint of the 1982 edition). ISBN 978-0-12-283522-3.
- Hodges, B. R., Imberger, J., Saggio, A., & Winters, K. B. (2000). Modeling basin-scale internal waves in a stratified lake. *Limnol. Oceanogr.*, 45, 1603–1620.
- Jackson, R. H., Lentz, S. J., & Straneo, F. (2018). The dynamics of shelf forcing in Greenlandic fjords. *J. Phys. Oceanogr.*, 48, 2799. https://doi.org/10.1175/JPO-D-18-0057.1
- Jakobsson, M., et al. (2024). The International Bathymetric Chart of the Arctic Ocean Version 5.0. *Sci. Data*. https://doi.org/10.1038/s41597-024-04278-w
- Johnson, H. L., & Garrett, C. (2006). What fraction of a Kelvin wave incident on a narrow strait is transmitted? *J. Phys. Oceanogr.*, 36, 945–954. https://doi.org/10.1175/JPO2896.1 (not read: cited as in Jackson et al., 2018)
- Kang, D., & Fringer, O. (2012). Energetics of barotropic and baroclinic tides in the Monterey Bay area. *J. Phys. Oceanogr.*, 42, 272. https://doi.org/10.1175/JPO-D-11-039.1
- Klingbeil, K., Mohammadi-Aragh, M., Gräwe, U., & Burchard, H. (2014). Quantification of spurious dissipation and mixing – discrete variance decay in a finite-volume framework. *Ocean Model.*, 81, 49–64. https://doi.org/10.1016/j.ocemod.2014.06.001
- Large, W. G., McWilliams, J. C., & Doney, S. C. (1994). Oceanic vertical mixing: a review and a model with a nonlocal boundary layer parameterization. *Rev. Geophys.*, 32, 363–403. https://doi.org/10.1029/94RG01872
- Nash, J. D., Alford, M. H., & Kunze, E. (2005). Estimating internal wave energy fluxes in the ocean. *J. Atmos. Oceanic Technol.*, 22, 1551.
- Støylen, E., & Fer, I. (2014). Tidally induced internal motion in an Arctic fjord. *Nonlin. Processes Geophys.*, 21, 87–100. https://doi.org/10.5194/npg-21-87-2014
- Yao, T. (1986). The response of currents in Trinity Bay, Newfoundland, to local wind forcing. *Atmos.–Ocean*, 24, 235–252. https://doi.org/10.1080/07055900.1986.9649249
- Data: IBCAO v5.2 (IBCAO Version 5.2 Compilation Group, 2026), https://doi.org/10.17043/ibcao-5.2; Norwegian Polar Institute, Kartdata Svalbard 1:100 000 (S100), https://doi.org/10.21334/npolar.2014.645336c7; UNIS Hydrographic Database, https://doi.org/10.21334/unis-hydrography (CC BY 4.0).
