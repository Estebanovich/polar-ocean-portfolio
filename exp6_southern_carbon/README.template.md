# Experiment 6 — Does anthropogenic carbon weaken the wind sensitivity of the Southern Ocean CO₂ sink? An online test in an idealized channel

<!-- Rendered by scripts/render_readme.py: every number below comes from figures/*.json, cases.yaml
or the analysis constants. Edit README.template.md, never README.md. -->

## Abstract

Stronger westerlies over the Southern Ocean bring more carbon-rich deep water to the surface and, in the pre-industrial ocean, release CO₂. Jarníková et al. (2025) argue, from an offline decomposition of an ocean biogeochemical model, that this wind response weakens as anthropogenic carbon accumulates near the surface and erodes the vertical gradient of dissolved inorganic carbon (DIC). We test that inference online, in the idealized re-entrant channel of the MITgcm tutorial ({{dx_coarse_km}} km, Gent–McWilliams eddies) coupled to `pkg/dic`. A control run at the {{y_first}} atmospheric CO₂ (GC) and a transient run under historical and SSP3-7.0 CO₂ (GT) share a spin-up of {{g0_years}} tracer years. From both, the zonal wind stress is stepped up by a factor of {{wind_factor}} in {{y1}}, {{y2}}, {{y3}} and {{y4}}, and the sensitivity S is the change of the air–sea flux over branch years {{win_lo}}–{{win_hi}}. Between the {{y1}} and {{y4}} branches the upper-200 m DIC gradient in the upwelling region falls by {{g200_drop}} %. The wind sensitivity of the transient run changes by ΔS = {{dS}} mol C m⁻² yr⁻¹, against a control drift of {{drift_gc}}, so the registered hypothesis H6a {{h6a_verdict}}. The offline decomposition of Jarníková et al. (eqs. 1 and 4, with a phosphate-based biological term) reproduces the measured change of the mixed-layer carbon budget within ±25 % in {{h6b_n}} of the four branches (H6b {{h6b_verdict}}). In a {{dx_fine_km}}-km eddy-permitting version, H6c {{h6c_verdict}}.

## 1. Question and design

The Southern Ocean takes up a large share of the anthropogenic CO₂ that enters the ocean, and its sink responds to the wind (Le Quéré et al. 2007; Lovenduski et al. 2008; Gruber et al. 2018). Stronger westerlies drive more upwelling of old, DIC-rich water: the natural carbon outgasses, while the anthropogenic carbon is taken up faster by the newly exposed water (Lovenduski et al. 2008). Jarníková et al. (2025) decompose this response in a global model into mixed-layer entrainment (their eq. 1), upwelling by the overturning (eq. 4) and the biological response (eqs. 2, 3 and 5), and infer that the outgassing term shrinks as the vertical DIC gradient is eroded by anthropogenic carbon. The inference is offline: the gradients of one state are combined with the circulation change of another.

We ask whether the same physics, run online, gives the same answer. Three hypotheses were registered before the transient run started (`notes/pd0_preregistro.md`):

- **H6a.** The wind sensitivity of the transient run moves toward uptake between the {{y1}} and {{y4}} branches, ΔS = S_GT({{y4}}) − S_GT({{y1}}) > 0, and by more than twice the drift of the control, |S_GC({{y4}}) − S_GC({{y1}})|.
- **H6b.** In each of the four transient branches, the offline sum (eq. 1 + eq. 4 + a biological term with phosphate in place of iron) is within ±25 % of the measured change of advection, mixing and biology in the mixed-layer DIC budget of the upwelling region.
- **H6c.** In the eddy-permitting model, the same ΔS is detectable (|ΔS| > 2σ over the members), has the sign of the coarse ΔS, and is within a factor [0.5, 2] of it.

Validation tests V6a–V6d have thresholds fixed in the same document; a failure is reported as a failure, and failures are classified as (a) our error, (b) missing physics, (c) a numerical choice or (d) an error in the reference.

## 2. Configuration

The channel is `tutorial_reentrant_channel` of MITgcm (Fig. 0): {{Lx_km}} km periodic in x, {{Ly_km}} km closed in y, {{nr}} levels, with the tutorial's zonal wind stress, its sloping southern boundary and its restoring of surface temperature. We add a seasonal cycle to the restoring temperature (Fig. 0a) and to the light, so that the winter mixed layer deepens and the entrainment term of Jarníková et al. has something to act on. Carbon is carried by `pkg/dic` with six tracers (DIC, alkalinity, phosphate, dissolved organic phosphorus, oxygen and preformed phosphate); there is no iron, no explicit ecosystem and no sea ice. Gas exchange uses a fixed wind speed, so a change of the wind stress acts on the carbon only through the circulation. The northern rows are a sponge that restores the tracers to pre-industrial profiles. The analysis region R is the southern {{R_rows}} rows (y < {{R_edge_km}} km), where the upper overturning cell upwells.

![Fig. 0](figures/fig0_schematic.png)

**Fig. 0.** Geometry and forcing. (a) Zonal wind stress and the restoring temperature with its seasonal range. (b) Depth, the region R and the pre-industrial sponge.

The coarse model ({{dx_coarse_km}} km, time step {{dt_s}} s) parameterizes eddies with Gent–McWilliams and Redi. Its spin-up (G0) accelerates the tracers by a factor of {{accel}} against the momentum (Bryan 1984) for {{g0_years}} tracer years, followed by {{g0s_years}} synchronous years (G0s), about the length Wang (2001) found sufficient after an asynchronous equilibrium. The drift in the last century of G0 is {{v6a_dflux}} mol C m⁻² yr⁻¹ in the channel-mean flux and {{v6a_ddic}} µmol kg⁻¹ in the DIC below the thermocline, so V6a {{v6a_verdict}}; over G0s it is {{g0s_dflux}} and {{g0s_ddic}}. At equilibrium the upper cell carries {{moc_R}} Sv out of R ({{moc_edge_frac}} % of its maximum, reached at the edge of R), the winter mixed layer in R reaches {{mld_R}} m in {{mld_month}}, and the channel takes up {{flux_g0}} mol C m⁻² yr⁻¹ at the pre-industrial CO₂. The share of the cell's growth that lies inside R is {{upw_frac}} %; the rest is noise on the plateau north of R, not upwelling.

GC holds the atmospheric CO₂ at its {{y_first}} value for {{gt_years}} years; GT follows the CMIP6 historical concentrations (Meinshausen et al. 2017) and then SSP3-7.0 (Meinshausen et al. 2020) to {{y_last}}. Each branch starts from the parent's state on 1 January of its year, multiplies the wind stress by {{wind_factor}} and runs {{branch_len}} years. The budget of the annual channel DIC closes to {{v6b_max}} of the surface flux (V6b {{v6b_verdict}}).

## 3. Method

**Sensitivity.** F is the channel-mean `DICCFLX` in mol C m⁻² yr⁻¹, positive into the ocean, and S(t_b) is the mean of F_branch − F_parent over branch years {{win_lo}}–{{win_hi}}. Following Lovenduski et al. (2008), S_nat = S_GC and S_ant = S_GT − S_GC.

**Gradients.** G₂₀₀ = [DIC(200 m) − DIC(surface)]/200 m and ΔDIC₁₀₀₀ = DIC(1000 m) − DIC(surface), over R, in the decade of the parent centred on the branch year (Jarníková et al. 2025, eqs. 1 and 4).

**Offline decomposition.** With the parent's gradients and the branch-minus-parent changes of the winter mixed-layer depth δMLD and of the upwelling δMOC (residual overturning in temperature layers, including the GM bolus), the terms in mol C yr⁻¹ are δMLD × G₂₀₀ × ρ₀ × A_R (eq. 1), δMOC × ΔDIC₁₀₀₀ × ρ₀ (eq. 4), and −R_CP × e × (the extra phosphate supply by the same two routes), with R_CP = {{r_cp}} and e the export fraction of the supplied phosphate in the parent. The measured counterpart is the branch-minus-parent change of advection + mixing + biology in the DIC budget of the mixed layer of R, from the model's own tendency diagnostics. A temperature term, the change of the flux from the branch-minus-parent surface temperature alone (carbonate chemistry of `pkg/dic`), acts through gas exchange rather than circulation; it is reported but not added to the sum.

## 4. Results

### 4.1 The gradient erodes

![Fig. 1](figures/fig1_sections.png)

**Fig. 1.** DIC in GT, zonal mean, in the decades centred on (a) {{y1}} and (b) {{y4}}; (c) G₂₀₀ by row, with R shaded.

Anthropogenic carbon enters at the surface and raises the DIC of the upper ocean, while the deep water that upwells in R changes little over the run, so the vertical contrast shrinks (Fig. 1). Over R, G₂₀₀ falls from {{g200_first}} to {{g200_last}} µmol kg⁻¹ m⁻¹ between the {{y1}} and {{y4}} decades.

![Fig. 2](figures/fig2_series.png)

**Fig. 2.** Channel-mean air–sea flux of GC, GT and the wind branches (dashed).

### 4.2 H6a: the wind sensitivity

![Fig. 3](figures/fig3_central.png)

**Fig. 3.** Wind sensitivity of each branch against the parent's G₂₀₀ in R: S_GT, S_nat = S_GC and S_ant = S_GT − S_GC.

| Branch | G₂₀₀ | S_GT | S_nat | S_ant |
|---|---|---|---|---|
| {{y1}} | {{g200_1}} | {{s_gt_1}} | {{s_nat_1}} | {{s_ant_1}} |
| {{y2}} | {{g200_2}} | {{s_gt_2}} | {{s_nat_2}} | {{s_ant_2}} |
| {{y3}} | {{g200_3}} | {{s_gt_3}} | {{s_nat_3}} | {{s_ant_3}} |
| {{y4}} | {{g200_4}} | {{s_gt_4}} | {{s_nat_4}} | {{s_ant_4}} |

**Table 1.** G₂₀₀ in µmol kg⁻¹ m⁻¹; S in mol C m⁻² yr⁻¹, positive = more uptake under the stronger wind.

ΔS = {{dS}} mol C m⁻² yr⁻¹ ({{dS_pct}} % of |S_GT({{y1}})|), against a control drift of {{drift_gc}}: H6a {{h6a_verdict}}. S_GT {{zero_cross}}.

The whole change is in S_ant. S_nat, the response of the control, is the same in every branch, as it should be: the control has the same circulation and the same atmospheric CO₂ in every branch year. Under the stronger wind the transient channel loses less carbon than the control from the first branch on, and in the last branch it takes up more than it did before the wind step. The direction of the change is the one Jarníková et al. (2025) infer; in this channel it is large enough to reverse the sign of the response.

### 4.3 H6b: the offline decomposition

![Fig. 4](figures/fig4_h6b.png)

**Fig. 4.** Offline terms against the measured change of the mixed-layer budget of R, per transient branch. The temperature term (hatched) is not part of the sum.

| Branch | eq. 1 (MLD) | eq. 4 (MOC) | bio (PO₄) | measured | offline / measured | temperature |
|---|---|---|---|---|---|---|
| {{y1}} | {{mld_1}} | {{moc_1}} | {{bio_1}} | {{meas_1}} | {{ratio_1}} | {{temp_1}} |
| {{y2}} | {{mld_2}} | {{moc_2}} | {{bio_2}} | {{meas_2}} | {{ratio_2}} | {{temp_2}} |
| {{y3}} | {{mld_3}} | {{moc_3}} | {{bio_3}} | {{meas_3}} | {{ratio_3}} | {{temp_3}} |
| {{y4}} | {{mld_4}} | {{moc_4}} | {{bio_4}} | {{meas_4}} | {{ratio_4}} | {{temp_4}} |

**Table 2.** Branch − parent, Gmol C yr⁻¹, positive = DIC added to the mixed layer of R.

The offline sum is within ±25 % of the measured change in {{h6b_n}} of the four branches: H6b {{h6b_verdict}}.

**Post hoc (not registered).** Compared term by term (Table 3), the physical and the biological parts fail in different ways.

| Branch | offline physics (eq. 1 + eq. 4) | measured advection + mixing | offline / measured | measured mixing | measured bio |
|---|---|---|---|---|---|
| {{y1}} | {{phys_off_1}} | {{phys_meas_1}} | {{phys_ratio_1}} | {{mix_meas_1}} | {{bio_meas_1}} |
| {{y2}} | {{phys_off_2}} | {{phys_meas_2}} | {{phys_ratio_2}} | {{mix_meas_2}} | {{bio_meas_2}} |
| {{y3}} | {{phys_off_3}} | {{phys_meas_3}} | {{phys_ratio_3}} | {{mix_meas_3}} | {{bio_meas_3}} |
| {{y4}} | {{phys_off_4}} | {{phys_meas_4}} | {{phys_ratio_4}} | {{mix_meas_4}} | {{bio_meas_4}} |

**Table 3.** As Table 2, with the measured change split into its physical and biological parts (post hoc).

- **Biology.** The phosphate-based term is the same in every branch, and so is the measured biological change: neither the circulation nor phosphate depends on CO₂. The offline term overestimates it in every branch.
- **Physics.** The offline sum of eqs. 1 and 4 declines as G₂₀₀ erodes, as does the measured supply, which is the mechanism Jarníková et al. propose. But it recovers only part of the measured supply in the early branches, and only part of its decline. Most of the gap is mixing: in the budget it falls from {{mix_meas_1}} to {{mix_meas_4}} Gmol C yr⁻¹ between the {{y1}} and {{y4}} branches, while eq. 1, which sees only the change of mixed-layer depth, gives {{mld_1}} and {{mld_4}}.
- **The sum.** The physical and biological changes have opposite signs and similar sizes, so modest errors in each become large errors in their sum.

The decomposition gets the direction of the change right but not its size. We classify the failure tentatively as (d), a process the offline decomposition does not represent (the change of diapycnal and eddy mixing under the stronger wind). An error of class (a) in our implementation of eq. 1 is not excluded.

### 4.4 H6c: eddies

![Fig. 5](figures/fig5_fine.png)

**Fig. 5.** Wind sensitivity of the {{dx_fine_km}}-km members and of the coarse model in the {{y1}} and {{y4}} branches.

The {{dx_fine_km}}-km model, without GM, starts from the coarse state of GT in {{y1}} and {{y4}} and runs {{n_members}} control/wind pairs of {{fine_len}} years each. Its ΔS is {{h6c_dS}} mol C m⁻² yr⁻¹ against a detection threshold of {{h6c_2sd}}, and {{h6c_ratio}} times the coarse ΔS: H6c {{h6c_verdict}}.

## 5. Limitations

- **Missing physics.** There is no iron, which limits Southern Ocean productivity, no explicit ecosystem (the biology of `pkg/dic` is a nutrient-restoring export, a null hypothesis for the biological term), no sea ice and no ozone forcing. The biological term of H6b therefore uses phosphate where Jarníková et al. use iron.
- **The sponge.** The northern sponge restores to pre-industrial profiles, so anthropogenic carbon enters the channel only through the surface; the lateral supply of anthropogenic carbon from the north is missing. A sensitivity run with an anthropogenic sponge was cut from the scope for time.
- **Resolution.** At {{dx_fine_km}} km the eddies are only partly resolved, and the fine branches are short; H6c tests the sign and the order of magnitude, not the eddy compensation itself.
- **The temperature term.** The temperature term uses an empirical conductance: the parent's air–sea flux over R divided by its mean air–sea pCO₂ difference. In the first branch that difference is small, so the conductance, and with it the temperature term, is poorly conditioned. The term is reported but not used in any verdict.
- **Idealization.** One channel, one wind profile and a fixed gas-exchange wind; the step in wind stress is a sensitivity experiment, not a scenario.

## 6. Reproducibility

`reproduce.sh --check` lists every step; `reproduce.sh --quick` regenerates the evaluation, figures and this README from existing runs. Configuration: `cases.yaml`; inputs: `gendata/`; pre-registration and decisions: `notes/` (in Spanish).

## References

Bryan, K., 1984: Accelerating the convergence to equilibrium of ocean-climate models. *J. Phys. Oceanogr.*, https://doi.org/10.1175/1520-0485(1984)014%3C0666:ATCTEO%3E2.0.CO;2

Gruber, N., et al., 2018: The variable Southern Ocean carbon sink. *Annu. Rev. Mar. Sci.*, https://doi.org/10.1146/annurev-marine-121916-063407

Jarníková, T., et al., 2025: Decreasing importance of carbon-climate feedbacks in the Southern Ocean in a warming climate. *Sci. Adv.*, https://doi.org/10.1126/sciadv.adr3589

Le Quéré, C., et al., 2007: Saturation of the Southern Ocean CO₂ sink due to recent climate change. *Science*, https://doi.org/10.1126/science.1136188

Lovenduski, N. S., et al., 2008: Toward a mechanistic understanding of the decadal trends in the Southern Ocean carbon sink. *Global Biogeochem. Cycles*, https://doi.org/10.1029/2007gb003139

Meinshausen, M., et al., 2017: Historical greenhouse gas concentrations for climate modelling (CMIP6). *Geosci. Model Dev.* <!-- final DOI not verified (Consensus gives only a discussion-comment DOI) -->

Meinshausen, M., et al., 2020: The shared socio-economic pathway (SSP) greenhouse gas concentrations and their extensions to 2500. *Geosci. Model Dev.*, https://doi.org/10.5194/gmd-13-3571-2020

Wang, D., 2001: A note on using the accelerated convergence method in climate models. *Tellus A*, https://doi.org/10.3402/tellusa.v53i1.12179
