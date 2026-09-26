# Vegetation and delta formation: experiment design

This project asks how hydraulic resistance from sparse, small vegetation changes sediment deposition and channel development in an idealised delta. It is an exploratory extension of a small second-year teaching practical, whose inverted plastic trees acted as qualitative vegetation proxies. A few-percent response can be relevant at that scale; the model does not impose a minimum effect size. It follows the sequence of the 2022 EM4 experiments described in *Effects of the Silurian greening on delta structure and delta formation*. The numerical geometry is an analogue of that experiment. It is not a calibrated reconstruction of the flume or a quantitative reconstruction of a Silurian landscape.

## Evidence and selected conditions

| Quantity | Evidence or status | Treatment |
|---|---|---|
| Water discharge | Essay reports 35 mL/s | 0.000035 m³/s; this is discharge, not velocity |
| Table tilt | Essay reports 2° | Distinguish table inclination from the local sediment-surface slope |
| Bare phase | Essay reports 30 minutes | 1,800 s |
| Half-delta treatment | Vegetation added to half of the existing delta for another 30 minutes | Compare with a bare continuation from the same state |
| River-margin treatment | Reset bed; river-margin vegetation for 30 minutes, then add delta vegetation for 30 minutes | Compare with a matching fresh-bed bare run |
| Box dimensions | Current EM4 manufacturer specification: 3.7 × 1.2 m | Nominal geometry; not a surveyed 2022 working area |
| Sediment density and sizes | Current manufacturer specifications; exact 2022 mixture unknown | Declare the selected mix and fractions in the configuration |
| Bed shape, basin stage, sediment supply and inlet width | Not measured in the essay | Scenario inputs, not inferred measurements |
| Vegetation resistance and placement | Plastic pieces described; no drag or frontal-area measurements | Prescribed spatial treatment and sensitivity range |

The original experiment included manual channel redirection and damming. The numerical comparisons should allow free channel development unless an intervention is explicitly prescribed in both members of a pair. Visual agreement cannot recover an undocumented intervention history.

The plastic media must be represented using its physical grain sizes and density. A model colour does not provide a universal, reversed mapping to natural sediment size. Matching geometry alone also does not establish dynamic similarity: water depth, Froude number and sediment mobility matter.

The input generator selects one 0.7 mm fraction at 1,600 kg/m³, a 0.06 m finite sediment layer and 0.4 porosity. Upstream erosion supplies sediment; no measured external sediment feed is available or imposed. It selects explicit horizontal viscosity and diffusivity of 10⁻⁵ m²/s as uncalibrated flume-scale scenarios, avoiding the engine's much larger general-purpose defaults. The default D-Morphology Van Rijn transport closure is not calibrated for angular plastic. Grain size, mixing and transport assumptions require sensitivity checks before quantitative interpretation.

Clear-water supply is prescribed explicitly rather than assumed from omitted inputs: `NeuBcSand=false` disables equilibrium suspended-sand inflow, and `sedfracbndplastic=0 kg/m³` supplies zero concentration on inflow. At the upstream boundary, `IBedCond=5` and a zero `transport excl pores plastic` table prescribe zero incoming bedload. These settings are necessary because native defaults can introduce external sediment. The native sediment ledger must confirm the realised influx before describing an executed case as a finite-bed experiment.

## Physical structure

For water depth \(h=\eta-z_b\), depth-averaged velocity \(\mathbf u\), water surface \(\eta\), and bed elevation \(z_b\), the depth-averaged balance is schematically

\[
\partial_t h+\nabla\!\cdot(h\mathbf u)=s_w,
\]
\[
\partial_t(h\mathbf u)+\nabla\!\cdot(h\mathbf u\otimes\mathbf u)
=-gh\nabla\eta-\frac{\boldsymbol\tau_b}{\rho}
-\frac{\mathbf F_v}{\rho}+\nabla\!\cdot(h\boldsymbol\nu_t\nabla\mathbf u).
\]

Here \(s_w\) is a water source per horizontal area, \(\boldsymbol\tau_b\) is bed friction, and \(\mathbf F_v\) is vegetation force per horizontal area. These equations describe the physical balances; the retained engine configuration determines the actual discretisation and closures.

An idealised rigid-stem drag law is

\[
\mathbf F_v=\tfrac12\rho C_D a h_v^*|\mathbf u|\mathbf u,
\qquad a=N d,\qquad h_v^*=\min(h,H_v),
\]

where \(N\) is stem count per square metre, \(d\) is stem diameter, \(H_v\) is plant height and \(C_D\) is a drag coefficient. Thus \(a\) has units m⁻¹. This equation explains the drag mechanism conceptually; it is not an exact transcription of the engine's resistance closure. The generated setup selects the native Baptist vegetation model using spatial stem-height, stem-diameter and stem-density fields. Uniform background bed roughness is held fixed between cases. Roughness, stem drag, root reinforcement and sediment cohesion are not interchangeable quantities. Increased resistance may slow flow inside a patch while accelerating or diverting it around the patch.

For sediment fraction \(i\), with depth-averaged suspended volume concentration \(c_i\), settling/deposition \(D_i\), erosion \(E_i\), and suspended diffusivity \(\epsilon_s\),

\[
\partial_t(hc_i)+\nabla\!\cdot(h\mathbf u c_i-h\epsilon_s\nabla c_i)=E_i-D_i.
\]

With porosity \(p\) and volumetric bed-load flux \(\mathbf q_{b,i}\), bed evolution obeys

\[
(1-p)\partial_tz_b+\nabla\!\cdot\sum_i\mathbf q_{b,i}
=\sum_i(D_i-E_i).
\]

Here erosion/deposition are solid-volume fluxes per unit area; mass-based engine outputs require division by the relevant grain density. Boundary transport and changes in suspended inventory must be included in the sediment budget. Morphological acceleration should be one for a direct minute-scale comparison, unless a separate sensitivity test establishes a defensible alternative.

The generated setup couples `Sedimentmodelnr=4` to cell-centred bed elevations (`BedlevType=1`, `BedlevMode=1`, `Conveyance2D=-1`). This is consequential: the pinned community source rejects other bed-level types for this sediment model because subsequent bed-level reconstruction would undo the morphological update. `Kmx=0` selects depth-averaged flow and `AngLat=0` omits Coriolis, appropriate to this laboratory-scale problem.

`BedlevMode=1` and sediment map output are native defaults, confirmed in the pinned source `m_flowparameters.f90` (`BLMODE_DFM=1`, `ibedlevmode=BLMODE_DFM`, map output `sed=1`). Their optional explicit keys are omitted because HYDROLIB 1.0.1 does not expose them; this does not change the selected native behaviour. HYDROLIB parsing checks its supported input schema, not every morphological closure or successful execution. The native run remains a separate gate.

Static XYZ fields use legacy `FILETYPE=7`, `METHOD=5` (triangulation), as implemented in the pinned `meteo1.f90` routine `timespaceinitialfield`. Method 4 invokes the polygon reader for these fields and is inappropriate for XYZ triples. This source-specific distinction matters even when a generic forcing-method table appears to assign different method numbers.

A whole-domain `massbalanceareaflume` sample field requests native mass-balance accounting. `MbaInterval` follows the map-output interval, with `MbaWriteNetCDF=1` and `MbaWriteCsv=1`. These three native keys are read by `unstruc_model.f90` but are absent from HYDROLIB 1.0.1. The test suite uses an explicit three-field schema extension while retaining strict validation of the actual input files; it does not remove these scientifically necessary options to satisfy the older library. The resulting ledger still needs native execution and interpretation of its inventory and exchange terms.

## Establishing the initial flow

The default run first establishes flow for 60 s on a fixed bed, with sediment transport, bed composition and bed-elevation updates all disabled. The active morphological duration is then 1,800 s, giving an absolute stop time of 1,860 s. The native `SedTransStt`, `CmpUpdStt` and `MorStt` offsets use seconds because `Tunit=S`. They are zero in restarted continuation stages; the initial conditioning is not repeated.

This defines a hydrodynamically conditioned numerical experiment, rather than reconstructing the undocumented wetting transient in the original lab experiment. It is a consequential initial-condition choice. In a 60 s diagnostic at a maximum timestep of 0.025 s, a dry start lost about 0.00603 kg of sediment, while the conditioned start's residual was about 0.00000108 kg. This supports a startup-transient explanation without identifying a unique source-code mechanism. The dry-start result is not used as evidence for a vegetation effect.

The raw native water-volume ledger omits the geometric change caused by morphological bed updates. Wet-cell water-surface elevation is generally retained during the operator-split bed update. Accordingly, raw water residual and bed-volume change are reported separately. Explaining their relationship does not establish physical water conservation for a moving bed.

## Comparisons and interpretation

The delivered 30-minute comparison starts four cases from the same fresh bed: bare, half-vegetated delta, its mirror, and river-margin plus delta vegetation. The mirrored patch tests whether the response follows the imposed side rather than an initial geometric preference.

The separate staged workflow compares a bare continuation against vegetation additions from the same saved delta state, closer to the sequence in the practical. Its restart continuity and six stages were checked at 60 seconds of active morphology per stage; those short checks are not presented as a completed 30-minute staged experiment. Compare equal physical durations and identical sediment supply; for supply sensitivity also compare equal cumulative sediment input.

Keep stage, discharge, sediment, initial bed, numerical resolution and output times matched. Vary vegetation resistance separately. Small initial-bed perturbations or repeated patterns help establish whether a difference persists beyond one channel configuration. A changed channel path in one realisation is not sufficient evidence of a robust vegetation response.

Report deposition volume, exported sediment, mean and maximum aggradation, delta footprint at a stated elevation threshold, deposition centroid and flow partition between patches. For channel behaviour use connected active-flow masks, counts at fixed transects and persistence between outputs. State depth/velocity thresholds and test their sensitivity. Short-lived shallow flow paths should not automatically count as persistent distributary channels; do not label every mask change an avulsion.

The essay suggests greater local aggradation, less lateral spreading and more stable channels with vegetation. Treat those as hypotheses. Dense vegetation can also deflect water and sediment into neighbouring channels and reduce deposition inside a patch. The literature therefore does not require a monotonic increase in deposition with plant resistance.

This initial treatment does **not** calculate mechanical root reinforcement. That extension needs an independently specified bank-erosion or sediment-stability law and separate drag-only, reinforcement-only and combined experiments. Fixed vegetation is appropriate for plastic inserts; biological growth and mortality would require additional assumptions.

## Checks required before interpreting a run

1. Confirm that the actual Delft3D FM engine runs the intended inputs, writes evolving hydrodynamic and morphological output, and reports the intended sediment and resistance settings. Input generation alone is not an executed simulation.
2. Check water conservation and the full sediment budget, accounting for initial bed inventory, boundary exchange, suspended material and bed change. Report residuals and normalisation; unexpected residuals require investigation.
3. Check finite values, non-negative depths and concentrations within documented numerical tolerances, wetting/drying behaviour, and boundary-stage/discharge consistency. Inspect velocity, Froude number and sediment mobility at flume scale.
4. Repeat a representative paired comparison with smaller cells and tighter time stepping. Quantify changes in the selected metrics. Do not reduce physical detail solely to obtain a passing run.
5. Separate numerical checks from empirical validation. No measured bed or transport time series are available from the essay, so the available comparison is qualitative.

Silurian interpretation is conditional: this is an early-land-plant-inspired resistance experiment, not a model of a Silurian forest. Neither modern marsh vegetation nor plastic branches constrain Cooksonia rooting depth, strength or drag. A result establishes how the chosen physical treatment behaves under the imposed conditions; it does not identify a unique ancient vegetation assemblage.

## Early-time numerical sensitivity

At 60 s of active morphology after the 60 s flow conditioning, reducing the maximum timestep from 0.025 to 0.0125 s changed positive and negative bed-change volumes by less than 0.03% in the bare and half-vegetated cases. Their 1 mm deposition footprints were unchanged. Reducing cell width from 2 to 1 cm changed those volumes by about 3.7% and the footprint by about 1.1%. The nominal inlet also intersects six coarse boundary cells versus ten fine cells, so its resolved width contributes to this spatial sensitivity.

Both grids retained the same direction of the small early vegetation contrast. This is a sensitivity check at one early time, not a full-duration grid-convergence study. The JSON records in `verification/` preserve individual metrics, so a small treatment effect can be assessed separately from changes in absolute totals.
