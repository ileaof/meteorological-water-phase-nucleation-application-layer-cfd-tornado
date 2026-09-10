# Saved-output tornado audit — 2026-09-05

The saved fields support a low-level concentration/alignment problem, with different manifestations in different runs. They do **not** support the simpler claim that the solver cannot generate rotation or stretching. The mature parent contains both. A separate, verified implementation limitation is that the low-memory projection's dynamic pressure does not enter saturation thermodynamics. Therefore absence of a visible condensation funnel is not an independent test of whether the velocity field contains a tornado.

No model, solver, configuration, or original output was changed. No new simulation was run. This audit added only `scripts/audit_saved_tornado.py`, this report, and generated diagnostics under `outputs/implementation_audit/`.

## Scope and reproducibility

Independently processed all 27 NPZ files containing u, v, w, including duplicate render exports and smoke runs; these are **not 27 independent experiments**. All processed successfully. Other NPZ files contain lookup/radar data. The NetCDF inventory consists of reference nucleation outputs, initial conditions, boundaries, and forcing; these were not treated as mature storm evolution.

Run from the repository using Python with NumPy and Matplotlib:

```text
python scripts/audit_saved_tornado.py
```

Outputs:

- `outputs/implementation_audit/audit.json`: input keys, geometry, pressure availability, reconstructed statistics, vertical profiles, boundary sensitivity, and downward ridge candidates.
- One CSV per velocity snapshot: every height, all three vorticity components, stretching, tilting, convergence, dilatation, rotation/updraft alignment, and pressure anomaly where saved.
- `parent_vertical_profiles.png`: inspected 0–2 km parent profiles.
- `legacy_alignment_history.csv`: existing evolution diagnostics transcribed from their JSON, not independently reconstructed in time.

The independent derivative routines pass analytic rotating/converging-flow and tilting-sign checks on a nonuniform vertical grid. Saved zeta agrees with independently reconstructed interior zeta to roundoff (maximum errors below 4e-16 s^-1). This establishes that the interior curl arithmetic in those exports is not the principal problem.

## Definitions and limits

Velocities are averaged from C-grid faces to centres when staggered. Derivatives use actual x, y, z coordinates and NumPy's second-order nonuniform-grid formulas. Some rendering exports store coordinates in kilometres but dx in metres; their ratio determines the conversion. Parent caches omit geometry; dx=dy=600 m, Lz=15 km and geometric stretch=1.05 are reconstructed from the producing scripts. Such caches lack sufficient provenance to verify every historical environment override or source revision.

With x eastward, y northward, and z upward:

```text
xi   = dw/dy - dv/dz
eta  = du/dz - dw/dx
zeta = dv/dx - du/dy                         [s^-1]
C    = -(du/dx + dv/dy)                     [s^-1; positive convergence]
S    = zeta * dw/dz                         [s^-2]
T    = xi * dw/dx + eta * dw/dy              [s^-2]
D    = -zeta * (du/dx + dv/dy + dw/dz)       [s^-2]
S+D  = zeta * C                             [s^-2]
```

These are relative-vorticity kinematic terms, with no assumed Coriolis parameter. They are **not a closed prognostic budget**: time tendency, advection, actual stress curl, boundary forcing, and the solver-consistent pressure/buoyancy curl would be needed for attribution. In anelastic flow, velocity divergence is not generally zero; stretching alone must not be called net concentration.

Statistics exclude the outer one-sixth in each horizontal direction, with additional 5%, 10%, 20%, and 25% mask checks. Those fractions are sensitivity tests, not established physical sponge widths, and cannot define a controlled cross-resolution comparison. No derivative at the lateral outer edge is used for the interior findings. Centred-velocity divergence is a diagnostic, not the native face-based mass-projection residual.

The main profile independently selects each level's strongest **positive** zeta and strongest updraft. Their offset is a warning about alignment, not proof they belong to the same storm or vortex. A second diagnostic follows a local cyclonic ridge downward from approximately 1.5 km within a fixed 1 km search radius per level. It is a geometric candidate, not a materially tracked vortex: permitted slope depends on vertical spacing, and competing peaks can still switch. No topology or persistence classification is asserted.

## Measured mature parent

Input: `outputs/parent_matured_120_48_2800.npz`, approximately t=2800 s, 600 m horizontal spacing, first cell centre 39.89 m, 17 centres below 2 km.

| Height | Cyclonic zeta maximum | w at that maximum | Offset from strongest updraft | Stretching | Tilting |
|---|---:|---:|---:|---:|---:|
| 39.9 m | 0.00555 s^-1 | 0.302 m/s | 1.34 km | +4.22e-5 s^-2 | -1.05e-6 s^-2 |
| 297.7 m | 0.00533 | 2.52 | 2.40 km | +4.35e-5 | -1.10e-5 |
| 705.7 m | 0.00976 | 4.62 | 5.53 km | +3.88e-5 | +1.11e-5 |
| 1488.3 m | 0.01076 | 14.40 | 0.60 km | +5.61e-5 | +4.37e-5 |

At the lowest cyclonic maximum C=+0.00711 s^-1, so net kinematic concentration is approximately +3.95e-5 s^-2. Thus **missing convergence or missing stretching everywhere is contradicted by this snapshot**. Low-level tilting at the selected core is negative through roughly 300 m. Stronger rotating ascent is elevated, and the independently selected maxima separate substantially below it.

The downward ridge from 1488 m reaches 39.9 m with zeta=0.00521 s^-1 and w=0.340 m/s; its endpoints shift about 2.47 km horizontally. This is compatible with a sloping, weak surface extension of a mesocyclone. It does not demonstrate an intense, narrow, vertically aligned tornado core, but it also rules out claiming that all rotation is disconnected from the lowest resolved layer.

A 100–300 m-wide feature spans only 0.17–0.5 horizontal cells on this parent. Its visible funnel geometry cannot be resolved by the parent output. First resolved cloud condensate (ql+qi >1e-5 kg/kg) is at 705.7 m anywhere in the domain; no cloud condensate meeting that threshold reaches the first level. This is a field measurement, not a diagnosis of the reason for cloud-base height.

## Fine nests: distinct failure modes

**Moore funnel export:** `outputs/moore_real_funnel/fields_L3.npz`, dx=45.74 m, first centre 62.09 m, 12 levels below 2 km. The all-domain low-level |zeta| maximum is 0.02768 s^-1, versus 0.002352 s^-1 inside the one-sixth mask: a factor of 11.8. At 10%, 20%, and 25% exclusions it is 0.00259, 0.00229, and 0.00206 s^-1. Large exterior peaks are not representative of the weak central circulation.

At 62 m the cyclonic maximum is only 0.000576 s^-1, with w=0.191 m/s, stretching=1.81e-6 s^-2, and tilting=-6.70e-7 s^-2. Around 463 m, w at the maximum is 0.185 m/s and stretching only 5.35e-7 s^-2. This snapshot supports weak low-level rotating ascent/concentration, despite a fine horizontal mesh. Its first centre and roughly 124 m first cell thickness also show that fine horizontal resolution is not fine surface-layer resolution. This export contains neither pressure nor condensate, so it cannot establish a pressure deficit or condensation funnel.

**Matched-domain v2 fine export:** `outputs/matched_domain_v2/fields_fine.npz`, dx=22.22 m, first centre 5.05 m, 38 levels below 2 km. At 5 m its strongest interior cyclonic rotation has zeta=0.01011 s^-1, w=-0.00716 m/s, C=-0.00159 s^-1, S=-1.43e-5 s^-2 and T=-1.82e-6 s^-2. Around 98 m, w=-0.185 m/s, C=-0.00439 s^-1 and S=-5.63e-5 s^-2. The selected near-ground circulation is descending and diverging, so it is not undergoing the convergence/stretching needed for intensification at that location and instant.

The peak switches to ascending air near 283 m, with offset falling from about 2.7 km to 0.37 km. This jump must not be interpreted as a single coherent column. Boundary sensitivity is still substantial: low-layer |zeta| falls from 0.206 at 5% exclusion to 0.0592 at 10% and 0.0193 at 25%. No tornado intensity is assigned from these maxima.

These two nests represent different experiments and times. Their numbers do not constitute a resolution-convergence result or a causal intervention.

## Pressure and rendering audit

Current `src/storm_dynamics/core.py:93` stores low-memory projection pressure in `st.p_dyn=phi/dt`. `src/meteorological_flow/state.py:91` uses `P_total=P_base+self.p` for thermodynamics. The dynamic field is deliberately separate because boundary-driven nest pressure previously destabilized thermodynamics; this audit does not change that decision.

The mature parent cache confirms this distinction: **p is identically zero**, while p_dyn ranges from -570.3 to +338.7 Pa. Dynamic pressure exists in the momentum projection, but its perturbation is absent from the saved thermodynamic pressure path. This is a concrete limitation on representing the pressure-induced lowering of saturation height in the low-memory formulation, not proof that fixing the coupling would produce a tornado.

Using the interior plane median as a reproducible pressure reference, p_dyn at the selected parent cyclonic peak is +69.4 Pa near 40 m and -284.9 Pa near 1.49 km. These values sample different selected locations and contain general storm pressure, not an isolated vortex pressure decomposition. They do not show a co-located near-surface pressure minimum at the selected low-level peak.

The matched-domain v2 producer explicitly saves p_dyn under the name `p`; its range is -5790 to +7695 Pa. Existing boundary/projection concerns make these unsuitable for a physical tornado-pressure threshold. The audit reports raw plane-relative anomalies without treating them as vortex suction. Velocity alone at one instant cannot uniquely recover physical perturbation pressure without the momentum forcing, density/base state, boundary conditions and appropriate tendency information.

`examples/render_tornado_3d.py` sums ql, qi, qr, qs, qg and qh into `cond`. This includes precipitation, so an isosurface of that field is not specifically a condensation funnel. `scratchpad/render_storm_3d.py` uses ql+qi but reconstructs missing grid metadata and uses a default LCL; its graphical output alone is not a dynamics validation.

Physical interpretation: tornado-scale concentration and its pressure deficit are distinct requirements beyond an organized convective cloud, and a visible funnel additionally depends on saturation. Relevant synthesis: [Supercell Tornadogenesis: Recent Progress in Our State of Understanding](https://journals.ametsoc.org/view/journals/bams/105/7/BAMS-D-23-0031.1.xml). Near-ground tilting/stretching attribution is also sensitive to parcel history and surface processes: [Transition of Near-Ground Vorticity Dynamics during Tornadogenesis](https://journals.ametsoc.org/view/journals/atsc/79/2/JAS-D-21-0181.1.xml).

## Existing diagnostics and time evidence

Current `surface_connection_report` selects a separate |zeta| maximum at each height and classifies connection largely by a surface/aloft speed ratio. It does not require spatial continuity of those maxima. Its reported rotational speed is the maximum departure from a local mean horizontal wind, which also responds to shear and strain; it is not exclusively azimuthal circulation. The previously documented minimum-cell radius inflation has now been fixed in source by rejecting under-resolved radii. The older REVIEW_REQUEST is therefore not an authoritative description of every current defect.

`vorticity_budget.py` labels its equation D(zeta)/Dt while also including explicit advection on the RHS; that sum corresponds to an Eulerian partial time tendency. Its tilting implementation has the correct curl-consistent sign. Missing pressure/density/stress inputs become zero contributions, which should be interpreted as unavailable terms, not measured physical absence. Its scalar-diffusion expression is explicitly a proxy rather than the exact curl of SGS momentum stresses.

The budget module claims periodic-aware horizontal gradients, but `Grid._central_x/_central_y` use one-sided edges, including on periodic grids. Its stretched-z central difference uses a two-point secant across neighbours, rather than the general second-order three-point formula. These affect diagnostic accuracy, especially at boundaries; the independent audit avoids using those operators. No evidence here establishes either discrepancy as the cause of failed tornadogenesis.

The existing freely evolving supercell history reports signed tilting alignment rising from 0.020 at 301 s to 0.205 at 2400 s, followed by a low-level |zeta| peak of 0.00953 s^-1 at 3300 s and decay to 0.00248 at 5400 s. Streamwise fraction reaches 0.638 at the end. This is evidence of transient organization, not sustained tornadic intensification. It is a separate run, and a domain-wide signed alignment can hide cancellation of positive and negative local production. Matching 3-D time sequences were not saved, so those values cannot be independently converted into trajectories, closed budgets, or time-resolved vortex coherence.

## Diagnosis and remaining evidence gap

The organized anvil is compatible with deep convection while the low-level requirements remain unmet. In the parent, rotation and positive concentration exist but are broad, sloping and poorly aligned with the strongest low-level updraft; the mesh cannot represent narrow funnel geometry. In the weak Moore nest, near-ground rotating ascent is too weak in the saved interior. In the matched fine snapshot, the selected surface cyclonic maximum is locally descending and diverging. Boundary peaks and existing surface-connection metrics can overstate how tornado-like these fields are. Separately, dynamic pressure is excluded from the low-memory thermodynamic saturation path.

These observations diagnose the saved states and identify a verified thermodynamic coupling limitation. They do **not** uniquely attribute failure to numerical diffusion, cold-pool strength, surface drag, or nesting, and do not prove that any particular solver modification would fix it.

Completing causal attribution requires a synchronized sequence from one identified run: native face u/v/w, coordinates and nest origins, time/dt, p and p_dyn with semantics, base density/pressure, theta/qv/ql/qi, buoyancy, momentum tendencies before/after advection/SGS/drag/projection/relaxation, and exact configuration/source revision. That would permit material tracking, native mass-residual checks and budget closure in the same physical low-level region. No rerun or solver instrumentation was performed in this audit.
