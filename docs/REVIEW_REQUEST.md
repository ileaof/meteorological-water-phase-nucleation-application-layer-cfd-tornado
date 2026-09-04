# Tornadogenesis CFD — defect synthesis for external review

**Date:** 2026-09-04
**Codebase:** idealized supercell/tornadogenesis solver (`src/storm_dynamics/`, anelastic, GPU/CuPy,
AMR nesting) + a real-data pipeline (`src/atmospheric_data/`).
**Goal of the study:** get a simulated tornado-like vortex and compare against the
Moore 2013-05-20 KTLX radar observation (V_rot 26 m/s, couplet separation ~253 m, ζ 0.205 s⁻¹).
**Status:** the model produces a surface-connected low-level mesocyclone. It has NOT produced a
tornado. Most quantitative claims made during the study turned out to be diagnostic artifacts.

What follows is (A) numerical defects found and measured, (B) claims that were asserted and later
disproved by measurement, (C) a recurring failure pattern, (D) physics results that survive,
(E) the questions I would most like a second opinion on.

---

## A. Numerical defects found (all measured, not inferred)

### A1. Low-memory pressure projection assumed SOLID WALLS on every non-periodic grid — FIXED
`project_anelastic_fft` did `u[0]=u[-1]=0; v[:,0]=v[:,-1]=0` for any `periodic_h=False`.
Correct for a nest (walls + sponge). Wrong for a limited-area parent, where the boundary normal
velocity carries the environmental inflow: it was zeroed every step, and the BC then restored a
value that had never been part of the projected field.

Scaled `|div(ρu)| / (ρU/dx)`, straddling the 64k-cell low-memory threshold:

| cells | solver | periodic | residual |
|---|---|---|---|
| 55,296 | direct | True | 2.8e−06 |
| 55,296 | direct | False | 1.2e−05 |
| 73,728 | lowmem | True | 5.9e−06 |
| 73,728 | lowmem | **False** | **3.5e−01** |

Fix: `lateral="wall"|"open"`. The `"open"` path keeps the boundary velocity and enforces the
Neumann solvability condition (a uniform normal-velocity offset; the boundary *profile* is
preserved). Limited-area parent went 6.36e−01 → **3.24e−05**, matching the periodic parent.

### A2. Sponge width specified in CELLS, so it shrinks physically at every refinement level — FIXED (opt-in)
`NestSpec.relax_width = 4` cells → 267 m at dx=67 m but **89 m at dx=22 m**. Damping across the
band scales as `rate × width / U`, so a finer nest damps ~`refine`× less while having more error
to absorb. Fix: `relax_width_m` pins a physical width. Default unchanged.

### A3. `follow_spec` dropped the entire sponge configuration on every nest move — FIXED
It rebuilt the moved spec with only `(refine, nz, z_stretch)`, silently reverting
`relax_width / relax_rate / relax_width_m` to defaults on the FIRST re-centring. Visible in logs:
a 267 m band read 89 m after one move. **Consequence:** an entire re-test run (317 min) that was
supposed to validate A2 never actually had A2 active beyond its first seconds.

### A4. Moving-nest tracker was steered by its own boundary artifact — FIXED (opt-in)
`tag_cells` thresholds at `frac × domain max`. When the tracked grid is itself a nest, its
sponge-generated edge vorticity both inflates the normaliser and drags the tagged cluster into the
wall — a feedback pushing the real vortex toward the opposite boundary. Fix: `border=` exclusion,
plumbed as `run_multilevel_nest(follow_border="auto")`.

### A5. Nest sponge applied AFTER the projection — FIXED (opt-in)
Order was `predictor → project → transport(… → relax_to_parent)`, so nudging the border
re-introduced divergence into the relaxation band. Relaxing first:

| region | relax after (default) | relax before |
|---|---|---|
| outermost face | 3.430e−01 | 3.402e−01 (unchanged) |
| cells 1–3 | 4.104e−04 | **1.468e−15** |
| cells 4–11 | 5.298e−05 | **1.265e−15** |
| interior | 7.03e−16 | 9.96e−16 |

### A6. Nest OUTERMOST FACE is still divergent — NOT FIXED
Cause: `_project` calls `_apply_nest_bcs()` *after* the solve, rewriting the faces the projection
just set. A nest legitimately carries net mass flux through its boundary (parent inflow one side,
outflow the other), so the parent's net-zero enforcement (A1) is NOT directly applicable — only
the *discrete* imbalance may be removed. **Open.**

### A7. `surface_connection_report` converts radius to CELLS with a floor
```python
R = max(3, int(radius_m / grid.dx))
```
A radius < 3·dx silently becomes 3 cells — a DIFFERENT physical disk per resolution. A nominal
"400 m" radius means **1800 m at dx=600** and **900 m at dx=300**: a 2× mismatch in exactly the
quantity under test in a resolution comparison. Also affected earlier nest runs
(`int(400/66.7)=5 → 333 m` vs `int(400/22.2)=18 → 400 m`). Worked around by choosing
radius ≥ 3·dx_coarsest; **not fixed in source.**

### A8. `p_dyn` / pressure deficit unusable on nests — known, unresolved
The nest projection absorbs the imposed-inflow imbalance, so `phi` carries a large
boundary-driven component (~1/dt). Measured local values on one 67 m field: −1418, −3250, −5161,
**+7412 Pa**. Any diagnostic keyed on a pressure threshold (e.g. "vortex = ζ peak co-located with
a p_dyn depression") is unreliable on nests. Cross-check against `−ρ v_θ²`.

### A9. Real-data pipeline was never actually limited-area — FIXED
Three independent, silent defects: (1) `_make_grid` built `periodic=True` under a Davies zone;
(2) fixing that was cosmetic because `StormSimulation` builds its OWN grid from
`build_storm_config`; (3) `run_multilevel_real_case` never applied the lateral relaxation at all
(only the single-level `run_case` did). Also, the engine had no open `y` boundary — only
`free_slip`/`wall` (both pin v=0) or `periodic` — so a Davies zone nudging v toward an
environment with v≠0 was fighting a wall.

### A10. ERA5 fetcher hardcoded 11 pressure levels — FIXED (27)
Only two below 1 km. (Note: this turned out NOT to be the SRH limiter — see B3.)

---

## B. Claims I asserted and then DISPROVED by measurement

These are listed because the pattern matters more than any single item: in each case a plausible
mechanism was stated as a finding before it was measured.

**B1. "Periodic laterals mean the storm ingests its own cold pool; there is no environmental
inflow."** Inferred from transit distance (parcel travels 42–50 km in a 72 km box).
**MEASURED FALSE:** on the periodic parent at t=2800, only **1.8 %** of the 500 m level is >1 K
colder than the environment (0.2 % >2 K), and inflow air 5–60 km upstream is unmodified
(dθ −0.23…+0.14 K, dq_v −0.00…+0.06 g/kg). This claim redirected the whole program for a day.

**B2. "Every nest above 64k cells is compromised."** Stated from a whole-domain max-norm.
**FALSE:** localizing by cell class shows the nest interior is machine-precision exact (7.0e−16)
and only the outermost face fails. Nested results were not invalidated.

**B3. "ERA5's coarse vertical sampling flattens the hodograph and kills SRH."** Falsifiable
prediction made, then tested: re-downloaded at 27 levels (5 below 1 km instead of 2).
**SRH 0–1 km went 70 → 80.** Prediction failed.

**B4. "Horizontal averaging over the domain washes out SRH."** **FALSE:** nearest single column
82 vs full-box average 75.

**B5. "The real-case storm decayed because it drifted into the Davies zone."** **FALSE:** on a
108 km domain the storm stayed 34.8–51.6 km from any boundary and decayed on an identical curve
(w_max 7.5 at t≈975 s in both domains). I had read a dead storm's drift position as its cause of
death.

**B6. "The divergent nest boundary drives the spurious edge vorticity."** **FALSIFIED:**

| nest | solver | face0 residual | ζ edge/interior |
|---|---|---|---|
| 36,000 cells | direct | 1.28e−04 | 323 |
| 230,400 | lowmem | 2.85e−01 | 778 |
| 518,400 | lowmem | 2.17e−01 | 217 |

No relationship. The edge vorticity now looks **intrinsic to Davies-zone nesting** (shear from
nudging toward an interpolated parent), i.e. expected behaviour that `border_frac` exists to
exclude, not a bug.

**B7. Headline intensity numbers.** Peak surface V_rot 12.28 m/s → retracted (moving-nest tracker
steered by sponge vorticity; same configuration gives 8.47 once fixed). With A7 also corrected,
a 600 m run that read V_sfc 12.95–14.45 on the broken radius reads **4.9–6.7** on the correct one.
**Treat every large V_rot in this study's history as suspect until re-measured.**

**What actually explained the low SRH (B3/B4):** the environment was sampled in the wrong PLACE.
SRH 0–1 km across the 143 ERA5 columns runs −29 / 70 (median) / **227** (max) with a monotone
E/SE gradient, while low-level moisture is nearly uniform (15.0 vs 15.6 g/kg) — a wind-shear
gradient. The storm's own column reads 82 because reanalysis smears the storm's circulation into
it. Sampling the inflow sector (~110 km E, 35 km S) gives SRH 0–1 = **158**, matching the
documented ERA5 value for this case (~146–152).

---

## C. The recurring failure pattern (the thing I would most like reviewed)

**Anything expressed in CELLS silently rescales with resolution, and every occurrence corrupted a
resolution comparison.** Three independent instances:

1. `relax_width` in cells → sponge shrinks 3× per refinement level (A2).
2. `R = max(3, int(radius_m/dx))` → different physical measurement disk per mesh (A7).
3. Interior margins / `border_frac` conventions varying across levels (earlier: per-level
   `radius_m` scaling made V_rot incomparable across the cascade).

Every one of these was invisible in the code review and only appeared when a physical quantity
was compared across two meshes. I suspect there are more.

**Second pattern:** diagnostics keyed on `p_dyn` are unreliable on nests (A8), and I built and
briefly trusted one before catching it.

---

## D. Physics results that appear to survive

These are internal idealized-vs-idealized comparisons and do not depend on the defects above:

- **Tilting geometry is the bottleneck, not source vorticity.** Cold-pool baroclinic horizontal
  vorticity generation measured 1.47e−4 s⁻², **22× the tilting rate** (6.7e−6). Only ~4.5 % is
  tilted to vertical. Alignment cos θ ≈ +0.04…+0.09, and **−0.146 at the surface** (anti-aligned).
- **Free evolution builds the right geometry.** In a freely-evolving (not forced) supercell,
  low-level tilting alignment rose 10× (0.02 → peak 0.20) and streamwise fraction rose monotonically
  0.40 → 0.64, with |ζ| peaking afterwards — tilting first, then stretching.
- **Surface-layer closure FORM is the switch, not the drag coefficient.** At dz₁ = 6.2 m,
  surface/aloft ratio: bulk 0.53 | log-law-C_d in bulk form 0.51 | **stress-divergence 0.69** |
  **stress-divergence + log-law 0.82 (connected)** | drag-off 0.88. The bulk sink rate
  `C_d|V|/dz₁` diverges as the mesh refines; the physical form `τ(z)=τ_s(1−z/h)` is mesh-independent.
- **Near-surface vertical resolution dominates horizontal.** Surface/aloft ratio 0.08 at
  dz₁ = 49.6 m vs 0.53 at dz₁ = 6.2 m (**6.6×**), against +22 % for a 9× horizontal refinement.

**Ruled out by measurement as the limiter:** mesh resolution (partially — now re-opened),
environmental SRH, updraft strength, cold-pool source strength.

---

## E. Open problems / what I would like a second opinion on

1. **A6** — what is the correct projection treatment for a nest boundary that legitimately carries
   net mass flux? Removing only the discrete imbalance, or something else?
2. **Is the nest edge vorticity really intrinsic?** (B6 suggests yes.) If so, what is the correct
   minimum trusted-interior margin, and is one-way Davies nesting the wrong tool for resolving a
   vortex whose inflow crosses the nest boundary?
3. **The resolution question is still unanswered.** A 125 m tornado RMW needs dx ≲ 22 m
   (~6 cells across the core). Every nested attempt at 22 m produced a contaminated measurement.
   Is there a formulation that avoids Davies nesting entirely for this?
4. **The grid is coarse at ground level.** The parent has **ONE cell below 100 m** (dz₁ 79.8 m).
   Measured cost of fixing it at nx=120: dz₁ 5.5 m → dt 0.14 s → 195 min maturation (affordable);
   dz₁ 1.4 m → 65 h (not). Is there a better approach than brute-force vertical CFL (implicit
   vertical diffusion/advection, or a surface-layer parameterization valid at dz₁ ~50 m)?
5. **Validation target.** The 26 m/s is a beam-limited radar observable (~250 m resolution volume,
   ~460 m AGL), not the tornado's wind (Moore 2013 was EF5, ≥89 m/s). Comparisons should go through
   a radar forward operator (exists: `atmospheric_data/radial.py`). Is that the right framing?
6. **Case mixing.** The 26 m/s target is real Moore 2013 but the productive runs use an analytic
   sounding (CAPE 2226, shear 42, SRH 0–3 = 648) vs the real environment (CAPE ~2300–3100,
   shear 22–28, SRH 0–3 ~213–243). All percentages against 26 m/s are therefore invalid.
7. **Would an inverse-problem formulation help?** Considered: parameter estimation (surface-layer
   closure coefficients) via derivative-free ensemble methods, since no adjoint exists. Concern:
   inverse methods absorb structural model error into parameters, and this model is still yielding
   new defects.

---

## F. Environment / reproduction

- GPU RTX 4050 6 GB, CuPy 14.2. Parent 120×120×48 (dx 600 m) matures to t=2800 s in ~9 min;
  240×240×48 (dx 300 m) in ~35 min; 2.76M cells uses ~3.0/6.1 GB.
- Test suite: 370 passed, 3 skipped.
- Everything above is committed with the measurements in the commit messages
  (`git log` around 2026-09-03/04).

---

## G. ADDENDUM (2026-09-04) — the validation target itself is an artifact

Building the radar forward operator forced an inspection of the cached KTLX extraction
(`outputs/nexrad_moore/ktlx_velocity.npz`, 0.5211° sweep, 2013-05-20 20:20:58Z, 65,907 gates).
Two independent defects in the observational target:

**G1. `V_rot = 26 m/s` is the NYQUIST VELOCITY, not a measurement.**

| region about Moore | gates | min v_r | max v_r |
|---|---|---|---|
| 1 km | 74 | −9.0 | +18.5 |
| 2 km | 290 | −25.5 | +26.0 |
| 3 km | 652 | **−26.000** | **+26.000** |
| 5 km | 1,840 | **−26.000** | **+26.000** |
| 8 km | 4,788 | **−26.000** | **+26.000** |
| 15 km | 17,584 | **−26.000** | **+26.000** |

Every region ≥3 km saturates at exactly ±26.000. Adjacent-gate velocity jumps confirm folding:
**9 neighbouring pairs differ by >40 m/s, 7 of them in the 45–52 m/s band** — i.e. ~2×Nyquist =
52 m/s, the classic aliasing discontinuity. So `(max − min)/2` over any region containing a
folded gate returns exactly 26 by construction.

**The true rotational velocity is a LOWER BOUND of 26 m/s and cannot be recovered from this
field without dealiasing.** Every percentage this study has quoted against "26 m/s" was
measured against an instrument ceiling.

**G2. The recorded couplet was extracted 21 km from Moore.**
The saved couplet (`V_rot 26.0`, separation 253 m) sits at 4.9 km range from KTLX. Moore is at
19.3 km range — **21.2 km away** from the extracted feature. The extraction script takes the
strongest inbound/outbound pair over the whole cropped field, and that pair is near-radar
clutter, not the Moore mesocyclone.

At Moore's actual location the geometry is:

| quantity | recorded in docs | **measured at Moore** |
|---|---|---|
| range from KTLX | (implied ~30–41 km) | **20.1 km** |
| beam diameter | 250 m assumed | **325 m** |
| sampling height | ~460 m AGL | **207 m AGL** |
| couplet separation | 253 m | **1015 m** |
| V_rot | 26 m/s | **≥26 m/s (aliased)** |

**Required before any quantitative validation:** dealias the Level II velocity field (Py-ART
`dealias_region_based` or similar), re-extract the couplet constrained to the Moore mesocyclone,
and re-derive the target with its true scan geometry. Until then no model/observation ratio is
meaningful.
