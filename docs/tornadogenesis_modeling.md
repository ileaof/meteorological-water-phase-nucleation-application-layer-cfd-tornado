# Physically-consistent supercell & tornadogenesis modelling

*How the model represents the storm chain, the diagnostic instruments, how to configure and run it,
and the limitations. Companion to [`TORNADOGENESIS_FINDINGS.md`](TORNADOGENESIS_FINDINGS.md)
(the honest quantitative results) and [`ROADMAP.md`](ROADMAP.md).*

The tornado is **never imposed** — no Rankine/Lamb–Oseen vortex, no funnel. Rotation emerges from the
governing equations: convective instability → deep updraft → tilting of horizontal (shear +
baroclinic) vorticity into the vertical → rotating updraft → mesocyclone → downdrafts + cold-pool
baroclinity → low-level convergence & stretching → surface-connected tornado-like vortex. A test
asserts no `impose_vortex` code path exists.

## 1. Formulation (unchanged core; audited)

| Aspect | Choice | Where |
|---|---|---|
| Equations | **anelastic** ∇·(ρ₀u)=0 (Boussinesq fallback; no compressible) | `core.py`, `pressure_*.py` |
| Prognostic | u,v,w, θ, qv, ql, qi, qr, qs, qg, qh, (TKE for tke15) | `state.py`, `core.py` |
| Grid | staggered Arakawa-C, geometric vertical stretch, periodic lateral / walled z | `grid.py` |
| Time | explicit, adaptive CFL (advective + diffusive + **sedimentation**) | `core.py:_dt` |
| Advection | conservative flux-form MUSCL/minmod (order 2) | `momentum.py`, `advection.py` |
| Pressure | exact anelastic projection; FFT+tridiag / Jacobi-CG low-memory for large grids | `core._project`, `pressure_fft.py`, `pressure_iterative.py` |
| Turbulence | Smagorinsky (Lilly) or TKE-1.5 Deardorff | `turbulence.py` |
| Microphysics | ice scheme: cond/evap/freeze/melt/dep/subl/autoconv/accr/sediment + latent heat; validated nucleation kernel (shifted equilibrium) | `precip_microphysics/*`, `nucleation_adapter.py` |
| Surface | bulk momentum drag + **sensible/latent heat fluxes** (opt-in) | `surface_drag.py`, `surface_fluxes.py` |
| Nesting | AMR concurrent multi-level, storm-following (Galilean), **two-way coupling** (opt-in), conservative restrict/prolong | `nesting.py`, `amr.py` |
| CPU/GPU | NumPy / CuPy via `grid.xp`, automatic CPU fallback | `backend.py` |

## 2. Diagnostic instruments (new)

- **Vorticity budget** (`vorticity_budget.py`): each term of Dζ/Dt as a field — advection,
  stretching ζ∂w/∂z, tilting ξ∂w/∂x+η∂w/∂y, **baroclinic** (1/ρ²)(∇ρ×∇p)_z, divergence,
  diffusion, surface friction — plus streamwise/crosswise decomposition, `budget_layer_summary`
  and `dominant_mechanism` (names baroclinic vs tilting vs stretching in a layer). Multi-layer UH
  and bulk Richardson number in `rotation.rotation_report` / `soundings`.
- **Vortex diagnostics** (`vortex_diagnostics.py`): circulation Γ=∮u·dl=∬ζdA, tangential velocity
  v_θ=−u sinφ+v cosφ about a detected centre (ζ-max snapped to p′-min), core radius, pressure
  deficit, swirl proxy.
- **Cold pool** (`coldpool.py`): θv′, buoyancy, intensity C=√(2∫−b dz), footprint, gust-front
  gradient + convergence, downdrafts (RFD not forced cold — reported).
- **Classification** (`classification.py`): resolution-aware ladder NO_DEEP_CONVECTION →
  ORDINARY_CONVECTION → SUPERCELL → LOW_LEVEL_MESOCYCLONE → TORNADO_LIKE_VORTEX →
  SURFACE_CONNECTED_TORNADO_LIKE_VORTEX, keyed on v_θ / Γ / Δp / UH (not raw ζ); persistence gates
  the tornado-like tiers.
- **Two-scale temperature gradient** (`micro_gradient.py`): the resolved **macro** |∇T| vs the
  sub-grid **micro** |∇T| via a documented closure |∇T|_micro = |∇T|_macro·(Δ/r_eff)^{2/3}
  (Obukhov–Corrsin inertial-convective scaling; r_eff=max(r_particle, Batchelor scale(ε))). Drives
  only the nucleation term — never the resolved energy equation. Named diagnostics:
  macro/micro gradient, local supersaturation, nucleation-rate + latent-heat proxies.
- **Radar operators** (existing): reflectivity dBZ + radial velocity V_r=V·r̂ (compare to NEXRAD).

## 3. Configure & run

Unified overlay (additive; both legacy schemas still work):

```python
from storm_dynamics.tornado_config import load_tornadogenesis_config, run_diagnostics
from storm_dynamics.core import StormSimulation
cfg = load_tornadogenesis_config("config/tornadogenesis.yaml")   # impose_vortex=false enforced
sim = StormSimulation(cfg.storm)
# ... integrate ...
diag = run_diagnostics(sim, cfg)     # rotation, vorticity budget, vortex, cold pool, class, macro/micro
```

Two-way coupling in the deep cascade (the Attempt-G lever): `run_multilevel_nest(..., two_way=True,
storm_motion=(cx,cy))`. Surface heat fluxes: `SurfaceFluxConfig(enabled=True, ...)`. cos² trigger:
`storm_dynamics.initiation.smooth_bubble` / `multi_bubble`.

## 4. Reproduce

```bash
python -m pytest tests/test_vorticity_budget.py tests/test_vortex_diagnostics.py \
  tests/test_coldpool.py tests/test_classification.py tests/test_micro_gradient.py \
  tests/test_surface_fluxes.py tests/test_two_way_coupling.py tests/test_initiation.py \
  tests/test_tornado_config.py tests/test_experiment_matrix.py -q      # all new instruments

python examples/experiment_matrix.py                 # control matrix (shear->rotation, ...)
python examples/tornado_diagnostic_figures.py        # -> docs/media/storm/tornado_diagnostics.png
# real-data supercell + two-way cascade: scratchpad/moore_twoway_ab.py (Attempt G)
```

## 5. Limitations (honest)

- At Δx ≥ ~30–50 m we resolve a tornado-**scale** circulation, not the true core — hence
  "TORNADO_LIKE_VORTEX"; the first cell centre height is reported and exact ground contact is not
  claimed when the mesh cannot resolve it.
- The macro→micro gradient is a **sub-grid closure**, not a resolved field; it feeds only the
  nucleation term.
- No compressible core, no terrain-following coordinates, β-plane Coriolis (f-plane only) — by
  design; lateral BCs are periodic (Davies relaxation exists but is not wired into the main loop).
- The quantitative Moore-2013 story (why the low-level V_rot plateaus, and how resolution / SRH /
  updraft / **two-way coupling** move it) is in `TORNADOGENESIS_FINDINGS.md`: two-way coupling is the
  lever that breaks the ~6 m/s one-way ceiling (8.1→11.6 at the nest, 6.1→15.0 at the parent).

## 6. Continuity

Development state and next steps: this file, `TORNADOGENESIS_FINDINGS.md`, `ROADMAP.md`,
`src/storm_dynamics/handoff.md`. All new diagnostics are pure/`xp`-generic with CPU==GPU parity
tests; new physics (two-way coupling, surface fluxes) is opt-in so default runs are byte-unchanged.
