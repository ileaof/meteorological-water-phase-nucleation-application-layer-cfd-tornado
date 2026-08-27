# Handoff — storm_dynamics (rotating supercell / tornadogenesis core)

**Status:** M1 (rotating supercell) delivered & verified; M2 (low-level rotation)
delivered & verified; **M3 phases 1 & 2 (one-way nesting) delivered** — phase 1
(static/frozen-parent) intensifies near-surface ζ ~2.4× over a short window;
phase 2 (concurrent/time-evolving parent boundary) sustains the nest as long as
the parent drives it. Remaining: storm-following moving nest, much higher
refinement, two-way / adaptive (AMR) nesting.

## Objective

Add the capability to **simulate idealised tornadogenesis**: a supercell that
develops a mid-level mesocyclone under shear and, under a curved hodograph with
surface friction and a cold pool, **near-surface rotation** — the physical proxy
for a tornado. Idealised only: no data assimilation, no real event, no
observational verification. The deliverable is a *rotational dynamical core*
validated against the classical supercell results, reusing the repo's already
validated physics.

## Decisions locked (user-approved)

1. **New package** `src/storm_dynamics/` (name confirmed by the user), a **fork**
   of the `meteorological_flow` dynamical core — not a destructive rewrite. The
   `_engine/` nucleation kernel and the scientific behaviour of
   `meteorological_flow` are untouched.
2. **Reused as a physics library, unchanged:** `grid`, `pressure_solver`
   (anelastic projection, already periodic), `advection.advect_center_massflux`,
   `buoyancy`, `diffusion`, `microphysics_coupling` + `precip_microphysics`,
   `base_state.weisman_klemp` (thermodynamics), `diagnostics`, `thermodynamics`,
   and the validated nucleation kernel.
3. **Rewritten (the rotational core only):** conservative flux-form staggered
   momentum advection, f-plane Coriolis, LES (Smagorinsky) closure replacing the
   Rayleigh drag + velocity clip, surface bulk-drag, curved hodograph, rotation
   diagnostics.
4. **`--validate` integrity:** a pre-existing line-ending artifact of the
   `_tornado` checkout had left the SHA-256-guarded ice reference script with LF
   endings (hash `755d7daf…`) while the guard expects the CRLF form (`c9fa9c01…`).
   Restored the canonical CRLF bytes on that one file (content unchanged — only
   EOLs) and added `.gitattributes` (`-text`) so it cannot re-normalise.
   `--validate` and `test_18_reference_preserved` are green again. **No science
   was changed.**

## File map (new)

```
src/storm_dynamics/
  __init__.py       package doc + scope statement
  config.py         StormConfig(sim: SimulationConfig, dyn: StormDynamicsConfig);
                    build_storm_config(...), storm_config_from_yaml(...), coriolis_f
  soundings.py      curved quarter-circle / unidirectional hodograph (item 5);
                    SRH (Davies-Jones), bulk shear, Bunkers storm motion
  momentum.py       conservative flux-form staggered momentum advection (item 1)
  coriolis.py       f-plane Coriolis on the perturbation wind (item 2)
  turbulence.py     Smagorinsky LES: strain -> K_m, momentum + scalar SGS (item 3)
  surface_drag.py   bulk aerodynamic drag on the lowest level (item 4)
  rotation.py       zeta, 3D vorticity, updraft helicity, trackers, report (item 7)
  core.py           StormSimulation: the anelastic time loop wiring 1-7 (item 6)
  plotting.py       rotation figures: mid-level/near-surface zeta slices (the split
                    couplet), w, hodograph+SRH, and the rotation time series
  handoff.md        this file
tests/test_storm_dynamics.py      13 fast unit tests (momentum/coriolis/LES/drag/SRH + short run)
tests/test_storm_milestones.py    M1 (+ M2) regression
examples/supercell_tornadogenesis.py   runnable demo (both scenarios), prints diagnostics
configs/storm_supercell.yaml           M1 declarative config
configs/storm_tornadogenesis.yaml      M2 declarative config
docs/storm_dynamics_guide.md           model, can/cannot claims, limits, references
pyproject.toml                         packages.find include += "storm_dynamics*"
```

## Build order (prompt items 1→7) — each verified before the next

1. **Momentum advection** — verified: uniform flow → exactly zero tendency;
   domain-integrated momentum conserved to ~1e-17 (machine); finite everywhere.
2. **Coriolis** — verified: `f(36°N)≈8.6e-5`; rotates the perturbation only; a
   pure base-state wind feels no force.
3. **LES** — verified: `K_m ≥ nu_background` at rest, grows with resolved strain.
4. **Surface drag** — verified: retards the lowest level only; disabled = no-op.
5. **Curved hodograph** — verified: right-mover **positive** SRH; curved ≫
   straight (SRH₀₋₃ ≈ 470 vs 130 m²/s²); shear₀₋₆ > 20 m/s; `v0(z)` populated.
6. **Microphysics / cold pool** — reused `MicrophysicsCoupler`; conservation holds.
7. **Rotation diagnostics** — ζ, updraft helicity, mid-level / near-surface
   trackers, SRH/shear from the hodograph.

## Verification (demonstration results)

**M1 (32×32×40, Δx≈1.25 km, 40 min):** w_max ≈ 20 m/s; ζ_max/ζ_min ≈ +0.050 /
−0.050 s⁻¹ → **storm splitting**; mid-level mesocyclone ≈ 1.1×10⁻² s⁻¹; updraft
helicity ≈ 230 m²/s² (peak ≈ 340). Conservation: mass |div| ≈ 1.6e-4, energy
rel-err ≈ 2e-3, water rel-err ≈ −1.5% (top damping-layer boundary flux — the same
channel the parent solver documents). The M1 regression test uses a faster
24×24×36 / 20 min run (~30 s) that still splits and forms a mesocyclone.

**M2 (curved hodograph + drag + cold pool, near-surface levels clustered,
dz₀≈100 m, 40 min):** near-surface ζ spins up and is **sustained** (~3.3×10⁻³
s⁻¹) while a straight-hodograph / no-drag control peaks lower (~2.5×10⁻³) and
**decays** (~1.4×10⁻³) — the curved+friction case holds ~2.4× the low-level
rotation. Updraft physical (w_max ≈ 18 m/s), conservation unchanged.
**Stability:** a strong low-level-SRH environment on a coarse uniform grid
over-amplifies the updraft (w blow-up to ~110 m/s); the shipped M2 config fixes
this with `z_stretch=1.05`, `U_max=18`, `z_turn=2000`, `C_s=0.22` — documented
demonstration choices, not tuning of the physics.

**Existing suite:** the baseline was 108 tests (one, `test_18_reference_preserved`,
had been red purely from the EOL artifact above; now green). storm_dynamics adds
13 unit tests + the milestone tests. `_engine/` untouched → `--validate` green.

## Key implementation notes / gotchas

- **Flux form needs a divergence-free transporting velocity.** Momentum is
  advected with the *projected* velocity, so the flux form equals the advective
  form (`u_j ∂u_i/∂x_j`) and conserves the integral. Advecting with the divergent
  predictor would not.
- **Environment persistence.** LES momentum diffusion and Coriolis act on the
  *perturbation* `(u−u0, v−v0)`, so the environmental hodograph is not mixed or
  inertially oscillated away; the strain/`K_m` still use the full field. Momentum
  *advection* uses the full field (vertical advection of environmental momentum is
  physical). Surface drag acts on the full lowest-level wind (friction is real).
- **No velocity clip.** The only bound is `dyn.v_guard` (150 m/s) — an extreme
  numerical guard, documented, that should never bite in a resolved run.
- **Anelastic residual.** `mass_continuity_residual_norm` ~1e-5 (projection
  enforces `div(ρ₀u)=0`); the *plain* `div(u)` is intentionally nonzero for the
  anelastic system (`= −w ∂ρ₀/∂z / ρ₀`). Use the anelastic (ρ₀-weighted) residual.
- **SRH sign convention:** `SRH = ∫[(v−cy)∂u/∂z − (u−cx)∂v/∂z]dz` (Davies-Jones
  with the leading minus folded in) → veering (clockwise) hodograph = positive SRH
  = right-mover. The quarter-circle turns to `(U_max, −U_max)` for this sign.

## Remaining work

- **M2 quantification / tuning:** confirm near-surface ζ magnitude and its
  location on the forward-flank cold-pool interface; may want a finer near-surface
  Δz (`z_stretch`) and a stronger/earlier cold pool.
- **M3 phases 1 & 2 done** (`nesting.py`, `examples/tornado_nest.py` [`--concurrent`],
  `tests/test_storm_nesting.py`): one-way nest — exact parent→nest trilinear
  interpolation, finer nest grid, stable/conserving nested integration with
  Davies-style border relaxation. **Phase 1** (frozen parent): near-surface ζ
  intensifies ~2.4× over a 120 s window at 3× refinement (Δx 1.3 km→0.44 km);
  valid only ~2–3 min before the frozen border decays the storm (ζ max then drifts
  to the edge — use `interior_near_surface_zeta` to read the physical interior
  vortex). **Phase 2** (`run_concurrent_nest`): the parent steps alongside the nest
  and feeds time-evolving boundaries, sustaining the nest as long as the parent
  drives it; a small SGS boost (`les_boost`) + tighter CFL keep the sharpening
  vortex stable. **Remaining M3:** a storm-following MOVING nest (a fixed nest
  loses features that advect out), much higher refinement toward O(10–100 m), and
  two-way / adaptive (AMR) nesting. Gotcha: the nest deep-copies `parent.dyn`
  (shared reference would let nest tweaks mutate the parent).
- **Figures done:** `plotting.py` + the example's `--plots` flag write the
  rotation slices (the split couplet), hodograph and time series. **Still
  optional:** NetCDF field output (the core returns a report dict + `history` +
  the `tracker`; it does not write `flow.nc` — reuse `meteorological_flow.io` if
  persisted 3-D fields are wanted).
- **TKE-1.5 (Deardorff) closure:** advertised as `les_model='tke15'` but not
  implemented — `strain_and_viscosity` raises `NotImplementedError` for it (only
  `smagorinsky` and `none` are implemented). A prognostic TKE equation is the
  remaining work.

## Done since the first handoff

- **Rotation figures** (`plotting.py` + example `--plots`): mid-level / near-surface
  ζ slices (the split couplet), w, hodograph+SRH, rotation time series. Plus an
  **animated GIF** of the evolving vorticity (`--animate`, Pillow — no ffmpeg;
  frame capture via `run(capture_frames=True)`) and a **history CSV** (`--csv`).
- **Kernel-nucleation coupling** (`build_storm_config(couple_nucleation=True)` /
  example `--kernel-nucleation`): the validated 2nd-order nucleation rate J is
  evaluated on the post-transport state each step and fed to the microphysics as
  the embryo source (eq39 pathway), exactly as `meteorological_flow` does. Off by
  default (the lookup-table build is slow). Verified: conserves water, stays
  finite. `report["kernel_coupled"]` records it.
- **GPU support** (`build_storm_config(device=...)` / example `--device`): the
  core is backend-agnostic (every hot-loop module uses `grid.xp`; zero direct
  `numpy` in `momentum`/`turbulence`/`coriolis`/`surface_drag`/`rotation`), so it
  runs on CPU (NumPy) or GPU (CuPy). Default `cpu` (preserves prior behaviour);
  `gpu` fails loudly if CuPy/CUDA absent; `auto` falls back. Verified GPU parity
  with CPU (bit-identical diagnostics) via `test_storm_gpu_matches_cpu` (skips
  without a GPU). `report["backend"]` records it. Note: at demo grid sizes the
  direct pressure solve stays on the host, so CPU is usually faster; GPU wins on
  large grids (CG on GPU).
