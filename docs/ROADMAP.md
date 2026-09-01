# Continuity & generalisation roadmap — `storm_dynamics`

**Purpose.** This is the **start-here document after a session loss.** It snapshots what
is built and verified, then lays out the remaining generalisation of this idealised
rotating-convection core toward a fuller atmospheric model — each step with concrete file
/ function pointers, a first action, and a "definition of done". It is a *plan*, honestly
scoped; nothing below is implemented unless the snapshot says so.

Companion docs: `src/storm_dynamics/handoff.md` (build log + gotchas),
`docs/amr_design.md` (the rigorous AMR-projection design + call sites),
`docs/storm_dynamics_guide.md` (what the model can/can't claim). Persistent memory index:
`…/memory/MEMORY.md` → `amr-composite-projection-state.md`.

---

## 0. Snapshot — what is DONE and VERIFIED (as of 2026-08-31)

**Rotational core (M1/M2), reused physics.** Conservative flux-form momentum advection,
f-plane Coriolis, Smagorinsky LES, surface drag, curved hodograph + SRH, rotation
diagnostics. Reuses grid, anelastic projection, transport, microphysics, nucleation kernel
unchanged. CPU/GPU (NumPy/CuPy), GPU parity tested. `--validate` green (engine untouched).

**Nesting (M3 phases 1/2/2b/3a).** One-way static, concurrent time-evolving boundary,
storm-following (storm-relative frame), approximate two-way injection. Strongest verified
case: w_max 7.1→13.6 m/s, near-surface ζ interior +3.95× (`examples/tornado_nest.py
--follow`, tunable `--u-max/--les-boost/--cfl`).

**AMR algorithms — all hard pieces built + verified** (`storm_dynamics/`):
- `amr.py` — Berger–Colella **refluxing** (mass drift 2e-16 vs 1.1e-4), conservative
  restriction/prolongation, free-stream (0.0).
- `poisson_mg.py` — geometric-multigrid Poisson (h-independent, 2nd order).
- `composite_poisson.py` — the **composite coarse-fine interface stencil** in 1-D/2-D/3-D
  (2nd order incl. corners, conservative) AND the **two-level MAC projection**
  (`div(ρ₀u)=0` across the interface, ~1e-13), now **anisotropic** (Δx≠Δy) + **stretched-z**.
  The unified operator is `solve_composite_hz`; the storm-array projection is
  `composite_project_massflux_hz`.
- `nesting.composite_project_two_level(parent, nest, spec)` — the **call site**: forms
  `ρ₀u*` on both levels, one composite solve, recovers `u=m/ρ₀`, writes back. Verified on
  real `FlowState`s (`test_composite_project_two_level_call_site_divergence_free`).
- `amr_port.py` — AMReX↔our-physics binding (drift 0.0 on WSL/pyAMReX).

**Key fact.** The composite projection **is** the anelastic projection in mass-flux
variables — no new *algorithm* is needed to make the storm's AMR pressure solve
consistent across the interface. Nothing remains on the projection itself; regridding, resolution, and
(much later) the leap to real forecasting.

**Composite projection IN the time loop (DONE 2026-08-31).** `run_concurrent_nest(
composite_projection=True)` (`nesting.py`) skips both per-level `project_anelastic`
calls and runs `composite_project_two_level` once per nest sub-step: `core._step` and
`NestedStormSimulation._step` are split into `_predictor → _project → _transport` so the
projection can be deferred; with `follow` the nest is reconciled to the ground frame for
the joint solve and shifted back after; the nest velocity sponge is dropped (the interface
coupling replaces it); `follow` footprints are snapped cell-aligned/matched-z
(`_snap_aligned_spec`). Verified `|div(ρ₀u)|` at the interface ~4e-18 every sub-step on a
stepping run, stable, w_max ≥ the sponge path, mass residual ~1e-15. `--composite` in
`examples/tornado_nest.py`.

**Tests.** `tests/test_storm_nesting.py` (composite 1-D/2-D/3-D, projection 2-D/3-D, wall
BC, mass-flux bridge, stretched metric, unified operator, full-storm projection,
anisotropic, call site, composite-in-the-time-loop — plus refluxing/restriction).
`tests/test_storm_dynamics.py`, `tests/test_storm_milestones.py`.
Run: `python -m pytest tests/test_storm_nesting.py -q` (expect 25 pass, 1 skip).

---

## 1. ~~NEXT (immediate) — composite projection in the time loop~~ ✅ DONE (2026-08-31)

Landed exactly as planned: `composite_projection=False → True` on `run_concurrent_nest`;
after both levels' predictors it skips the per-level projections and calls
`composite_project_two_level` once per **nest sub-step** (so the interface stays consistent
while the nest sub-cycles); the nest velocity sponge is dropped; `follow` footprints are
snapped cell-aligned + matched-z; frames are reconciled to the ground frame for the joint
solve. See §0 "Composite projection IN the time loop". The test
`test_composite_projection_in_time_loop` (`tests/test_storm_nesting.py`)
asserts the interface divergence ~machine precision over the window and w_max ≥ sponge
path.

**NEXT: §2b remaining** — momentum (velocity) refluxing between levels (`amr.TwoLevelReflux`
per pair; the multi-level up-feedback is scalar-only so far), a 3rd level (Δx≈50 m = true
funnel scale), and combining the multi-level driver with `follow`/`regrid` for a long-lived
storm. ★  (§2a regridding + §2b increments 1 & 2 — recursive nesting + the concurrent
multi-level driver `run_multilevel_nest` — landed 2026-09-01, verified; demonstrated at
Δx=149 m with ζ_max ~5× the single-level value.)

**Known issue before scale-up (2026-08-31).** Composite mode is *verified* for
correctness (interface `|div(ρ₀u)|` ~machine precision, mass residual ~1e-15) but
*not yet stable at long windows on the render-scale grid* (parent 24×24×40, nest
24×24×40 @ r=3): the velocity sponge is dropped in composite mode (the interface
coupling replaces it) and, with no lateral damping channel left, the nest
over-intensifies — w grows smoothly (not a noise spike, max ~5 cells from the
border, mid-level) to 100 m/s @ window 300 (u_max 18), 75 @ 210 (les_boost 1.5,
cfl 0.15), 44 @ 180 (u_max 14). The small-grid in-loop test (14×14×16, 18 s) is
fine. *Next stability item:* a border treatment compatible with the composite
coupling (absorbing/hyperviscous band that does NOT nudge interface faces back
toward the coarse state), then re-validate the w_max trajectory before using
composite mode for presentation figures. Until then the sponge-path figures
remain the demonstrated-physical ones.

---

## 2. AMR completion (the path to a resolved funnel)

**2a. Adaptive / dynamic regridding (Berger–Oliger).** Detect where to refine each N steps
(criterion on |ζ|, |∇w|, or updraft-helicity), cluster tagged cells into a new nest
footprint, re-create/move the nest during integration, restrict/prolong state across the
regrid.

- **Increment 1 — detection + clustering primitives (DONE 2026-09-01).**
  `nesting.tag_cells(state, grid, field='uh'|'zeta', frac)` → boolean (nx,ny) tag map;
  `nesting.cluster_to_box(tags, margin)` → `(i0,j0,ncx,ncy)` bounding box (single-cluster
  Berger–Rigoutsos surrogate); `nesting.regrid_spec(parent, refine, field, frac, margin)`
  → an **aligned** `NestSpec` (matched-z) covering the tagged vortex — a *data-driven*
  footprint that replaces the constant-C follow. Verified: a synthetic solid-body-rotation
  blob is tagged and the returned aligned footprint contains it
  (`test_adaptive_regridding_primitives_track_the_vortex`).
- **Increment 2 — time-loop re-centring (DONE 2026-09-01).** `regrid_nest(old_nest,
  parent, new_spec)` re-creates the nest at a shifted **aligned** footprint, preserving the
  old nest's fine field in the overlap by an **exact integer fine-cell shift** (`_shift_copy`,
  no interpolation/smoothing) and filling the newly-exposed strip from the parent; carries
  `t/step/tracker/history`. `run_concurrent_nest(regrid_interval=N, regrid_field, regrid_frac)`
  (ground frame, `follow=False`): every N parent steps it re-centres the fixed-size nest on
  the parent's tagged vortex — **data-driven storm-following** replacing the constant-C
  follow. Example flag `--regrid-interval`. Verified: the overlap is preserved exactly
  (`test_regrid_nest_preserves_fine_structure_in_overlap`) and a nest starting off the
  vortex hops to re-centre on it (`test_regrid_interval_recentres_nest_on_the_vortex`).
  *Remaining within 2a (later):* variable footprint SIZE (not just re-centring; regrid_nest
  degrades gracefully to a parent re-interpolation if the size changes), and a
  `conservative_restrict`-checked conservation report across each regrid.

*Done (2a):* the nest tracks the vortex data-driven over a run, preserving fine structure
across each re-centre.

**2b. Higher refinement + subcycling toward O(10–100 m).** Multiple nest levels (r=3–4 per
level), each subcycling in time (finer dt), refluxing between levels (`amr.TwoLevelReflux`
generalised to N levels). This is what actually lets a low-level funnel condense to the
surface — the physics gap behind the "no funnel-to-ground" caveat.

- **Increment 1 — recursive nesting (DONE 2026-09-01).** `NestedStormSimulation` **composes**:
  its `parent` may itself be a `NestedStormSimulation`, so a nest-of-a-nest gives
  `Δx = parent.dx / r²` (1.3 km → ~440 m → ~150 m; a 3rd level → ~50 m). Verified it builds,
  matches z through the stack, and runs stably
  (`test_recursive_nesting_second_level_refines_further`). Wired into the renderer:
  `examples/render_tornado_3d.py --levels 2` stacks a second level over the finer updraft
  and renders it (`dx = parent.dx/9`).
- **Increment 2 — concurrent multi-level driver (DONE 2026-09-01).** `run_multilevel_nest(
  parent, specs, window)` builds the chain (each level's parent is the level above) and, per
  parent step, drives every finer level recursively: each **sub-cycles in time under the
  level above** (its own finer dt), with time-evolving boundaries blended down and a
  **conservative scalar restriction up** onto the parent overlap (`restrict_up`). `specs`
  entries may be callables `(grid_above)->NestSpec`. Verified: a 2-level stack sub-cycles the
  finest level (Δx=parent.dx/r²) over the window, stable/finite
  (`test_multilevel_concurrent_driver_sustains_finest_level`). Demonstrated at Δx=149 m
  (`render_tornado_3d.py --levels 2`): ζ_max 6.7×10⁻² s⁻¹ — ~5× the single-level value; the
  near-surface rotation is now a resolved dense mass, not a sparse gap
  (`docs/media/storm/nest_3d_L2_column.png`).
  **Momentum feedback up** (`restrict_velocity`, opt-in `run_multilevel_nest(
  restrict_momentum=True)`): conservative face average-down of the fine staggered velocity
  onto the coarse overlap faces (mass-flux-preserving), complementing the scalar
  `conservative_restrict`. Verified (`test_restrict_velocity_face_average_down`).
  *Remaining in 2b:* a true **flux-register interface reflux** of the momentum *flux*
  (`amr.TwoLevelReflux` per pair — `restrict_velocity` restricts the overlap, not the
  interface flux), a **3rd level** (Δx≈50 m, true funnel scale), and combining the
  multi-level driver with `follow`/`regrid` so the stack tracks a long-lived storm.
  **Memory note (investigated 2026-09-01):** the finest level's pressure solve is a
  *direct sparse LU*, so a 3rd level at 48³ OOM'd.  The obvious fix — switch large grids to
  the existing **CG** — does **NOT** work: Jacobi-preconditioned CG diverges on the
  stretched anelastic operator (NaN on the periodic pure-Neumann system, residual ~8 after
  5000 iters on the wall system; the storm always uses stretched grids, so this CG path was
  never actually exercised before). Verified/guarded by
  `test_pressure_cg_does_not_converge_on_stretched_operator`; `_pressure_method` keeps
  stretched grids on `direct`.  The **real low-memory fix** is in §3f. Until then keep the
  finer footprints SMALL (`render_tornado_3d.py --sub-half-frac 0.22 --sub-window-frac 0.4`).

*Done (2b):* a resolved near-ground vortex that connects vertically to the meso (ζ continuous
through the depth; the 3-D render then shows a funnel, not a gap).

**2c. Framework path (optional, for scale).** The pure-NumPy AMR is a reference; a
production build wants AMReX `MLMG` + `FluxRegister` (not in the default pyAMReX —
`docs/amr_design.md`), or Chombo/p4est, with MPI domain decomposition. Port the verified
operators onto framework data structures (`amr_port.py` is the seam).

**Remaining composite generalisations.** Non-cell-aligned nests (interface off integer
coarse faces → shifted interface interpolation in `_interface_flux_2d`); vertical
refinement (z-interface, currently matched-z only).

---

## 3. From idealised core → atmospheric MODEL (the big arc, honestly huge)

This is the leap the phrase "modelo de previsão atmosférica" implies. Each item is a
multi-month effort; listed so the direction is explicit.

**3a. Real initial & boundary conditions.** Ingest real soundings / reanalysis (ERA5,
RAP/HRRR) for the base state and, for a limited-area run, time-dependent lateral BCs from a
driving model (relaxation/Davies zone — the nest machinery already does the LBC pattern).
*Seam:* `base_state`, `soundings.py`, the BC application in `core.py`/`bc`.

**3b. Data assimilation.** 3D/4D-Var or EnKF to fit initial conditions to observations —
the actual difference between a *simulation* and a *forecast*. Standalone module consuming
the model as the forward operator; large.

**3c. Physics completeness.** Implement the advertised **TKE-1.5 (Deardorff)** closure
(`turbulence.py` raises `NotImplementedError` for `tke15`); add radiation, a land-surface /
surface-flux scheme, and a more complete microphysics/coupling path; optionally a
non-hydrostatic **compressible** core alongside the anelastic one for steep terrain / strong
vertical accelerations.

**3d. Geometry & Coriolis.** f-plane → β-plane / full variable f; map projections and a
curved-earth metric for domains large enough that sphericity matters.

**3e. Verification.** Beyond idealised conservation: standard benchmarks (density current /
Straka, Weisman–Klemp supercell intercomparison, warm-bubble) and, once forecasting,
observational verification (against radar/METAR). This is how "indicative" becomes
"quantitative".

**3f. Performance & scaling — the low-memory pressure solver (★ the fine-nest memory fix).**
The direct sparse-LU pressure solve stores a full factorisation and OOMs at ~48³; CG with
the current Jacobi preconditioner diverges on the stretched anelastic operator (see §2b
memory note). The clean fix exploits the storm's structure — **uniform-x,y + stretched-z
only**, so the horizontal Laplacian diagonalises under a transform:
- **Parent (periodic x,y):** FFT in x,y → for each horizontal wavenumber `(kx,ky)` a
  **tridiagonal** system in z (the variable-dz anelastic operator minus `k²`), solved by
  Thomas. O(N log N) time, **O(N) memory** (no stored factorisation), exact.
- **Nest (wall x,y):** the same with a **DST/DCT** (Neumann/Dirichlet) instead of the FFT.
- **Anelastic ρ0(z):** enters only the z-tridiagonal coefficients — still tridiagonal.
**Increment 1 — the solver (DONE 2026-09-01).** `storm_dynamics/pressure_fft.py`:
`project_anelastic_fft(u,v,w, rho0_c, rho0_wface, dx,dy, dzc,dzf, periodic_h)` makes
`div(rho0 u)=0` via FFT (periodic parent) / DCT-II (wall nest) in x,y + a batched Thomas
sweep in z (the `(0,0)`/mean mode is a *pinned* vertical Neumann solve, not zeroed — that
was the one subtlety).  Verified: on a 48²×40 stretched grid `div` → ~2e-17 (periodic) /
4e-17 (wall) in ~0.01 s, **no LU factorisation** — the direct-splu OOM is gone
(`test_fft_tridiag_anelastic_projector_divergence_free`).  NumPy + `scipy.fft`, self-contained,
does not touch `meteorological_flow.PressureSolver`.
**Generality — the tradeoff (important).** FFT+tridiag is **exact, not an approximation**,
but only for the **separable** case, and it *does* lose generality relative to the direct/
iterative solvers.  It requires: (a) **uniform horizontal spacing** (constant dx,dy — the
transform diagonalises a constant-coefficient horizontal Laplacian); (b) coefficients
**homogeneous in x,y** (anelastic ρ0=ρ0(z) only; the z-tridiagonal must be the same for every
column — dz(z) and ρ0(z) may vary in z, not in x,y); (c) **separable homogeneous BCs**
(periodic → FFT, uniform wall → DCT/DST).  The **current model and all its AMR nests sit
exactly in this class**, so there is *no* loss today — machine precision.  It would NOT
handle **terrain-following coordinates / orography, horizontal grid stretching, map
projections, x,y-varying reference states, or immersed/irregular boundaries** — the general
regimes the forecast-model generalisation (§3d, §3c) introduces.  Design: keep FFT+tridiag as
the fast low-memory solver for the separable case (route by structure), and keep the
**direct (small grids) + a multigrid-preconditioned CG (the general, non-separable fallback)**
— §3f increment 3.  Not a silent default: the solver is selected by grid structure, and a
non-separable grid must fall back, never silently use FFT+tridiag.

**Increment 2 — wire it in (NEXT).** Use it for the anelastic projection of large/fine nests
(a `NestedStormSimulation` option, or a `_pressure_method` that routes large *separable* grids
here — guarded by a `separable(grid)` check), replacing the per-level
`PressureSolver.project_anelastic`; verify a real fine nest steps stably. Then the 3rd level
(~50 m) runs without OOM.
**Increment 3 — the general low-memory fallback (DONE 2026-09-01).**
`storm_dynamics/pressure_iterative.py`: `project_anelastic_iterative(...)` assembles the
sparse anelastic Poisson operator and solves it with a **Jacobi-preconditioned CG** — general
(any sparse SPD operator: terrain, horizontal stretching, x,y-varying coefficients change only
the assembly, not the solver) and low-memory (no factorisation at all, below even ILU). Two
subtleties made CG converge, both bugs during dev: (1) the Laplacian is **negative** definite,
so we assemble/solve the **negated** operator `-L φ = -f` — CG needs SPD; on the wrong-sign
operator it stalls (this is exactly why the old `PressureSolver` Jacobi-CG "diverges" on the
stretched grid); (2) the preconditioner must be SPD — `scipy` **ILU is not symmetric** and
breaks CG (stalls at ~1e-3), so use a **Jacobi (diagonal)** preconditioner. Verified: on a
48²×40 stretched grid `div` → ~1e-11 in <1 s, periodic AND wall, and it **agrees with the
independent FFT+tridiag projector to ~2e-9** — the two cross-validate
(`test_iterative_cg_anelastic_projector_matches_fft_and_is_divergence_free`). scipy.sparse,
self-contained, does not touch `PressureSolver`. (Optional future upgrade: an algebraic-
multigrid preconditioner — `pyamg`, or generalise `poisson_mg.py` — cuts iterations on very
large stiff operators; Jacobi-CG already suffices at storm/nest scale. Beyond: multi-GPU / MPI
domain decomposition.)

---

## 4. How to resume after a session loss

1. **Read** this file, then `src/storm_dynamics/handoff.md` and `docs/amr_design.md`
   (Milestone 2 has the exact composite-projection call sites + recipe).
2. **Recall memory:** `…/memory/MEMORY.md` → `amr-composite-projection-state.md`
   (AMR state), `storm-dynamics-debugging-technique.md` (localise truncation residual by
   cell class — the technique that cracked the 2-D/3-D bugs).
3. **Verify the baseline is green:**
   `python -m pytest tests/test_storm_nesting.py -q` (expect ~24 pass, 1 skip = pyAMReX)
   and `python -m storm_dynamics.composite_poisson` (prints the 1-D/2-D/3-D + projection +
   final-assembly convergence table).
4. **See the current results:**
   `python examples/tornado_nest.py --parent-duration 1200 --u-max 18 --refine 3 --follow
   --window 300 --les-boost 1.4 --cfl 0.20 --plots --animate --device gpu`
   and `python examples/render_tornado_3d.py --device gpu --mode streamlines`.
5. **Pick up at §2a** (adaptive/dynamic regridding) — §1 (composite projection in the
   time loop) is DONE; `test_composite_projection_in_time_loop` guards it.

**Hard constraints (do not break).** Never modify `src/met_water_nucleation/_engine/`
(SHA-256 guarded; `--validate` must stay green). Do not change `meteorological_flow`
scientific behaviour or break its tests. Preserve conservation. Keep the idealised /
not-a-forecast caveats honest in every user-facing artefact.

**Git.** Work on `main` (or a branch); commit per verified step; `outputs/` is gitignored
(figures regenerable), tracked media live in `docs/media/storm/`.
