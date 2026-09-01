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

**NEXT: §2b increment 2** — a concurrent multi-level driver (each level sub-cycles under the
one above, with inter-level refluxing) so the finest level is *sustained*, not a short
static snapshot. ★  (§2a adaptive regridding + §2b increment 1 = recursive nesting landed
2026-09-01, both verified; `--levels 2` in the renderer reaches ~150 m = funnel scale.)

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
- **Increment 2 — concurrent multi-level driver + inter-level refluxing (NEXT).** Today the
  extra level runs *static* (phase-1, frozen L1 boundary) for a short window. Generalise the
  concurrent driver so each level sub-cycles in time under the level above (recursive
  `run_concurrent_nest`), feeding time-evolving boundaries down and **refluxing mass/momentum
  up** across each interface (`amr.TwoLevelReflux` per level pair). Then the finest level is
  *sustained*, not just a short snapshot.

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

**3f. Performance & scaling.** Multi-GPU / MPI domain decomposition; the direct pressure
solve → scalable multigrid/CG at scale (`poisson_mg.py` is the kernel).

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
