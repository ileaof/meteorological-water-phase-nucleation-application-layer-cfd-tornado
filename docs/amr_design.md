# AMR design — rigorous two-way adaptive mesh refinement for `storm_dynamics`

This document is the **honest engineering plan** for the remaining M3 work:
flux-conservative **refluxing**, a **multilevel Poisson** (pressure) solve, and
**adaptive** (dynamic) refinement. It is a *design*, not a claimed
implementation — full AMR is a separate, multi-week project, normally built on an
existing framework (AMReX, Chombo, PARAMESH, p4est). What is already implemented
in `storm_dynamics.nesting` is called out explicitly so the boundary is clear.

> **Scope check.** Delivered today (verified): one-way nesting (phases 1/2/2b),
> approximate two-way *injection* feedback (phase 3a), and **conservative
> restriction** (average-down, exact to machine precision — the first rigorous
> conservation piece). Everything below — refluxing, multilevel Poisson, adaptive
> regridding — is **not** implemented; this is the plan.

---

## 1. Where we are vs. what full AMR needs

The current solver is a **single structured grid** (`Grid` + `FlowState` + one
explicit predictor/projection `_step`). Nesting is bolted on: a second `Grid`
driven at its border by the parent, and (phase 3a) an injection feedback. That
demonstrates the mechanism but is neither conservative at the interface nor
adaptive.

Full AMR requires three things the current architecture does not have:

| Piece | What it guarantees | Current gap |
|---|---|---|
| **Refluxing** (Berger–Colella) | discrete conservation across the coarse–fine boundary | injection/relaxation only; interface fluxes not reconciled |
| **Multilevel Poisson** | the anelastic constraint `div(ρ₀u)=0` holds on the *composite* grid | pressure solved per-grid; no coarse–fine coupling |
| **Adaptive regridding** (Berger–Oliger) | refinement follows the vortex automatically | one static (or storm-following) nest, fixed refinement |

They are coupled: the projection makes the update non-local, so refluxing the
advective fluxes alone does **not** make the anelastic system conservative — the
Poisson solve must also be a composite (multilevel) solve. This coupling is the
core reason AMR here is a project, not a patch.

## 2. Conservation: restriction (done) → refluxing (to do)

### 2.1 Conservative restriction — implemented ✅

`conservative_restrict(nest, parent, spec)` replaces each covered coarse cell by
the block-mean of the `refine×refine` fine cells in its z-column. For a
**cell-aligned, matched-z** nest (`NestSpec.aligned`) this preserves the scalar
integral over the overlap **exactly** (test:
`test_conservative_restriction_preserves_overlap_integral`, error `0`). This is
the AMR "average-down" operator for the covered region.

### 2.2 Refluxing — the plan

Restriction fixes the *covered* cells; **refluxing** fixes the *coarse cells
adjacent to the interface*, whose flux across the coarse–fine face was computed
from coarse data and must be corrected to equal the sum of the fine fluxes there.

Algorithm (per coarse step, with the fine level sub-cycled `r` times):

1. Give each level a **flux register** on the coarse–fine faces.
2. On the coarse step, add the coarse face flux `F_c·Δt_c·A_c` to the register.
3. On each fine sub-step, subtract the matching fine face fluxes
   `Σ F_f·Δt_f·A_f`.
4. At sync, the register holds the mismatch `δF`; apply `∓δF / V_c` to the coarse
   cells on each side of the interface. Now the discrete flux divergence is
   single-valued across the boundary → conservation.

Code impact: the advection (`advection.advect_center_massflux`,
`momentum.py`) must **expose the face fluxes** at the interface (today they are
computed and discarded). A `FluxRegister` object per coarse–fine boundary
accumulates and applies the correction at sync points. This is mechanical but
touches every transported field.

## 3. Multilevel Poisson (the hard part)

The anelastic projection solves `∇·(1/ρ₀ ∇p') = (1/Δt)∇·(ρ₀u*)`. With a hierarchy
this must hold on the **composite** grid: fine where refined, coarse elsewhere,
with matched normal-gradient (flux) at the coarse–fine face. The current
`PressureSolver` is single-grid (cached `splu` / CG). Full AMR needs a
**geometric multigrid** on the composite grid (cell-centred, with coarse–fine
stencil modifications — a "MLMG"-style solve):

- V-cycles over the level hierarchy; the coarse–fine interface needs the standard
  cell-centred AMR stencil (ghost-cell interpolation + flux matching).
- Refluxing of the **pressure-gradient** correction at the interface, consistent
  with the velocity refluxing above, or the projected velocity is not
  divergence-free across the boundary.

This is the component that genuinely cannot be retrofitted onto the current
single-grid solver — it is a new solver. It is also the strongest argument for
**building on a framework** (AMReX ships a production `MLMG` cell-centred
multigrid with exactly this coarse–fine handling).

## 4. Adaptive regridding (Berger–Oliger)

- **Tagging.** Flag cells where a criterion fires: `|ζ|`, updraft `w`,
  `|∇θ|`/`|∇q|`, or an error estimator (Richardson extrapolation between levels).
- **Clustering.** Group tagged cells into a small set of rectangular patches
  (Berger–Rigoutsos), with a buffer so features don't leave a patch between
  regrids.
- **Regridding cadence.** Every `k` coarse steps: re-tag, re-cluster, create new
  patches, **conservatively interpolate** (prolong) coarse→new-fine, copy where
  patches overlap old ones, destroy stale patches.
- **Nested time-stepping.** Advance level `ℓ`, then recursively sub-cycle level
  `ℓ+1` `r` times, then sync (restrict + reflux + composite-project).
- **Data model.** Replace the single `Grid`/`FlowState` with a **level hierarchy**
  (`GridHierarchy` of `LevelData`), patch metadata, and a level-recursive
  `advance(level)` in place of the flat `_step`.

## 5. Required refactors of the current code

| Module | Change |
|---|---|
| `grid.py` | add a `GridHierarchy` (levels, patch boxes, coarse–fine boundaries) |
| `state.py` | `LevelData` (per-level fields) replacing the single `FlowState` |
| `advection.py`, `momentum.py` | return interface face fluxes for the `FluxRegister` |
| `pressure_solver.py` | a composite **multigrid** solve (MLMG), not single-grid |
| `core.py` | level-recursive `advance(level)` (sub-cycling + sync + reflux) |
| `nesting.py` | promote to a hierarchy manager; regridding, prolong/restrict, registers |

## 6. Milestones & effort (honest)

1. **Refluxing on a static 2-level hierarchy** (advection only; frozen Poisson):
   expose fluxes, `FluxRegister`, sync — verifiable by exact tracer conservation
   across the interface. **✅ Prototyped & verified** in pure NumPy
   (`storm_dynamics.amr.TwoLevelReflux`, `python -m storm_dynamics.amr`): total-mass
   drift over 40 steps is `1.1e-4` *without* refluxing and `2.0e-16` *with* it
   (test `test_amr_refluxing_conserves_across_interface`). This reference logic
   ports directly onto a framework `FluxRegister`.
2. **Composite multigrid Poisson** on the static 2-level grid — *~3–5 weeks*
   (the crux); verify divergence-free across the interface. **Solver kernel done**
   (`storm_dynamics.poisson_mg`): a geometric-multigrid V-cycle (red-black GS,
   full-weighting restriction, bilinear prolongation) that converges h-independently
   (~8 V-cycles to 1e-10 at any resolution) and is 2nd-order accurate on a
   manufactured solution — this is the "bring-our-own-solver" kernel needed because
   pyAMReX exposes no MLMG. **Composite coarse-fine interface stencil done in 1-D**
   (`storm_dynamics.composite_poisson`): a fine patch replacing part of a periodic
   coarse grid, coupled with a **2nd-order ghost** (quadratic through two fine + one
   coarse cell, at the ghost point symmetric to the fine boundary about the
   interface) and a **single-valued (conservative)** interface flux. Verified
   2nd-order on a manufactured solution (error falls ~4x per 2x refinement across
   the interface, `test_composite_poisson_1d_second_order_across_interface`).
   **2-D and 3-D done too** (`solve_2d`, `solve_3d`): a rectangular fine patch / fine
   box in a periodic coarse grid, with a tangential coarse interpolation for the ghost
   (linear in 2-D, **bilinear** in 3-D) and the interface flux oriented
   `d(phi)/d(+axis)` and single-valued — verified 2nd-order *including the patch
   corners (2-D) and box edges/corners (3-D)* (`test_composite_poisson_2d_*`,
   `test_composite_poisson_3d_*`, ratios ~4.0). (The bug that made a first 2-D attempt
   blow up was a flux-orientation sign on the L/B (−axis) edges, found by localising
   the truncation residual to the edge cells; the 3-D generalisation then worked first
   try.) This composite operator is now **wired into a two-level MAC projection**
   (`project_divergence_2d`, `test_composite_projection_2d_divergence_free_across_interface`):
   a face-flux velocity `u*` on the composite grid → `f = div(u*)` → `solve_2d` gives
   `p` → `u = u* − grad(p)`. Because the divergence, gradient and the composite
   Laplacian all use the **same** single-valued interface flux (`L = div·grad`),
   `div(u) = f − L p` vanishes to the solve tolerance (~1e-13) **including at the
   coarse-fine interface** — verified on a random `u*` with a divergence/gradient
   built independently of the solver (self-validating). **3-D projection done too**
   (`project_divergence_3d`, `test_composite_projection_3d_*`): `max|div u|` ~1e-13 at
   the interface (a subtle bug — a `mean_iface` closure over the *uncorrected* face
   field left the coarse interface-adjacent cells divergent while fine cells were
   already clean — was localised by splitting interior-coarse vs interface-adjacent
   and fixed by threading the field through).

   **This is the anelastic projection — no new algorithm remains.** The storm's
   projection is `u = u* − grad(p)/ρ0` with `div(ρ0 u)=0`; in mass-flux variables
   `m = ρ0 u` that is exactly `m = m* − grad(p)`, `div(m)=0` — the projection above
   with the face field read as the anelastic mass flux `ρ0 u`. The density weight
   cancels in the divergence constraint and only re-enters when recovering `u = m/ρ0`.

   **Remaining step 2 — integration into `NestedStormSimulation` (plumbing, not math).**
   Today the parent (`StormSimulation._step`, `core.py:213`) and the nest
   (`NestedStormSimulation._step`, `nesting.py:341`) each call
   `pressure.project_anelastic(...)` *independently* on their own grid. The composite
   projection replaces those two calls with **one** solve over both levels' mass
   fluxes:
   1. after the momentum predictor, gather the two levels' staggered face mass fluxes
      `ρ0 u*` (parent coarse cells outside the nest footprint + nest fine cells) into
      the composite layout (`ci0..cj1` = the nest's placement in the parent);
   2. `f = div(ρ0 u*)` (composite divergence) → `solve_2d/solve_3d` → `p`;
   3. correct `u = u* − grad(p)/ρ0` on both levels with the shared single-valued
      interface flux, write back into each `FlowState`.
   The three genuine engineering items are (a) the **solid-wall vertical BC** (and the
   nest's lateral wall/relaxation BC) in place of the periodic wrap the reference
   solver assumes — a Neumann row on the outer boundary, the interface stencil
   unchanged; (b) extracting/reinserting the staggered `u,v,w` faces from the two
   `FlowState`s; (c) the **stretched physical z-grid** (`hf`, `hc` become per-level
   metric factors). None of these change the verified interface stencil; they are the
   production wiring that makes the two-level solve act on the real storm.

   **Item (a) done** (`solve_2d(..., periodic=False)`, `project_divergence_2d(...,
   periodic=False)`): the solid-wall Neumann BC — a boundary coarse cell drops its
   outward face (`dphi/dn=0`), the outer boundary velocity carries no flux and gets no
   pressure correction, and the interior interface stencil is untouched. Verified
   2nd-order on `phi=cos(pi x)cos(pi y)` (which satisfies `dphi/dn=0` on the walls;
   ratio 4.00) and the wall projection is divergence-free across the interface to
   ~1e-13 (`test_composite_solid_wall_bc_second_order_and_projection`).

   **Item (b) done** (`composite_project_massflux_2d`): the face-array bridge that reads
   the storm's staggered C-grid **mass fluxes** in their native convention (`u:(nc+1,nc)`,
   `v:(nc,nc+1)` on the parent; `(nfx+1,nfy)`/`(nfx,nfy+1)` on the nest), maps them to the
   composite decomposition, projects (`div` → `solve_2d` → correct), writes the corrected
   fluxes back, and **refluxes** the parent's interface faces to the single-valued fine
   mean. Verified by an *independent* recomputation of `div(m)` straight from the
   written-back arrays: ~1e-13 across the interface, with the solid-wall BC (and periodic)
   — `test_composite_massflux_bridge_storm_arrays_divergence_free`. It works in mass-flux
   variables `m = ρ0 u`, so it is exactly the anelastic constraint; the caller forms
   `ρ0 u*` (as `PressureSolver.project_anelastic` already does) and recovers `u = m/ρ0`.
   Remaining: item (c) the stretched physical z-grid metric (per-level `hf`,`hc`), and
   calling this bridge from `NestedStormSimulation._step` in place of the two independent
   `project_anelastic` calls (the 3-D analogue extends `composite_project_massflux_2d` the
   same way `solve_2d`→`solve_3d` did).
3. **Nested time-stepping + sync** wiring both together — *~2 weeks*.
4. **Adaptive regridding** (tag/cluster/regrid, prolong/destroy) — *~2–3 weeks*.
5. **Hardening**: multi-level (3+), moving/merging patches, load balancing if
   parallel — open-ended.

Total: a **multi-month** effort from scratch; markedly less on a framework.

## 7. Recommendation

Do **not** hand-roll composite multigrid + AMR bookkeeping. Build the hierarchy on
**AMReX** (Python via `pyamrex`, or C++): it provides the box/patch data model,
`FluxRegister`, the cell-centred `MLMG` multigrid, regridding and
tag/cluster — exactly the three hard pieces above — and is battle-tested for
atmospheric/CFD codes. The physics we already have (advection, LES, buoyancy,
microphysics, the nucleation kernel) ports as the per-level right-hand side; the
AMR infrastructure comes from the framework. That turns "several months from
scratch" into "port the physics onto a proven AMR core."

### Getting pyAMReX built (WSL2, reproducible)

`scripts/build_pyamrex_wsl.sh` is the turnkey recipe, with every blocker found on a
fresh WSL2/Ubuntu baked in:

- **Python** — pyAMReX ≥ 26.8 needs **Python ≥ 3.11**; use a conda (Miniforge)
  **3.12** env. Ubuntu 22.04's system Python 3.10 is too old; the Miniforge *base*
  is 3.14, too new. 3.12 is the sweet spot.
- **Memory** — the pybind11 bindings are memory-hungry; a full-parallel build
  **OOM-kills `cc1plus`** on a 7 GB WSL. Build with few jobs (`JOBS=2`), or raise
  WSL RAM in `%UserProfile%\.wslconfig` (`[wsl2]` → `memory=12GB`) then
  `wsl --shutdown`.
- **Toolchain present out of the box on WSL2/Ubuntu 22.04**: g++ 11.4, make;
  install `cmake`/`ninja` via `pip` (no sudo). The GPU (e.g. RTX 4050) is visible
  in WSL via `nvidia-smi`, so an `AMReX_GPU_BACKEND=CUDA` build is possible once the
  CPU build is working.

### pyAMReX validation (what the default build exposes)

Built and imported (pyAMReX `26.08`, CPU, 3D) and probed on WSL:

- ✅ **Data model**: `MultiFab`, `Geometry`, `BoxArray`, `DistributionMapping`
  work — a `set_val(3)` / `sum` / `mult` round-trip is exact.
- ✅ **AMR hierarchy**: `AmrCore`, `AmrMesh`, `AmrInfo`, `AmrParGDB` are present —
  the regridding / level-hierarchy backbone the port builds on. A `Poisson` class
  is also exposed.
- ❌ **Not exposed by pyAMReX at all** (confirmed by inspecting the pyAMReX
  sources — there is **no `LinearSolvers` binding directory** and **no `MLMG` /
  `MLLinOp` / `FluxRegister` binding** anywhere under `src/`): the composite
  multigrid solve and the flux register. AMReX's C++ library *has* them
  (`AMReX_LINEAR_SOLVERS=ON`), but pyAMReX does **not** wrap them for Python. **A
  rebuild does not help** — there is nothing to enable. To use MLMG/FluxRegister
  from the framework you must either (a) write that part in **C++/AMReX** (the
  Python bindings are a deliberate subset), or (b) add the pybind bindings to
  pyAMReX upstream.

**Revised path (important).** The Python (pyAMReX) route gives us the **data model
+ AMR hierarchy (`AmrCore`) + our physics** — proven by the scaffold — but **not**
the framework's Poisson solver or refluxing. So the pragmatic port is: keep
pyAMReX for the hierarchy/data and **bring our own** solver and reflux — we already
have flux-conservative refluxing, conservative restriction and prolongation in
`storm_dynamics.amr` (verified), and a composite **multigrid Poisson we write
ourselves** (design-doc Milestone 2) closes the last gap **without needing
pyAMReX's MLMG**. The alternative, for a production code, is to write the AMR app
in C++/AMReX. Either way, the earlier "just rebuild pyAMReX with the solver
bindings" step is a dead end.

### Port scaffold — validated ✅

`storm_dynamics.amr_port.AmrexAdvect3D` is the first step of the physics port: the
tracer lives in an **AMReX `MultiFab`**, AMReX does the ghost exchange
(`fill_boundary` with periodicity), and **our** flux-form NumPy stencil (accessing
the MultiFab via `to_numpy()`, layout `(nx+2g, ny+2g, nz+2g, ncomp)`) computes the
update. Verified on WSL (`python -m storm_dynamics.amr_port`): a 32³ periodic
advection over 50 steps has **total-mass drift `0.0`** — the "AMReX infrastructure
+ our RHS" binding the full port rests on. The two-level, refluxed version reuses
the verified `storm_dynamics.amr` operators once the `FluxRegister`/`MLMG` bindings
are enabled. (Note: access the MultiFab via `mf.to_numpy()`; the per-`MFIter`
`array(mfi).to_numpy()` path segfaulted in this build — use the whole-box view.)

## 8. Verification plan

- **Conservative restriction**: overlap integral preserved (done, error `0`).
- **Refluxing**: a passive tracer is conserved to machine precision across a
  static coarse–fine interface over many steps.
- **Composite projection**: `max|div(ρ₀u)|` at the interface ~ solver tolerance.
- **Free-stream preservation**: a uniform flow stays uniform through a refined
  patch (no spurious sources at the boundary).
- **Regridding**: total mass/energy unchanged across a regrid (conservative
  prolongation).

## References

- Berger, M. J., & Oliger, J. (1984). *Adaptive mesh refinement for hyperbolic
  partial differential equations.* J. Comput. Phys., 53, 484–512.
- Berger, M. J., & Colella, P. (1989). *Local adaptive mesh refinement for shock
  hydrodynamics.* J. Comput. Phys., 82, 64–84.
- Berger, M. J., & Rigoutsos, I. (1991). *An algorithm for point clustering and
  grid generation.* IEEE Trans. SMC, 21, 1278–1286.
- Almgren, A. S., et al. (1998). *A conservative adaptive projection method for the
  variable-density incompressible Navier–Stokes equations.* J. Comput. Phys., 142,
  1–46 (AMR + projection — the closest analogue to our anelastic case).
- Zhang, W., et al. (2019). *AMReX: a framework for block-structured AMR.*
  J. Open Source Softw., 4(37), 1370.
