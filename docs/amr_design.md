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
   (the crux); verify divergence-free across the interface.
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
