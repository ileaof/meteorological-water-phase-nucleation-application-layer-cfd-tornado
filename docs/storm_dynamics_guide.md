# storm_dynamics — idealised rotating supercell / tornadogenesis core

> A styled, self-contained HTML edition of this guide is at
> [`docs/MANUAL_storm_dynamics.html`](MANUAL_storm_dynamics.html) (same format as
> `docs/MANUAL.html`).

`storm_dynamics` is a **fork of the dynamical core** of `meteorological_flow` that
adds the physics a deep-convective storm needs in order to **rotate** — and,
under a curved hodograph with surface friction and a cold pool, to spin up
**near-surface vertical vorticity**, the proxy for tornadogenesis. It reuses the
repository's already-validated physics (thermodynamics, bulk microphysics,
staggered grid, anelastic pressure projection, conservative scalar transport, the
nucleation kernel) unchanged, and rewrites **only** the momentum/dynamics.

> **Read this first — scope and honesty.** This is **idealised simulation, not
> operational forecasting.** There is **no data assimilation, no real-event
> initial or boundary condition, and no observational verification.** The
> deliverable is a *rotational dynamical core* validated against the **classical**
> idealised-supercell results (Klemp & Wilhelmson 1978; Weisman & Klemp 1982;
> Rotunno & Klemp 1985). Nothing here forecasts a real tornado, a real storm, or
> a real event, and the vortex is **under-resolved** at demonstration resolution.

---

## Why the demonstration core cannot rotate

The `meteorological_flow` storm (`meteorological_flow/simulation.py`) makes three
choices that are correct for a *demonstration of nucleation-coupled convection*
but fatal to rotation:

1. **Momentum advection is switched off** (a documented v1 simplification). With
   no `(u·∇)u`, there is no **tilting** of horizontal (shear) vorticity into the
   vertical and no **stretching** of vertical vorticity by the updraft — the two
   terms of the vertical-vorticity equation that *make* a mesocyclone.
2. **Rayleigh drag** (`gamma_damp`) linearly relaxes the velocity every step —
   it would damp the very rotation we want.
3. A hard **velocity clip** at ±120 m/s caps peak winds unphysically.

`storm_dynamics` fixes all three.

## What the fork adds (and what it reuses)

| # | Piece | Module | Reused from the repo |
|---|-------|--------|----------------------|
| 1 | **Conservative flux-form momentum advection** (staggered C-grid, MUSCL/minmod) — the enabling term (tilting + stretching) | `momentum.py` | reconstruction style of `advection.advect_center_massflux`; `_minmod` |
| 2 | **f-plane Coriolis** on the perturbation wind (~36 °N) | `coriolis.py` | — |
| 3 | **LES (Smagorinsky) subgrid closure** replacing Rayleigh drag + clip | `turbulence.py` | grid stencils, `thermodynamics` |
| 4 | **Surface bulk-drag law** (corner-flow friction) | `surface_drag.py` | — |
| 5 | **Curved (quarter-circle) hodograph** → storm-relative helicity | `soundings.py` | `base_state.weisman_klemp` thermodynamics |
| 6 | **Evaporative cold pool** (baroclinic vorticity source) + microphysics; optional **nucleation-kernel coupling** | `core.py` | `microphysics_coupling.MicrophysicsCoupler`, `precip_microphysics`, `nucleation_adapter` + validated kernel |
| 7 | **Rotation diagnostics** (ζ, updraft helicity, SRH, trackers) | `rotation.py` | grid stencils, `diagnostics` budgets |

The grid, the anelastic Chorin projection (`pressure_solver.project_anelastic`,
already periodic-capable), conservative scalar transport, moist buoyancy, the
bulk microphysics and the validated nucleation kernel are imported **unchanged**.
The immutable `met_water_nucleation._engine` is never touched, so
`met-water-nucleation --validate` (tests [1]–[21] + the SHA-256 checksum guard)
stays green.

## Numerical formulation

- **Core:** anelastic (`div(ρ₀u)=0`), reference density `ρ₀(z)` from the sounding
  — appropriate for a 10–16 km storm column (Boussinesq is stretched past
  validity there).
- **Time stepping:** explicit predictor + pressure projection, adaptive CFL that
  now includes the advected momentum.
- **Momentum advection:** flux/divergence form on each staggered control volume,
  2nd-order MUSCL (minmod) reconstruction. Because the transporting velocity is
  the projection's discretely divergence-free field, the flux form equals the
  advective form while telescoping to the boundary flux — domain-integrated
  momentum is conserved (to the boundary flux) at any intensity, **with no
  velocity clip**. The only remaining bound is an *extreme numerical guard*
  (`v_guard`, default 150 m/s), documented as a safety rail that should never bite
  in a resolved run.
- **LES closure:** `K_m = (C_s Δ)² |S| · f_stab(Ri)`, `Δ = (ΔxΔyΔz)^{1/3}`, with a
  Lilly stability correction that damps mixing in statically stable layers;
  scalars use `K_h = K_m / Pr_t`. Momentum diffusion is applied to the
  *perturbation* so the environmental hodograph is not mixed away.
- **Boundaries:** periodic lateral (x, y) — the mean-wind shear is ingested
  through them; free-slip/drag bottom; damping-layer top.
- **Conservation:** water and energy close to the same order as the demonstration
  storm; the continuity residual comes from the projection (~0), not from
  limiters. Water/energy are conserved up to the **boundary flux** through the
  damping top — the expected (and documented) non-conservation channel, exactly
  as in the parent solver.

## Milestones

### M1 — rotating supercell (storm splitting + mid-level mesocyclone)
Unidirectional shear over the Weisman-Klemp sounding. A saturated warm bubble
grows into a deep updraft that **splits** into left- and right-moving cells
(cyclonic **and** anticyclonic vertical vorticity of comparable magnitude) with a
**mid-level mesocyclone** (ζ at ~3–6 km above the mesocyclone threshold). This is
the Klemp-Wilhelmson / Weisman-Klemp sanity test that the rotation dynamics work.

*Demonstration result (32×32×40, Δx≈1.25 km, 40 min):* w_max ≈ 20 m/s;
ζ_max/ζ_min ≈ +0.050 / −0.050 s⁻¹ (splitting); mid-level mesocyclone ≈
1.1×10⁻² s⁻¹; updraft helicity ≈ 230 m²/s² (peak ≈ 340).

### M2 — low-level rotation (tornadogenesis proxy)
Curved (quarter-circle) hodograph → strong storm-relative helicity, plus the
**surface drag law** and the **evaporative cold pool** (the baroclinic vorticity
source). Near-surface vertical vorticity develops on the cold-pool / forward-flank
interface — the target of the project. This is a *proxy*: the ~100 m tornado
vortex itself is not resolved at demonstration Δx.

*Demonstration result (curved hodograph + drag, vertical levels clustered near
the surface, dz₀≈100 m, 40 min):* the near-surface ζ **spins up and is sustained**
(~3.3×10⁻³ s⁻¹, roughly steady to the end), whereas an identical run with a
straight hodograph and no surface drag peaks lower (~2.5×10⁻³) and **decays** to
~1.4×10⁻³ — i.e. the curved hodograph + friction hold ~2.4× the low-level rotation
of the control. The updraft stays physical (w_max ≈ 18 m/s) and water/mass
conservation is unchanged. The near-surface ζ magnitude is a **proxy**, not a
tornado wind speed — the vortex is under-resolved (see the caveats below).

> **Stability note.** A strong low-level-SRH environment on a coarse grid will
> over-amplify the updraft (grid-scale w blow-up) unless the near-surface layer is
> resolved. The shipped M2 config therefore clusters vertical levels
> (`z_stretch`), moderates the hodograph (`U_max≈18`, `z_turn≈2 km`) and uses a
> slightly larger Smagorinsky constant. These are documented demonstration
> choices, not physical tuning of the result.

### M3 — fine vortex *(phases 1 & 2: one-way nesting delivered)*
Nested refinement to resolve the low-level vortex at finer scale. **Full AMR**
(adaptive, block-structured, two-way) remains a separate project; delivered here
is **one-way nesting** (`storm_dynamics.nesting`), the classical idealised-tornado
approach: mature the storm on the coarse **parent**, interpolate the updraft /
low-level-rotation region onto a finer **nest** (trilinear, exact for linear
fields), and integrate the nest — reusing the whole solver (momentum advection,
LES, drag, projection, microphysics) — with its border relaxed toward the parent
(Davies-style nudging).

**Phase 1 — static / frozen-parent boundary.** The border is nudged toward the
parent captured at the nest start time. *Result (parent Δx≈1.3 km → nest Δx≈0.44
km, 3×, 120 s):* the nest inherits and sustains the updraft and **intensifies the
near-surface ζ ~2.4×** over the window while conserving — the vortex sharpening
under refinement. It is valid only for a **short window** (~2–3 min); beyond that
the frozen border stops feeding fresh inflow, the storm decays, and the ζ maximum
drifts onto the nest edge (a sponge artefact — read the *interior*-masked ζ, which
excludes the relaxation band). Run: `examples/tornado_nest.py`.

**Phase 2 — concurrent / time-evolving boundary** (`run_concurrent_nest`, example
`--concurrent`). The parent **steps alongside** the nest; each parent step its
state is re-interpolated onto the nest and the nest sub-cycles (finer dt) with the
border relaxed toward the parent target **interpolated linearly in time**. Fresh
inflow keeps entering, so the nest is **sustained as long as the parent drives it**
(no frozen-boundary decay) and stays stable (a modest extra Smagorinsky boost +
tighter CFL guard the sharpening vortex). A *fixed* nest still loses a cell that
translates out of its region.

**Phase 2b — storm-following nest** (`run_concurrent_nest(follow=True)`, example
`--follow`). The nest runs in the **storm-relative frame**: a Galilean shift by the
storm motion **C** (Bunkers) makes the cell quasi-stationary, and the *sampled*
parent region slides at C so the storm stays centred in the fixed nest. Vorticity
and w are invariant under the constant shift, so the diagnostics are unchanged.
*Result (3× nest, C≈(9, −18) m/s, 400 s):* the updraft is not just sustained but
**intensifies from ~7 m/s (coarse parent) to ~19 m/s (fine nest), growing through
the whole window**, with excellent conservation (water ≈ −0.2%, |div| ≈ 5×10⁻⁵) —
the finer grid resolving a much stronger updraft that a fixed nest would have lost.

**Phase 3a — approximate two-way feedback** (`run_concurrent_nest(two_way=True)`,
example `--two-way`). Closes the loop: after each parent step the nest's finer
solution is blended back (converted to the ground frame) onto the parent cells it
overlaps, with a taper that fades to 0 at the overlap edge. Now the parent is
**improved by the nest** — *result:* with feedback the parent updraft ends at
~9 m/s vs ~6 m/s in an identical run without it (the nest's stronger, better-
resolved updraft propagates back), the loop stays **stable**, and water conserves
to ~−0.1%. This is **injection** feedback (sample fine at coarse points), **not**
rigorous flux-conservative refluxing.

> **Honest scope / remaining M3.** Delivered: phases 1, 2, 2b (one-way) and **3a
> (approximate two-way)**. What's left for *full* AMR is a genuine separate
> project: **rigorous conservation at the coarse–fine interface** (Berger–Colella
> refluxing) with a **multilevel pressure (Poisson) solve**, and **adaptive,
> dynamic** refinement (block-structured patches created/moved/destroyed by a
> refinement criterion, with nested time-stepping) — typically built on a
> framework (AMReX, Chombo, p4est). The storm-following frame also uses a
> *constant* C (a real storm's motion varies, so the near-surface ζ maximum can sit
> near the nest edge — read the *interior*-masked ζ), and reaching O(10–100 m)
> needs much higher refinement (`--refine 5–6` / finer parent, far costlier).

## What this model **can** claim

- The rotational **dynamics** are present and correct in kind: tilting and
  stretching produce a mid-level mesocyclone and storm splitting from a sheared
  environment (M1), reproducing the classical idealised-supercell behaviour.
- Under a curved hodograph with surface friction and an evaporatively driven cold
  pool, **near-surface vertical vorticity** develops (M2) — the qualitative
  tornadogenesis pathway.
- Conservation (water, energy, mass continuity) holds to the parent solver's
  levels.

## What this model **cannot** claim

- **It does not forecast** any real storm, tornado, or event. No assimilation, no
  observed initial state, no verification against observations.
- **It does not resolve the tornado vortex.** At Δx of ~1 km the near-surface ζ is
  a *proxy*; peak ζ, updraft-helicity and wind magnitudes are **indicative, not
  quantitative**. Quantitative low-level rotation needs O(100 m) or finer grids
  (M3) and is explicitly future work.
- **Idealisations:** f-plane (constant Coriolis), analytic sounding, single
  warm-bubble trigger, single-moment bulk microphysics, LES closure with a fixed
  Smagorinsky constant. These are the standard idealised-modelling choices, not
  deficiencies to be hidden — but they bound what the numbers mean.

## Non-goals (explicit)

Data assimilation, real-event initial/boundary conditions, operational
forecasting, and verification against observations are **out of scope** for this
phase. The deliverable is an **idealised rotating-storm dynamics** simulator built
on the repository's already-validated physics.

## How to run

```bash
# M1 — rotating supercell (storm splitting); --plots writes the rotation figures
PYTHONPATH=src python examples/supercell_tornadogenesis.py --scenario supercell --plots

# M2 — low-level rotation (curved hodograph + drag + cold pool)
PYTHONPATH=src python examples/supercell_tornadogenesis.py --scenario tornadogenesis \
    --nx 40 --ny 40 --nz 48 --duration 3600 --plots

# optional: couple the validated nucleation kernel as the microphysics embryo
# source (eq39 pathway), exactly as meteorological_flow does (builds a lookup
# table -> slower); off by default
PYTHONPATH=src python examples/supercell_tornadogenesis.py --scenario supercell \
    --kernel-nucleation

# M3 phase 1 — static nested-grid refinement of the low-level vortex
PYTHONPATH=src python examples/tornado_nest.py --refine 3 --window 120 --plots --animate

# declarative configs
python -c "from storm_dynamics.config import storm_config_from_yaml as L; \
           from storm_dynamics.core import StormSimulation as S; \
           print(S(L('configs/storm_supercell.yaml')).run()['rotation'])"

# tests
PYTHONPATH=src python -m pytest tests/test_storm_dynamics.py tests/test_storm_milestones.py -q
```

## Compute backend (CPU / GPU)

The core is backend-agnostic — every hot-loop module (`momentum`, `turbulence`,
`coriolis`, `surface_drag`, `rotation`) is written against `grid.xp`, so it runs
on **CPU (NumPy)** or **GPU (CuPy)** with identical results. Select with
`--device` (example) or `device=` / `device:` (config):

```bash
python examples/supercell_tornadogenesis.py --scenario supercell --device cpu   # default
python examples/supercell_tornadogenesis.py --scenario supercell --device gpu   # NVIDIA GPU + CuPy
python examples/supercell_tornadogenesis.py --scenario supercell --device auto  # GPU if available, else CPU
```

- `cpu` (default) — preserves the original behaviour; no GPU probing.
- `gpu` — requires an NVIDIA GPU + CuPy (`pip install .[gpu]`); **fails loudly** if absent.
- `auto` — tries GPU, falls back to CPU with a logged reason.

GPU parity is regression-tested (`test_storm_gpu_matches_cpu`, skipped when no GPU
is present). **Note:** at demonstration grid sizes (≤64k cells) the direct sparse
pressure solve stays on the host and CPU is usually *faster*; GPU pays off only at
larger grids (>64k cells), where the pressure Poisson solve uses CG on the GPU.

## References

- Klemp, J. B., & Wilhelmson, R. B. (1978). *The simulation of three-dimensional
  convective storm dynamics.* J. Atmos. Sci., 35, 1070–1096.
- Weisman, M. L., & Klemp, J. B. (1982). *The dependence of numerically simulated
  convective storms on vertical wind shear and buoyancy.* Mon. Wea. Rev., 110,
  504–520.
- Rotunno, R., & Klemp, J. B. (1982). *The influence of the shear-induced pressure
  gradient on thunderstorm motion.* Mon. Wea. Rev., 110, 136–151.
- Rotunno, R., & Klemp, J. B. (1985). *On the rotation and propagation of simulated
  supercell thunderstorms.* J. Atmos. Sci., 42, 271–292.
- Davies-Jones, R. (1984). *Streamwise vorticity: The origin of updraft rotation
  in supercell storms.* J. Atmos. Sci., 41, 2991–3006.
- Bryan, G. H., & Fritsch, J. M. (2002). *A benchmark simulation for moist
  nonhydrostatic numerical models.* Mon. Wea. Rev., 130, 2917–2928 (the CM1-style
  nonhydrostatic core taken as a reference formulation).
- Lilly, D. K. (1962). *On the numerical simulation of buoyant convection.* Tellus,
  14, 148–172 (Smagorinsky/Lilly SGS closure).
- Bunkers, M. J., et al. (2000). *Predicting supercell motion using a new
  hodograph technique.* Wea. Forecasting, 15, 61–79.
