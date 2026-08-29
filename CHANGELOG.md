# Changelog

All notable changes to this repository's organization are documented here.
Scientific behaviour is unchanged across the 1.0.0 reorganization (the engine
moved byte-identically; SHA-256 checksums preserved).

## [Unreleased] — storm_dynamics (rotating supercell / tornadogenesis core)

### Added
- New package `src/storm_dynamics/` — an **idealised rotating deep-convection**
  core (supercell → mesocyclone → near-surface rotation), forked from the
  `meteorological_flow` dynamical core. Adds: conservative flux-form staggered
  momentum advection (tilting + stretching of vorticity), f-plane Coriolis, a
  Smagorinsky LES closure (replacing the demonstration Rayleigh drag + velocity
  clip), a surface bulk-drag law, a curved (quarter-circle) hodograph, and
  rotation diagnostics (ζ, updraft helicity, storm-relative helicity, trackers).
  Reuses the grid, anelastic pressure projection, conservative scalar transport,
  moist buoyancy, bulk microphysics (evaporative cold pool) and the validated
  nucleation kernel **unchanged**.
- `tests/test_storm_dynamics.py` (16 unit tests) and `tests/test_storm_milestones.py`
  (M1 storm-splitting + mid-level mesocyclone; M2 sustained low-level rotation).
- Optional in-storm coupling of the validated nucleation kernel as the microphysics
  embryo source (`build_storm_config(couple_nucleation=True)` / example
  `--kernel-nucleation`), exactly as `meteorological_flow` does; off by default.
- CPU/GPU compute-backend selection for the storm core (`build_storm_config(device=)`,
  config `device:`, example `--device cpu|gpu|auto`; default `cpu`). The core is
  backend-agnostic (all hot loops use `grid.xp`); GPU parity with CPU is
  regression-tested (`test_storm_gpu_matches_cpu`, skipped without a GPU).
- **M3 one-way nested-grid refinement** (`storm_dynamics.nesting`,
  `examples/tornado_nest.py`, `tests/test_storm_nesting.py`): mature the storm on the
  coarse parent, interpolate the vortex region onto a finer nest (exact trilinear
  interpolation), and integrate the nest — reusing the whole solver — with
  Davies-style border relaxation.
  - **Phase 1** (frozen parent boundary): the finer grid intensifies the
    near-surface ζ ~2.4× over a short window (Δx 1.3 km→0.44 km, 3× finer) while
    conserving; valid only ~2–3 min before the frozen border decays the storm.
  - **Phase 2** (`run_concurrent_nest`, example `--concurrent`): the parent steps
    alongside the nest and feeds time-evolving boundaries, sustaining the nest as
    long as the parent drives it (a modest SGS boost + tighter CFL keep it stable).
  - **Phase 2b** (`run_concurrent_nest(follow=True)`, example `--follow`):
    storm-following nest — a storm-relative (Galilean) frame keeps the cell centred
    while the sampled parent region slides at the storm motion C, so the updraft is
    sustained AND intensified (~7→19 m/s over 400 s, growing the whole window) with
    water ≈ −0.2% — the finer grid resolving a much stronger updraft a fixed nest
    would have lost.
  - **Phase 3a** (`run_concurrent_nest(two_way=True)`, example `--two-way`):
    approximate two-way feedback — the nest's finer solution is blended back onto the
    parent overlap (converted to the ground frame, tapered at the edge), so the parent
    is improved by the nest (parent updraft ~6→9 m/s vs a no-feedback control), the
    loop stays stable, water ≈ −0.1%. Injection feedback, **not** rigorous refluxing.
  - **Conservative restriction** (`conservative_restrict`, `NestSpec.aligned`): the
    first rigorous conservation piece — average-down of a cell-aligned, matched-z
    nest preserves the overlap scalar integral **exactly** (machine precision; test
    `test_conservative_restriction_preserves_overlap_integral`).
  - `interior_near_surface_zeta` reports the physical interior vortex, excluding the
    boundary sponge. **Not** implemented (see `docs/amr_design.md` for the plan):
    flux-conservative refluxing, a multilevel-Poisson composite solve, and
    adaptive/dynamic regridding — full AMR, a separate multi-month project.
- `docs/amr_design.md` — the rigorous engineering plan for full two-way adaptive
  AMR (refluxing, multilevel Poisson, Berger–Oliger regridding, required refactors,
  milestones/effort, framework recommendation, verification plan).
- `storm_dynamics/composite_poisson.py` — **composite two-level Poisson interface
  stencil** (the AMR-projection crux), in **1-D, 2-D and 3-D**: a fine patch/box in a
  periodic coarse grid, coupled with a 2nd-order ghost (quadratic normal + tangential
  coarse interpolation — linear in 2-D, **bilinear** in 3-D) and a single-valued
  (conservative) interface flux oriented `d(phi)/d(+axis)`, assembled sparsely and
  solved directly. Verified 2nd-order on a manufactured solution in every dimension,
  **including the patch corners (2-D) and box edges/corners (3-D)** (error ∝ h²;
  ratios ~4.0). A first 2-D attempt blew up (1st-order/non-conservative); the fix was
  a flux-orientation sign on the L/B (−axis) edges, found by localising the truncation
  residual to the edge cells — the 3-D generalisation then worked first try.
  **Wired into a two-level MAC projection** in 2-D and 3-D (`project_divergence_2d`,
  `project_divergence_3d`): a random face-flux velocity is made discretely
  divergence-free by `div(u*)` → `solve_2d/3d` → `u = u* − grad(p)`; because the
  divergence, gradient and Laplacian share the same single-valued interface flux
  (`L = div·grad`), `max|div u|` falls to the solve tolerance (~1e-13) **including at
  the coarse-fine interface cells**. The divergence/gradient are built independently
  of the solver, so the tests are self-validating
  (`test_composite_projection_{2d,3d}_divergence_free_across_interface`).
  **This is the anelastic projection across a refinement interface** — in mass-flux
  variables `m = ρ0 u`, `u = u* − grad(p)/ρ0` with `div(ρ0 u)=0` is exactly
  `m = m* − grad(p)`, `div(m)=0` (the density weight cancels in the constraint), so no
  new algorithm is needed to make the storm anelastic. Integrating it into
  `NestedStormSimulation` (replacing the two independent `project_anelastic` calls with
  one composite solve over both levels' mass fluxes) is plumbing — wall BCs, face
  extraction, the stretched grid — documented in `docs/amr_design.md`.
  **Plumbing item (a) done:** `solve_2d`/`project_divergence_2d` take `periodic=False`
  for the **solid-wall (Neumann) BC** the storm/nest use — a boundary coarse cell drops
  its outward face, the interior interface stencil unchanged. Verified 2nd-order on
  `cos(pi x)cos(pi y)` (ratio 4.00) and the wall projection is divergence-free across
  the interface (~1e-13); `test_composite_solid_wall_bc_second_order_and_projection`.
  **Plumbing item (b) done:** `composite_project_massflux_2d` — the face-array bridge
  that reads the storm's staggered C-grid mass fluxes in their native convention
  (`u:(nc+1,nc)`, `v:(nc,nc+1)` parent; `(nfx+1,nfy)`/`(nfx,nfy+1)` nest), projects,
  writes back, and refluxes the parent's interface faces to the single-valued fine mean.
  `div(m)` recomputed independently from the written-back arrays is ~1e-13 across the
  interface (wall + periodic); `test_composite_massflux_bridge_storm_arrays_divergence_free`.
  In mass-flux variables `m=ρ0 u` this is exactly the anelastic constraint.
  **Plumbing item (c) done:** `manufactured_error_metric_z` — the stretched vertical
  metric (variable-dz finite volume, walls) composed with the horizontal composite
  interface (the storm nest refines horizontally only and shares the parent's stretched
  z). On `cos(2 pi x) cos(pi z/Lz)`: uniform z is clean 2nd order (ratio 4.01), moderate
  stretching is supraconvergent (~1.8-2, the standard non-uniform-FV order);
  `test_composite_stretched_vertical_metric_second_order`. **All three step-2 plumbing
  pieces (a) wall BC, (b) face bridge, (c) stretched metric are verified**.
- **Final assembly** (`solve_composite_hz`, `composite_project_massflux_hz`): the full
  unified operator for the storm's nest geometry — the horizontal composite interface at
  every z-level + a variable-dz finite-volume vertical coupling per column (matched
  stretched z, wall BCs, `m = ρ0 u`). `solve_composite_hz` is verified 2nd-order
  (`manufactured_error_hz`; s=1.0 ratio 4.10, s=1.05 4.03). `composite_project_massflux_hz`
  projects the storm's full 3-D staggered C-grid mass fluxes (parent + nest,
  `u:(nc+1,nc,nz)`, `v:(nc,nc+1,nz)`, `w:(nc,nc,nz+1)`) so `div(m)=0` across the interface
  to ~1e-13 for the nest (walls) and the parent (periodic horizontal), verified by an
  independent recomputation from the written-back arrays
  (`test_composite_hz_unified_operator_second_order`,
  `test_composite_projection_hz_full_storm_divergence_free`).
- **Call site** (`storm_dynamics.nesting.composite_project_two_level`): the composite
  two-level anelastic pressure projection wired to real parent+nest `FlowState`s — forms
  `m* = ρ0 u*` on each level, maps the nest footprint, calls `composite_project_massflux_hz`
  (physical `hx`, parent `periodic_h`), recovers `u = m/ρ0`, writes back. Verified on real
  staggered arrays with a stretched anelastic density profile: `div(ρ0 u)` → machine
  precision across the interface (`test_composite_project_two_level_call_site_divergence_free`).
  Requires a cell-aligned, matched-z nest (`NestSpec.aligned`) in a square parent. **The
  AMR pressure projection across a refinement interface is now complete end-to-end**
  (algorithm → operator → storm-array projection → call site), all verified. Enabling it as
  the default in `run_concurrent_nest` (it replaces the nest-boundary relaxation with the
  interface coupling) is opt-in; anisotropic/non-aligned nests are the only further
  generalisations. See `docs/amr_design.md`.
- `examples/composite_projection_demo.py` — runnable end-to-end demo of the call site:
  matures a small square parent, builds a cell-aligned nest, kicks both levels with a
  divergent perturbation, runs `composite_project_two_level`, and shows the anelastic
  divergence collapse from ~1e-1 to ~1e-16 everywhere including the interface. Documented
  in the README M3 section.
- `storm_dynamics/poisson_mg.py` — **geometric-multigrid Poisson** (the AMR
  projection kernel, since pyAMReX exposes no `MLMG`): 2-D cell-centred periodic
  V-cycle (red-black GS, full-weighting restriction, bilinear prolongation).
  Verified: h-independent convergence (~8 V-cycles to 1e-10 at n=64/128/256) and
  2nd-order accuracy (error ∝ h²) on a manufactured solution. The composite
  (coarse-fine) coupling on top of this kernel is the remaining AMR-projection step.
- `storm_dynamics/amr_port.py` — **AMR port scaffold**: the field on an AMReX
  `MultiFab` (framework data model + ghost exchange via `fill_boundary`) stepped by
  our flux-form NumPy physics. Verified on WSL/pyAMReX: a 32³ periodic advection has
  total-mass drift `0.0` — the "AMReX infrastructure + our RHS" binding the full
  port rests on. Import-safe without pyAMReX (test skips); runs in the WSL amr312
  env.
- `scripts/build_pyamrex_wsl.sh` — turnkey, reproducible pyAMReX build for
  WSL2/Ubuntu (the M3 "full AMR" framework path). Validated on WSL: pyAMReX `26.08`
  builds and imports (Miniforge **Python 3.12**; `JOBS=2` to avoid the pybind OOM
  on a 7 GB WSL). The default build exposes `MultiFab`/`Geometry` and the AMR
  hierarchy (`AmrCore`/`AmrMesh`) but **not** `MLMG`/`FluxRegister` — the port needs
  those bindings enabled (see `docs/amr_design.md`).
- `storm_dynamics/amr.py` — **AMR Milestone 1**: a pure-NumPy reference
  implementation of Berger–Colella **refluxing** on a static 2-level hierarchy.
  Verified: total-mass drift over 40 steps is `1.1e-4` without refluxing and
  `2.0e-16` with it (`test_amr_refluxing_conserves_across_interface`) — the flux
  register restores exact conservation across the coarse–fine interface. Ports
  directly onto a framework `FluxRegister` (AMReX/Chombo).
- `examples/supercell_tornadogenesis.py` (runnable, prints the rotation diagnostics;
  `--plots` writes the rotation figures).
- `storm_dynamics/plotting.py` — rotation figures: mid-level / near-surface ζ
  slices (the cyclonic/anticyclonic split couplet), vertical velocity, the
  environmental hodograph with SRH/shear, and the rotation time series; plus an
  **animated GIF** of the evolving vorticity (`animate_rotation`, Pillow — no
  ffmpeg) and a history-CSV writer. Example flags `--animate` / `--fps` / `--csv`
  (frame capture via `StormSimulation.run(capture_frames=True)`).
- `configs/storm_supercell.yaml`, `configs/storm_tornadogenesis.yaml` (declarative).
- `docs/storm_dynamics_guide.md` (model, what it can/cannot claim, resolution
  limits, references), `docs/MANUAL_storm_dynamics.html` (styled self-contained
  HTML manual matching `docs/MANUAL.html`), and `src/storm_dynamics/handoff.md`.
- `.gitattributes` pinning the two SHA-256-guarded `_engine/` reference scripts to
  their canonical CRLF byte form (see Fixed).

### Fixed
- Restored CRLF line endings on the checksum-guarded ice reference script
  `_engine/Nucleation_model_H2O_vapour_solid_Sim_2026_paper.py` (content
  byte-identical; only the checkout's EOL normalisation had changed, breaking the
  SHA-256 guard and `test_18_reference_preserved` in this working copy). Guarded
  by `.gitattributes` so it cannot recur. `--validate` green again.

### Preserved
- `_engine/` science untouched; `meteorological_flow` behaviour unchanged (only
  new files added). `--validate` passes; the existing suite is unchanged (aside
  from the pre-existing, environment-only `test_b3_xarray_roundtrip` xarray/scipy
  NetCDF incompatibility, which is not related to this work).

### Changed
- `pyproject.toml`: `packages.find` include += `storm_dynamics*`.

## [1.0.0] — 2026-08-20 — package reorganization

### Added
- Installable `src/`-layout package `met_water_nucleation` with a facade that
  re-exports the validated engine API.
- `pyproject.toml` with metadata, required deps (numpy/scipy/matplotlib),
  optional extras (`netcdf`, `grib`, `pandas`, `io`, `dev`) and the
  `met-water-nucleation` console entry point.
- `python -m met_water_nucleation` module entry point.
- Root `met_h2o_nucleation.py` backward-compatibility shim (DeprecationWarning).
- `tests/conftest.py` (CWD-independent import bootstrap).
- `docs/architecture.md`, `docs/migration-guide.md`, `MIGRATION_MANIFEST.md`.
- `configs/` (declarative scenario YAMLs), `scripts/` (run_validation,
  regenerate_outputs), `legacy/`, `references/`, `data/` section READMEs.
- `outputs/README.md` naming convention (`outputs/<scenario>/<run-id>/`).

### Changed
- Engine bundle (core + 2 reference scripts + `met` + `hca`) relocated as a
  unit to `src/met_water_nucleation/_engine/` (byte-identical `git mv`).
- Core's own ecosystem (its `README.md`, `MANUAL_*.html` and historical `out_*`
  outputs) relocated into the core's new folder alongside the core `.py`,
  restoring the original self-contained bundle (still untracked/gitignored).
- Examples → `examples/` (imports → package; write to `outputs/<scenario>/`).
- Tests → `tests/` (import → package; round-trip artifact → system temp dir).
- Docs → `docs/` (`README.md` → `docs/index.md`; new root `README.md`).
- `out_met_nucleation/` → `outputs/`.
- `.gitignore` updated (per-scenario output subdirs, orphaned core dir).

### Preserved
- All five engine files retain their pre-reorg SHA-256 checksums; `--validate`
  still passes (core [1]–[21], ice SHA-256 unchanged).
- 24/24 tests pass before and after.
- Committed flat reference outputs kept tracked and unchanged on disk.

### Resolved decisions (2026-08-20)
- `LICENSE`: **MIT** added at the repo root (`pyproject.toml` license → MIT +
  SPDX classifier). The integrity-guarded core/reference models remain
  read-only; the MIT licence permits modification but editing guarded files
  invalidates `--validate` (noted in `LICENSE`).
- Orphaned `unified_h2o_nucleation_climate/` directory (the core's own untracked
  docs + historical `out_*` outputs): **relocated into the core's new folder**
  `src/met_water_nucleation/_engine/unified_h2o_nucleation_climate/` (option 2
  from `legacy/README.md`), restoring the core's self-contained ecosystem
  alongside the core `.py`. Still untracked/gitignored as before; the empty old
  shell at repo root was removed.