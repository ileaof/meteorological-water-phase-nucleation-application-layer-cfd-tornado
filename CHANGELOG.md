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
- `examples/supercell_tornadogenesis.py` (runnable, prints the rotation diagnostics;
  `--plots` writes the rotation figures).
- `storm_dynamics/plotting.py` — rotation figures: mid-level / near-surface ζ
  slices (the cyclonic/anticyclonic split couplet), vertical velocity, the
  environmental hodograph with SRH/shear, and the rotation time series.
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