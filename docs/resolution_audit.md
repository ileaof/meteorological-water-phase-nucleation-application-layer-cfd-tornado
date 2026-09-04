# Resolution-dependence audit

**Phase 1 of the REVIEW_REQUEST.md programme.** Every quantity in the codebase that converts a
physical length to a cell count, or that is *stored* as a cell count while meaning a length.

**Classification key**

| code | meaning |
|---|---|
| **D** | intrinsically discrete — correctly cell-based |
| **L** | represents a physical LENGTH — must be requested in metres |
| **S** | minimum numerical support — a floor on cells, not a physical scale |
| **A** | ambiguous — needs clarification |

**Rule enforced from now on** (`storm_dynamics/scales.py`): a physical scale is requested in
metres; the discretisation is *reported*, never hidden; a scale that cannot be represented is an
explicit failure, never a silent substitution; and a represented region is never *larger* than
requested (floor, not round).

---

## Audited occurrences

| # | location | expression | class | status |
|---|---|---|---|---|
| 1 | `vortex_diagnostics.surface_connection_report` | `R = max(3, int(radius_m/dx))` | **L + S conflated** | **FIXED** — now `scales.cells_for_length`; under-resolved ⇒ NaN + `valid=False`, never a substituted radius |
| 2 | `vortex_diagnostics.surface_connection_report` | `nb = max(2, int(border_frac*nx))` | **A** | **IMPROVED** — `border_m` (metres) added and preferred; `border_frac` retained for compatibility but is a *domain fraction*, not a fixed length, and is now labelled as such in the returned `interior_margin` block |
| 3 | `vortex_diagnostics.find_vortex_center` | `b = max(1, int(border_frac*nx))` | **A** | open — domain-fraction margin; mesh-independent only at fixed domain size |
| 4 | `vortex_diagnostics.pressure_deficit` | `nb = max(1, int(ambient_frac*shape))` | **A** | open — ambient ring as a domain fraction |
| 5 | `vortex_diagnostics.circulation` | `r <= radius_m` on physical coords | **L** | correct — already physical; but a radius < dx encloses ~1 cell (support not checked) |
| 6 | `vortex_diagnostics.tangential_profile` | `r_max_m` | **L** | correct — physical |
| 7 | `nesting.NestSpec.relax_width` | `4` cells | **L stored as cells** | **FIXED (opt-in)** — `relax_width_m` pins a physical width via `effective_relax_width`; default unchanged for reproducibility |
| 8 | `nesting.effective_relax_width` | `max(2, int(round(relax_width_m/dx)))` | **L + S** | correct — reports cells for a requested metre width; `min 2` is support |
| 9 | `nesting.interior_near_surface_zeta` | `margin = relax_width + 2` | **S on top of L** | **FIXED** — now uses `effective_relax_width` |
| 10 | `nesting.build_nest_grid` | `max(4, int(round(Lx/dxp*refine)))` | **D** | correct — a grid dimension is intrinsically discrete |
| 11 | `nesting.cluster_to_box` | `margin` in cells | **A** | open — tagging margin; cells are arguably right (it is a stencil), but undocumented |
| 12 | `nesting.tag_cells` | `border` in cells | **S** | acceptable — an exclusion band tied to the sponge, which is itself now physical |
| 13 | `nesting.regrid_spec` | `min_cells=6` | **S** | correct — support floor for a nest footprint |
| 14 | `limited_area.lateral_relaxation_weight` | `width=8` cells | **L stored as cells** | **OPEN** — same defect class as (7), on the *parent* Davies zone. Not yet converted; parent dx is fixed within a run so it has not yet corrupted a comparison, but it will the moment a parent-resolution study is run |
| 15 | `atmospheric_data/driver._grid_dims` | `max(16, int(round(Lz/dx)))`, `min(n, max_n)` | **D + S** | correct — grid dimensions, with an explicit cap |
| 16 | `atmospheric_data/driver.run_multilevel_real_case` | `max(2, int(round(dx_parent/dx_nest)))` | **D** | correct — refinement ratio is intrinsically integer |
| 17 | `microphysics_coupling.sediment` | `nsub = max(1, ceil(v*dt/dz))` | **D** | correct — sub-cycle count from a CFL condition |
| 18 | `precip_microphysics/sedimentation` | `n = max(1, ceil(vt*dt/dz))` | **D** | correct — same |
| 19 | `storm_dynamics/config.ForcingConfig.radius_m` | metres | **L** | correct — already physical |
| 20 | `forcing.py` | `r / fc.radius_m` on physical coords | **L** | correct |

---

## Measured impact of the defect in row 1

A nominal **400 m** V_rot sampling radius, before and after:

| dx [m] | old `max(3,int(400/dx))·dx` | new represented | new status |
|---|---|---|---|
| 600.0 | **1800.0** (4.5× enlargement) | 0.0 | **under_resolved** |
| 300.0 | **900.0** (2.25×) | 300.0 | **under_resolved** (1 cell < 3) |
| 66.7 | 333.5 | 333.5 | ok (−16.6 %) |
| 22.2 | 399.6 | 399.6 | ok (−0.1 %) |

Two consequences for past results:

* A **600 m vs 300 m** comparison at a nominal 400 m radius would have measured **1800 m against
  900 m** — a 2× mismatch biased toward the coarse mesh. This is now refused outright.
* The **matched-domain** nest runs (dx 66.7 and 22.2 m) compared **333.5 m against 399.6 m**, a
  16.5 % mismatch, not the "identical 400 m radius" previously claimed. Those runs were already
  unusable for boundary-contamination reasons, so no verdict changes — but the claim was wrong.

## Remaining open items

* **Row 14** (`limited_area` Davies width in cells) is the same defect class as row 7 and is not
  yet converted. It has not yet corrupted a result because parent dx is constant within a run.
* **Rows 3, 4, 11** are domain-fraction or stencil margins that are mesh-independent only at
  fixed domain size. They need a decision: physical metres, or an explicit "fraction of domain"
  contract documented at the call site.
* **Row 5**: `circulation` is physically correct but does not check numerical support; a radius
  below a few cells yields a meaningless integral without saying so.
