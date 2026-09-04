> # ⚠ READ FIRST — validity notice (2026-09-04)
>
> **Every percentage in this document quoted against "the observed 26 m/s" is INVALID**, on both
> sides of the ratio:
>
> * **Denominator** — 26 m/s was the radar's *Nyquist ceiling*, and the couplet was extracted
>   21 km from Moore. The corrected target is **V_rot 39.5 m/s [34, 45] at 208 m through a 327 m
>   beam at 20.25 km** (see §2).
> * **Numerator** — the model V_rot values were produced by diagnostics since found defective:
>   a sampling radius that silently rescaled with the mesh, a moving nest steered by its own
>   sponge vorticity, and a tangential profile whose fixed 24 bins made `v_theta` at `r=0` the
>   raw meridional wind of a single cell. Peak values of 12.28 and 8.47 m/s are **retracted**,
>   and both `outputs/uniform_resolution` tornado classifications are **retracted**.
>
> What survives is the **mechanistic** content — internally consistent idealized-vs-idealized
> comparisons: the tilting-geometry bottleneck, free evolution building streamwise geometry, the
> surface-layer closure *form* being the switch, and near-surface vertical resolution dominating
> horizontal. Those do not depend on the observation.
>
> Model/observation comparison must go through `atmospheric_data/radar_operator.py`. See
> `docs/REVIEW_REQUEST.md` for the full defect synthesis.

# Tornadogenesis: findings, honest limits, and the sustainment lever

*Moore, OK — EF5 of 20 May 2013 — as the real-data test case.*
*Last updated 2026-09-02. Companion to [`ROADMAP.md`](ROADMAP.md) and [`REAL_CASE_DATA.md`](REAL_CASE_DATA.md).*

This note records what we learned trying to reproduce a **real tornado from real data**. The
short version: the package reproduces essentially all of severe-storm phenomenology, but the
**tornado-scale low-level vortex is the one thing it does not yet capture** — and we now know
*why*, quantitatively, and which lever moves it. Honest-failure reporting is a project value;
these are negative results as much as positive ones.

---

## 1. Why the tornado is *the* hard problem

Everything else we simulate — density currents (Straka), gravity waves, convection, the storm
cell, even the **mid-level mesocyclone** — is robust and near single-scale, and a well-built
anelastic LES core reproduces it. The tornado is uniquely hard because it concentrates **every**
difficulty at once:

1. **Scale separation of ~5 orders of magnitude, simultaneously.** Supercell ~50–100 km ·
   mesocyclone ~3–5 km · tornado ~0.3–1 km · suction/corner-flow vortices ~10–50 m. Resolving the
   tornado needs Δx ~10–30 m, but the domain must still hold the whole supercell that *generates
   and sustains* it. That is 10⁷–10⁹ cells — the regime of the leading LES groups (Orf's El Reno
   run at ~30 m on a supercomputer). Even our finest nest at **46 m under-resolves the observed
   ~0.3 km TVS couplet**.
2. **The vortex is the final link of a fragile causal chain**, not a forced feature: supercell →
   mid-level mesocyclone (tilting of *environmental* horizontal vorticity) → **low-level
   mesocyclone**. The low level comes from a *different* source — **baroclinic** vorticity
   generation along the storm's own cold-pool/forward-flank gradient, then tilted and stretched by
   the low-level updraft. This hinges on the **microphysics** (rain evaporation → cold-pool
   strength): too cold → outflow-dominated, no tornado; too warm → no baroclinic vorticity. A
   narrow window.
3. **It is a threshold / bifurcation phenomenon.** In nature, *most* supercells produce no
   tornado; a small change flips the outcome. Reproducing it deterministically from imperfect
   initial data is at the research frontier (why the VORTEX field campaigns exist).
4. **The input data are too coarse to hold the seed.** ERA5 (~31 km) smooths the low-level shear:
   we measured **SRH ≈ 146–152 m² s⁻²** from ERA5 versus **~247** from the real KOUN radiosonde.
   The environment we can *download* has already lost part of the low-level helicity that makes the
   tornado — a limitation of the reanalysis, shared by every forecaster.

**One-line assessment:** the tornado is the only phenomenon that demands LES resolution, a fragile
multi-stage causal chain, accurate microphysics, *and* a correctly-hit threshold — all at once.

---

## 2. The observed target (real NEXRAD KTLX Level II)

Read on WSL2 with Py-ART + nexradaws (no credentials) — `deploy/wsl2_nexrad_moore.py`,
`deploy/wsl2_extract_velocity.py`. Scan **KTLX 2013-05-20 20:20:58 UTC, 0.5° sweep**:

> ## ⚠ RETRACTED AND CORRECTED (2026-09-04)
>
> **Every number in the original version of this table was wrong, and every percentage quoted
> against it elsewhere in this document is therefore invalid.** Two independent defects:
>
> 1. **`V_rot = 26 m/s` was the radar's NYQUIST VELOCITY**, not a measurement. The sweep's own
>    metadata gives `nyquist = 26.12 m/s`; the field saturated at exactly ±26.000 m/s in *every*
>    sub-region, and 9 adjacent-gate pairs differed by >40 m/s (the 2×Nyquist folding signature).
>    `(max − min)/2` over any region containing a folded gate returns 26 **by construction**.
> 2. **The couplet was extracted 21 km from Moore.** The old script took the strongest
>    inbound/outbound pair over the whole cropped field; that pair sits 4.9 km from KTLX (near-radar
>    clutter) while Moore is at 19.3 km.
>
> Corrected by dealiasing the Level II velocity (Py-ART region-based) and constraining the couplet
> to the Moore mesocyclone — `deploy/wsl2_dealias_moore.py`, stable across 3/5/8 km search radii:

| Quantity | OLD (retracted) | **CORRECTED** |
|---|---|---|
| Reflectivity core | ~70 dBZ | ~70 dBZ (unchanged) |
| Velocity couplet Δv | 52 m s⁻¹ (= 2×Nyquist) | **79.0 m s⁻¹** |
| Rotational velocity V_rot | 26 m s⁻¹ (= Nyquist) | **39.5 m s⁻¹, interval [34, 45]** |
| Couplet separation | ~253 m (wrong feature) | **584 m** (3 rays × 1 gate) |
| Range from KTLX | ~30–41 km implied | **20.25 km** |
| Beam diameter | 250 m assumed | **327 m** |
| Sample height | ~460 m "AGL" | **208 m above the antenna** |
| Elevation | "0.5°" | **0.5211°** |
| Angular velocity Ω = Δv/sep | mislabelled 0.205 s⁻¹ | **0.135 s⁻¹** |
| Vertical vorticity ζ = 2Ω | — | **0.271 s⁻¹** |

> **Unit error, additionally:** the old `0.205 s⁻¹` row was labelled "azimuthal shear ≈ vertical
> vorticity" and compared against the model's ζ. For solid-body rotation `Δv/sep = Ω = ζ/2`, so
> that row compared unlike quantities and was off by a **factor of two** on top of being derived
> from the wrong couplet.
>
> **Uncertainty:** the dominant term is the *estimator choice* (top-1 gate 39.49 vs 3×3 median
> 33.49), **not** the dealiasing — a fold round-trip (`deploy/wsl2_fold_roundtrip.py`) recovers a
> known V_rot exactly up to ~45 m/s and only fails at 60, so the interval stays asymmetric on the
> high side.
>
> **How to compare a model against this:** never against a grid-point wind. Push the model through
> `atmospheric_data/radar_operator.py` at `RadarSpec(elevation_deg=0.5211)`, range 20.25 km, and
> compare observable to observable. The model-mesh penalty vanishes by dx ≈ 100 m
> (`scratchpad/mesh_recovery_curve.py`); the residual 0.700 recovery is the irreducible beam
> penalty.

Figures: `docs/media/storm/ktlx_moore_2013_radar.png`, `docs/media/storm/moore_sim_vs_radar.png`.

---

## 3. What the simulation actually produced (honest, quantitative)

Real ERA5 environment for 2013-05-20 20 UTC (CAPE ~1900–2050 J kg⁻¹, CIN −120 to −130,
0–6 km shear 27 m s⁻¹, SRH 146–152) as the base state; an idealised trigger; the GPU AMR cascade.

### Attempt A — 3-level AMR cascade to 46 m (`scratchpad/moore_real_funnel_gpu.py`)

| | Observed (KTLX) | Simulated |
|---|---|---|
| Low-level (≈460 m) V_rot | 26 m s⁻¹ | **3.2 m s⁻¹** |
| Low-level Δv | 52 m s⁻¹ | 6.3 m s⁻¹ |
| Low-level vorticity | 0.205 s⁻¹ | **0.026 s⁻¹** |
| Peak vorticity in the column | (TVS, near surface) | 0.093 s⁻¹ **at z ≈ 14.6 km** |

**Correction of an earlier overstatement.** An earlier note celebrated "ζ ≈ 9.3×10⁻² s⁻¹ at 46 m —
the funnel scale." That peak is at the **~14.6 km anvil/storm top**, *not* a surface funnel. The
tornado-relevant **low-level** rotation is only **0.026 s⁻¹ (~13 % of observed)**, with the
interior signal partly a nest-boundary artifact. The README and metrics were corrected
(`outputs/nexrad_moore/sim_vs_radar_metrics.json`).

### Attempt B — sustained-maturation + storm-following nest (`scratchpad/moore_sustained_gpu.py`)

Larger domain (100 km), 28 min maturation, then a storm-following 521 m nest. Result: the
**single warm bubble decays** — w_max **9.5 → 1.1 m s⁻¹** over the 28 min — and the nest, applied
to an already-dying storm, measured low-level ζ **9.5×10⁻⁴ s⁻¹, V_rot 2.3 m s⁻¹**. Worse than
Attempt A, precisely *because the storm died before nesting*.

**This is the central diagnosis:** with a *real* environment the failure is not (only) nest
resolution — **the parent storm does not sustain**. Two coupled causes:

- **The real CIN cap (−120 to −130 J kg⁻¹) suppresses re-development.** After the bubble's one
  impulse, nothing lifts new parcels through the cap, so the updraft detrains and decays. The
  classic Weisman–Klemp supercell sustains because its idealised environment has ~zero CIN; the
  real cap is the killer for a one-shot trigger.
- **The parent Δx (1.5 km) under-resolves the updraft** (~2–4 km wide → only 1–3 cells), so it is
  over-diffused and the updraft–shear feedback that keeps a supercell alive is weak.

The real Moore storm broke its cap by **sustained mesoscale ascent** — dryline convergence plus
afternoon boundary-layer heating — which a single thermal cannot represent.

---

## 4. The lever: sustained mesoscale-ascent forcing

Implemented as an **additive, opt-in** feature (`src/storm_dynamics/forcing.py`,
`MesoForcingConfig` in `config.py`, wired into `core._predictor`; off by default so idealised runs
and every existing test are untouched — `tests/test_forcing.py`, 4 tests).

A smooth low-level **heating (+ moistening) cylinder** held for `duration_s` continuously lifts
parcels through the CIN cap — a dryline/convergence proxy — so a supercell can *establish* from a
real capped environment; afterwards the forcing is removed and the storm must sustain on its own
dynamics. Knobs: `heat_rate_K_s`, `moist_rate_kgkg_s`, `radius_m`, `z_top_m`, `duration_s`,
`center`. Attempt C (`scratchpad/moore_forced_gpu.py`) pairs it with a finer parent (Δx 750 m).

### Attempt C — sustained forcing + finer parent (the lever, tested)

Forcing on for the first 25 min (heat 0.006 K s⁻¹, moisten 3×10⁻⁶ kg kg⁻¹ s⁻¹, r 7 km, below
2.5 km), Δx 750 m parent, then removed; a storm-following **250 m** nest over the updraft.

**The sustainment worked — this is the qualitative win.** The updraft
climbed to w_max **16.8 m s⁻¹** (t≈500 s), settled to a sustained ~7 m s⁻¹, and — the decisive
test — **re-intensified on its own *after* the forcing switched off** (w_max 7.2 → 8.7 m s⁻¹ over
t=1300→1700 s, forcing off at 1500 s). It did **not** die as the single bubble did (Attempt B:
9.5 → 1.1). The parent low-level vorticity climbed **monotonically** the entire run
(6.2×10⁻⁴ → 6.2×10⁻³ s⁻¹), and the mid-level mesocyclone reached 1.0×10⁻² s⁻¹ — a real,
organising low-level mesocyclone rather than a decaying pulse or a boundary artifact.

| Attempt | trigger | nest Δx | low-level **V_rot** | low-level ζ | storm fate |
|---|---|---|---|---|---|
| A | single bubble | 46 m | 3.2 m s⁻¹ | 0.026 s⁻¹ | pulses |
| B | single bubble | 521 m | 2.3 m s⁻¹ | 9.5×10⁻⁴ s⁻¹ | **decays** (w 9.5→1.1) |
| C | sustained forcing | 250 m | 4.9 m s⁻¹ | 6.6×10⁻³ s⁻¹ | **sustains** (self-holds after forcing off) |
| **D** | **forcing + cascade** | **28 m** | **6.0 m s⁻¹** | 1.17×10⁻² s⁻¹ | sustains; surface-connected meso |
| *observed* | *real* | ~250 m (0.5° beam) | **26 m s⁻¹** | 0.205 s⁻¹ | *tornado* |

*(ζ is Δx-dependent — A's 0.026 at 46 m is not comparable across rows; **V_rot**, a velocity, is
the resolution-fair metric. C and D come from a storm that **lives without the crutch**; D (below)
is the strongest and is surface-connected.)*

**Honest read:** V_rot ~4.9 m s⁻¹ is the best low-level rotation we have obtained from real data,
but it is still only **~19 % of the observed 26 m s⁻¹** — *not* a tornado. The lever is validated
(sustained supercell + organising low-level meso), but the remaining gap is now squarely the
things in §5: the 250 m nest still far under-resolves the 0.3 km TVS, the 600 s nest window is
short for full stretching intensification, and ERA5's SRH (152) is well below the real KOUN (247).
The path forward is to **cascade the sustained 250 m nest down to < 50 m at the storm base with a
longer window** — now that we finally have a living parent to nest into.

### Attempt D — deep cascade to 28 m, storm-relative (the refinement study)

`scratchpad/moore_cascade_gpu.py`. The sustained forced supercell, run **storm-relative**
(Bunkers motion C = (14.2, 3.8) m s⁻¹ subtracted from the base wind — a Galilean shift that leaves
tilting/stretching invariant but holds the storm quasi-stationary so the fine nests keep it), then
refined **750 → 251 → 84 → 28 m** (three ×3 refinements) with the concurrent multi-level driver
(`run_multilevel_nest`: time-sub-cycled coarse→fine boundaries + conservative restriction back up,
`restrict_momentum=True`), centred on the low-level mesocyclone, 180 s window. ~86 min on the GPU.

**Correct structure, at last.** The finest level shows a **deep, surface-connected column** of
rotation — peak |ζ| at ~500 m–1 km and still ~half that at 50 m — i.e. a genuine *low-level*
mesocyclone reaching the ground, not the anvil-top artifact of Attempt A:

```
 z ≈ 50 m   ζ = 6.0e-3   ####################
 z ≈ 490 m  ζ = 1.17e-2  ########################################   <- low-level meso peak
 z ≈ 1020 m ζ = 1.09e-2  #####################################
 z ≈ 1670 m ζ = 9.3e-3   ###############################
 z ≈ 2460 m ζ = 1.0e-2   ##################################
 z ≈ 3410 m ζ = 7.0e-3   ########################
 z ≳ 4.5 km              (drops off — the column is a low/mid-level feature)
```

**Low-level V_rot = 6.0 m s⁻¹** (ζ 1.17×10⁻² at 28 m) — the best yet, but still **~23 % of the
observed 26 m s⁻¹**.

**The key numerical finding — resolution is *not* the dominant limiter.** Refining the nest 9×
(250 m → 28 m) raised V_rot only **4.9 → 6.0 m s⁻¹ (+22 %)**. If grid resolution were the main gap
to the observed tornado, a 9× refinement would have moved it far more. The diminishing return says
the ceiling is the **source circulation** the stretching has to work on, which is set by the
*environment and the storm's own low-level vorticity budget*, not the mesh:

- **ERA5 SRH ≈ 152 vs the real KOUN ≈ 247** — less streamwise horizontal vorticity available to
  tilt into the vertical.
- **Modest sustained updraft** (w ~ 6 m s⁻¹, peaking ~21) — less vertical stretching of that
  vorticity.
- **The baroclinic (cold-pool) low-level vorticity source** — governed by rain-evaporation
  microphysics — is likely under-represented; without a strong forward-flank baroclinic zone the
  low-level vortex has a weak parent circulation regardless of Δx.

So the honest ranking of levers flips: **environment SRH + cold-pool microphysics now rank above
resolution.** We have built a correctly-structured, surface-connected low-level mesocyclone from
real data — a real qualitative milestone — but reaching tornado intensity needs a *stronger source
circulation*, not merely a finer grid.

### Attempt E — the real KOUN sounding (SRH 254 vs 152): the SRH lever, tested and *falsified*

`scratchpad/moore_koun_cascade_gpu.py`. Attempt D's ranking put *real low-level SRH* as the top
lever. So we replaced the ERA5 environment with the actual **KOUN 2013-05-21 00Z radiosonde**
(CAPE 1353, CIN −247, shear 32, **SRH 254**) — vs ERA5's (CAPE 1885, CIN −130, shear 27, SRH 152)
— and re-ran the identical storm-relative cascade to 28 m (the stronger −247 cap handled by a
slightly stronger forcing).

The KOUN environment did give a **stronger, steadier parent updraft** (w_max held ~13–14 m s⁻¹ vs
ERA5's ~6). But the **low-level rotation did not increase**. Measured consistently (max low-level
V_rot over z < 1.5 km, same window for both):

| Env | SRH | parent updraft | low-level **V_rot** (fair) | low-level peak |ζ| |
|---|---|---|---|---|
| ERA5 (D) | 152 | w ~6 | **6.7 m s⁻¹** | 1.4×10⁻² |
| **KOUN (E)** | **254** | **w ~13–14** | **6.0 m s⁻¹** | 2.2×10⁻² (more surface-peaked) |

*(The run's headline "V_rot 2.1" was a single-level (500 m) sampling artifact — 500 m is a local
dip in KOUN's profile between a surface max and a 1 km peak. Re-measured apples-to-apples, both
environments give ~6 m s⁻¹.)*

**Result: SRH is NOT the bottleneck.** Raising the real environmental SRH by **+67 % (152 → 254)
did not raise the low-level rotational velocity** — both sit at ~6 m s⁻¹, ~23 % of the observed 26.
KOUN's vorticity profile is more surface-concentrated (peak 2.2×10⁻² at 1 km, 1.5×10⁻² at 50 m),
so the helicity does sharpen the *structure* — but the tornado-relevant *velocity* is unchanged.
This **falsifies** Attempt D's top-lever hypothesis and points the finger squarely at the two
things both runs share.

### Attempt F — resolve the updraft on a fine 250 m parent (Alternative 1)

`scratchpad/moore_fineparent_gpu.py`. Attempts C–E kept the storm-generating **parent** coarse
(750 m–1.5 km) and only refined the *nest* at the end, so the parent updraft stayed weak (w ~6–14).
The old idealised model reached w ~25–56 m s⁻¹ by *resolving the updraft* (convergence study M8;
Bryan 2003: Δx ≲ 250 m). So we matured the real KOUN storm on a **fine 250 m parent**
(120×120×48, 6.9×10⁵ cells) from the start, then cascaded to 28 m. *(This required a real fix:
the direct `PressureSolver` — built eagerly, a host `splu` of an N×N Laplacian — hung >300 s at
this size and is **unused** on low-memory grids; it is now skipped, and the parent builds in 2 s.)*

**Resolving the updraft worked — and it moved the storm-scale rotation.** The parent updraft
climbed to **w_max 24.7 m s⁻¹** (~2× the coarse runs), and the **parent's own low-level V_rot
reached 8.1 m s⁻¹** — the strongest low-level rotation in the whole study (+21 % over D's 6.7).

**But the one-way nest cascade lost it.** Put through the 250 → 83 → 28 m cascade, the finest level
came back down to **V_rot 6.4 m s⁻¹** — the same ~6 ceiling as D (6.7) and E (6.0):

| Level | resolution | low-level V_rot |
|---|---|---|
| **Fine parent (one grid)** | 250 m | **8.1 m s⁻¹** |
| ↓ one-way cascade → finest | 28 m | 6.4 m s⁻¹ |

**The clue.** When the storm is resolved as *one continuous grid* (the 250 m parent) it reaches
8.1; the moment it is split into parent + nest with **one-way (downscale-only) coupling**, the nest
**re-equilibrates to its boundary forcing** and resets to ~6.4. So part of the 28 m ceiling is not
physics — it is a **numerical artifact of the artificial scale separation**: the fine nest does not
inherit and amplify the parent's vortex, because the fine→coarse (upscale) feedback is severed.
This is a distinct, testable limiter from resolution and SRH: **two-way scale coupling.**

### Attempt G — two-way coupling A/B (the lever, confirmed)

`scratchpad/moore_twoway_ab.py`. From one matured fine (250 m) parent, the **same** single-level
nest (refine 3 → 83 m) run **one-way vs two-way** (`run_concurrent_nest` `two_way=True`, the
fine→coarse injection). Two-way raised the nest low-level **V_rot 8.1 → 11.6 m s⁻¹** (and closer to
the surface, 122 m vs 208 m) and more than **doubled the parent's, 6.1 → 15.0** — the vortex↔updraft
loop intensifying the whole storm, breaking the ~6 m s⁻¹ ceiling. **Confirms** the Attempt-F
hypothesis: one-way nesting was a numerical cap. (Two-way is now wired into `run_multilevel_nest`.)

### Attempt H — KOUN two-way *deep* cascade + the vorticity budget (the measured verdict)

`scratchpad/moore_twoway_deep_gpu.py`. Real KOUN env → forced supercell → storm-relative deep
cascade to 28 m with `two_way=True`, then the **vorticity budget** read on the result. Two findings:

- **Two-way in the *deep* multi-level cascade did not reproduce G's boost** (low-level V_rot 4.9,
  class LOW_LEVEL_MESOCYCLONE; the run was also ~4× faster). Applying the aggressive rate-0.5
  injection at *every* level pair of a 4-level cascade appears to over-smooth the fine vortex,
  unlike the single coarse–fine pair in G — a controlled deep A/B (matched one-way vs two-way, or a
  gentler rate / top-interface-only) is the clean next test.
- **The budget names the real bottleneck.** The cold pool is strong (θv′ −3.5 K, C 14.6 m s⁻¹, 38 %
  area), and the **hydrostatic baroclinic generation of *horizontal* vorticity is 1.47×10⁻⁴ s⁻² —
  22× the tilting rate (6.7×10⁻⁶) and stretching (6.4×10⁻⁶)** at low levels. So the cold-pool
  source is **abundant**; only ~4.5 % of it is **tilted** into the vertical. **The limiter is the
  tilting efficiency — the geometric alignment of the baroclinic horizontal vorticity with the
  low-level updraft gradient — not the cold-pool source.** This implicates *storm structure*
  (rear-flank downdraft / occlusion, the idealised single-updraft trigger) over microphysics.
  *(Caveat: the direct vertical-vorticity baroclinic term B_z reads ~0 — it vanishes under
  hydrostatic balance and the low-memory solver does not persist `state.p` on nests; the horizontal
  baroclinic diagnostic above, computed from ρ, is the correct source measure —
  `vorticity_budget.baroclinic_horizontal_generation`.)*

### The decisive diagnosis — it is GEOMETRIC (misalignment), not a missing ingredient

Factorising the tilting term (`vorticity_budget.tilting_efficiency`), since
`tilting = ω_h · ∇_h w = |ω_h| |∇_h w| cos θ`, on the Attempt-H fields:

| z | ω_h (available) | \|∇_h w\| (tilting agent) | **alignment cos θ** | tilting |
|---|---|---|---|---|
| 51 m | 1.1×10⁻² | 1.3×10⁻³ | **−0.146** | 6.2×10⁻⁶ |
| 264 m | 7.8×10⁻³ | 3.3×10⁻³ | **+0.042** | 8.1×10⁻⁶ |
| 499 m | 5.3×10⁻³ | 2.2×10⁻³ | **+0.086** | 6.7×10⁻⁶ |
| 1043 m | 5.7×10⁻³ | 2.5×10⁻³ | **+0.275** | 6.3×10⁻⁶ |

**Neither ingredient is missing.** Horizontal vorticity is abundant (fed by the cold pool at
1.6×10⁻⁴ s⁻²) and a low-level updraft gradient exists — they are **geometrically misaligned**
(cos θ ≈ 0.04–0.09), and at the surface they are **anti-aligned (−0.146)**, so tilting there
generates *anticyclonic* vorticity. Alignment **improves with height** (+0.275 at 1 km) — precisely
why this model produces a mid-level mesocyclone but no surface vortex. The low-level vorticity is
also only **~49 % streamwise** (0.32 at 1 km), where tornadic supercells need it strongly streamwise.

**This single result explains every earlier negative:** refining the mesh (D) cannot fix an
alignment; more ambient SRH (E) adds vorticity that is still misaligned; a stronger updraft (F)
tilts a misaligned field no better; and two-way coupling (G) helped precisely because the upscale
feedback let the fine vortex *reorganise the local flow*. **The missing ingredient is the storm
structure — the rear-flank downdraft / occlusion process that reorients baroclinic vortex lines into
a streamwise, updraft-aligned geometry** — not resolution, not the environment, not microphysics.

### Attempt I — the falsifiable test: a freely-evolving supercell DOES reorient the vorticity

`scratchpad/supercell_alignment_evolution.py`. The misalignment diagnosis makes a prediction: a
storm allowed to evolve *freely* (idealised strong environment — CAPE 2225, shear 41, SRH 648;
warm-bubble trigger, **no held forcing**; 90 min at dx 600 m, storm-relative) should, as it occludes,
**raise the low-level alignment and streamwise fraction** — before V_rot. It does
(`docs/media/storm/supercell_alignment_evolution.png`):

| time | 5 min | 20 | 30 | 35 | 45 | 55 | 70 | 90 |
|---|---|---|---|---|---|---|---|---|
| **alignment cos θ** | +0.02 | +0.09 | +0.15 | **+0.20** | +0.14 | −0.04 | +0.03 | +0.07 |
| **streamwise frac** | 0.40 | 0.36 | 0.39 | 0.42 | 0.47 | 0.51 | 0.58 | **0.64** |
| **low-level \|ζ\| (10⁻³)** | 2.8 | 3.0 | 4.2 | 5.6 | **9.5** | 8.5 | 6.7 | 2.5 |

The low-level **alignment rose 10× to a peak of 0.20** (vs the forced runs stuck at 0.09) during the
occlusion phase, the **streamwise fraction rose monotonically 0.40 → 0.64** (a tornadic-supercell
signature), and the **low-level vorticity peaked at 9.5×10⁻³ — 3.4× its start — lagging the
alignment peak** (tilting first, then the vortex builds and stretching amplifies it), before the
mesocyclone cycled (normal cyclic mesocyclogenesis). This is at the bare 600 m parent, **no nest**.

**Verdict: affirmative.** The model *can* build the streamwise, updraft-aligned low-level geometry —
it just needs a **freely-evolving RFD/occlusion**, which the forced/held/short runs (A–H) never
allowed. (*The run script's automated one-line verdict compared the final alignment to the
mid-point and so read the cyclic downswing as "no rise" — a bad metric, corrected to the peak-rise +
streamwise-trend measure. The honest signal is the peak and the monotone streamwise rise above.*)
The remaining gap is purely **resolution at the occlusion**: a storm-following nest dropped into
*this* freely-evolved storm at its low-level-meso peak (~45–55 min) is the path to the tornado-scale
vortex — now with the geometry already correct.

### Attempt J — resolving the vortex to 22 m: a low-level mesocyclone, but ELEVATED

`scratchpad/tornado_occlusion_gpu.py`. The freely-evolved supercell (Attempt I) matured to its
low-level-meso peak (t=2800 s), then a deep one-way cascade **600 → 200 → 67 → 22 m** centred on the
occluding mesocyclone. The clean run (1.9 km finest domain so the vortex stays interior — a first,
1.2 km attempt let the vortex drift to the nest edge) gives the vertical V_rot profile
(`docs/media/storm/tornado_occlusion_profile.png`):

| z | 40 m | 208 m | 492 m | 942 m | 1488 m |
|---|---|---|---|---|---|
| **V_rot** | **2.7** | 4.5 | 4.9 | 4.8 | **7.7** |
| \|ζ\| (10⁻²) | 0.22 | 2.2 | 1.2 | 1.7 | 2.5 |

**The vortex is elevated.** V_rot *increases* with height (2.7 → 7.7 from the surface to 1.5 km); a
genuine, resolved **low-level mesocyclone** exists (V_rot ~4.9, ζ 2.2×10⁻² at 200–500 m, alignment
+0.13) but it **does not connect to the ground** — the surface V_rot is only 2.7 and the surface
vortex report is Vθ 0.9, Δp ≈ 0 (class **LOW_LEVEL_MESOCYCLONE**). Resolving to 22 m did not descend
the vortex.

---

## 5. Honest bottom line — the full arc

Across A–J, from **real** and **idealised** data, the model reproduces the tornadogenesis chain and,
crucially, we now understand each step **quantitatively**:

1. **Deep convection, splitting supercell, mid-level mesocyclone** — robustly, responding to CAPE and
   shear (the experiment matrix: remove shear and the mesocyclone collapses ~17×).
2. **Sustained supercell from real data** — the sustained-ascent forcing (§4) breaks the real CIN cap
   where a single bubble decays.
3. **A surface-connected low-level mesocyclone** — the storm-relative AMR cascade builds it.
4. **The ~6 m/s low-level-rotation ceiling is understood by elimination + measurement:** not
   resolution (D, +22 %), not environmental SRH (E, ~0 %), not updraft strength (F). Two-way
   coupling (G) breaks it (8.1→11.6) because it restores the upscale feedback. The vorticity budget
   (H) then showed the cold-pool baroclinic **source is abundant** (22× tilting) but the **tilting
   is inefficient** — the low-level horizontal vorticity is **geometrically misaligned** with the
   updraft gradient (cos θ ≈ 0.05, anti-aligned at the surface).
5. **A freely-evolving storm fixes the geometry** (I): as it occludes, alignment rises 10× (to 0.20)
   and the streamwise fraction rises monotonically (0.40 → 0.64), building the low-level ζ 3.4×.
6. **But the resolved vortex stays ELEVATED** (J): at 22 m the low-level mesocyclone is real and
   well-aligned, yet V_rot *grows* with height (2.7 m/s surface → 7.7 at 1.5 km) — **it does not
   collapse to the surface into a tornado.**

**The single remaining, well-localised gap is the surface connection** — the corner-flow collapse in
which the occlusion downdraft transports the mesocyclone's angular momentum to the ground and
surface drag drives the convergent intensification (Rotunno–Klemp; the tornado "corner flow"). That
is a sub-100 m, near-surface process; at these grid spacings and integration times it does not
spin up. This is exactly where the science frontier sits, and the diagnosis is now specific and
measured rather than a guess.

**Remaining levers, evidence-ranked:** (1) the **surface / corner-flow layer** (measured below);
(2) a **storm-following fine nest** at the occlusion — now implemented
(`nesting.follow_spec` + `run_multilevel_nest(follow_interval=…)`, same-size box re-centred on the
tracked rotation along a filtered trajectory); (3) longer freely-evolving integration for a stronger
occlusion cycle. Horizontal mesh resolution, environmental SRH, updraft strength, and the cold-pool
*source* have each been **ruled out by measurement.**

### Surface sensitivity — near-surface RESOLUTION dominates, drag does not

`examples/surface_sensitivity.py` (freely-evolving supercell; `surface_connection_report`'s
surface/aloft V_rot ratio; > 0.8 ⇒ surface-connected):

| case | C_d (applied) | first cell dz₁ | V_sfc | **sfc/aloft** | connected |
|---|---|---|---|---|---|
| baseline | 0.0120 | 49.6 m | 0.81 | 0.08 | ✗ |
| no drag | 0 | 49.6 m | 2.75 | 0.28 | ✗ |
| rough (2× C_d) | 0.0240 | 49.6 m | 0.82 | 0.08 | ✗ |
| smooth (⅓ C_d) | 0.0040 | 49.6 m | 0.91 | 0.09 | ✗ |
| coarse + log-law | 0.0042 | 49.6 m | 0.91 | 0.09 | ✗ |
| fine + bulk drag | 0.0120 | **6.2 m** | 1.25 | 0.53 | ✗ |
| fine + **log-law** drag | 0.0094 | **6.2 m** | 1.28 | 0.51 | ✗ |
| **fine + no drag** | 0 | **6.2 m** | **3.25** | **0.88** | **✓** |
| drag + heat/moisture fluxes | 0.0120 | 49.6 m | 0.81 | 0.08 | ✗ |

**This is the first SURFACE-CONNECTED vortex in the whole study** — and it isolates the last gap to
two ingredients, one expected and one not:

- **Near-surface vertical resolution is necessary.** Taking the first cell centre from ~50 m to
  ~6 m raises the ratio 6.6× (0.08 → 0.53); at 50 m even removing drag only reaches 0.28. The
  corner-flow layer simply is not representable by a 50 m first cell.
- **But the *drag parameterisation* is what blocks the connection.** Within the fine-mesh group
  (all dz₁ = 6.2 m — a **clean one-variable comparison**) the ratio is 0.53 with bulk drag, 0.51
  with the height-consistent log-law drag, and **0.88 with drag off** — the only case that meets
  the surface-connection criterion, and it also has the strongest near-surface convergence
  (5.6×10⁻³ vs 3.4×10⁻³). Roughness magnitude is nearly irrelevant (0.004–0.024 all ≈ 0.08 at 50 m),
  and surface heat/moisture fluxes change nothing.

**Interpretation (honest).** In nature surface drag *drives* the corner-flow inflow that concentrates
angular momentum (Rotunno–Klemp). In this model the drag is applied as an implicit damping of the
lowest cell's **full** horizontal wind — including its **tangential** component — so it removes
near-surface angular momentum without generating the compensating **radial** inflow. It is therefore
a net *sink* of near-surface rotation rather than the source of the corner flow. Making C_d
height-consistent (log-law) does not fix this, because the defect is the *form* of the closure, not
its coefficient.

### The fix: a surface-layer STRESS-DIVERGENCE closure — surface connection with drag ON

The defect above is the closure's *form*: the bulk sink rate is ``C_d|V|/dz₁``, which **diverges as
the mesh is refined** (measured: the lowest cell retains 0.994 at dz₁ = 19.6 m but only 0.936 at
dz₁ = 1.8 m — 10× more momentum stripped, mostly *tangential*). The physical form spreads the stress
over a **physical** depth: ``tau(z) = tau_s(1 - z/h)`` ⇒ ``du/dt = -tau_s/h``, uniform through the
layer and **mesh-independent** (measured: 0.9984 at both resolutions).
Implemented as `surface_drag.apply_surface_stress_divergence`
(`SurfaceDragConfig.stress_divergence`, `surface_layer_depth_m`; opt-in, default byte-identical).

At the refined mesh (all dz₁ = 6.2 m — a clean one-variable comparison of the **closure**):

| closure | C_d applied | V_sfc | **sfc/aloft** | connected |
|---|---|---|---|---|
| bulk (lowest-cell damping) | 0.0120 | 1.25 | 0.53 | ✗ |
| log-law C_d, bulk form | 0.0094 | 1.28 | 0.51 | ✗ |
| **stress-divergence** | 0.0120 | 3.02 | 0.69 | ✗ |
| **stress-divergence + log-law** | 0.0094 | **3.07** | **0.82** | **✓** |
| *(drag off — reference)* | 0 | 3.25 | 0.88 | ✓ |

**Fixing the form — not removing the drag — restores the surface connection.** V_sfc rises
1.25 → 3.07 (2.4×, essentially the drag-free reference) with the drag still physically active, and
the near-surface convergence rises 3.4×10⁻³ → 5.3×10⁻³. The remaining ingredient beyond the form is
the height-consistent C_d: `stress-divergence + log-law` (0.82) clears the criterion while
`stress-divergence` with the too-large fixed C_d (0.69) does not.

**Honest scope.** This closes the *structural* gap — the vortex is no longer elevated, and the
surface-connection criterion is met with realistic surface friction. It is **not** a tornado: these
are short runs on a modest horizontal mesh and V_sfc ≈ 3 m s⁻¹, far from the observed 26. The
`SURFACE_CONNECTED_TORNADO_LIKE_VORTEX` tier additionally requires the intensity, pressure-deficit
and persistence criteria in `classification.py`. What is now established is that the model *can*
carry rotation to the lowest resolved level with physical drag — the mechanism that was blocked.

### Attempt K — everything combined: the profile INVERTS (surface-intensified)

`scratchpad/tornado_intensity_gpu.py`. All the arc's results in one run: the **freely-evolving**
supercell (correct geometry, I) + **stress-divergence + log-law** drag (the closure fix) + a nest
that resolves the **surface layer vertically** (dz₁ = 5.1 m, vs 40 m on the parent) + the **moving
nest** (30 re-centrings, so the vortex cannot drift out). The vertical CFL forbids a fine
near-surface parent for a 45-min maturation (dz₁ ≈ 4 m ⇒ dt ≈ 0.05 s), so the surface layer is
refined only inside the nest — verified to work across the interface.

| z | 5.1 m | 52.8 m | 203 m | 522 m | 968 m | 1462 m |
|---|---|---|---|---|---|---|
| **V_rot** | **7.76** | 7.54 | 6.94 | 6.40 | 5.56 | 6.71 |
| \|ζ\| (10⁻²) | 3.60 | 3.66 | 3.89 | 4.28 | 4.34 | 4.95 |

**The profile has inverted.** V_rot is now **maximum at the lowest resolved level and decreases
upward** — against Attempt J's elevated profile (2.7 at 40 m rising to 7.7 at 1.5 km). Surface V_rot
is **2.9× Attempt J's** and the largest near-surface value in the study; the criterion reports
surface-connected. The structural chain — free evolution → aligned streamwise geometry → physical
surface closure → resolved corner-flow layer → rotation carried to the ground — is now complete.

**Still not a tornado, and two honest caveats:**
- **Intensity**: V_sfc ≈ 7.8 m s⁻¹ is ~30 % of the observed 26; the classifier still returns
  `LOW_LEVEL_MESOCYCLONE` (the intensity/pressure-deficit/persistence criteria are unmet), on a
  67 m horizontal mesh over a 200 s window.
- **A diagnostic bug was caught and fixed here.** The run first reported V_rot = 37.9 at the surface
  (ratio 1.00) — an artifact: `surface_connection_report` used a `nx//6` interior margin that reaches
  into the nest's boundary-relaxation zone. Independent checks contradicted it (|ζ| barely changed
  with height, Δp = 0, `vortex_report` gave V_θ = 12). The margin is now `border_frac = 0.2` and the
  honest value is **7.76**. *(Δp was unusable on nests at the time — see the next section, which
  fixes it.)*

*(Caveat: comparisons **across** the coarse/fine groups are indicative — `z_stretch` redistributes
all levels, so V_aloft differs (2.4–3.7 vs ~10). Comparisons **within** the fine group are clean.)*

---

### The pressure defect — a diagnostic that silently read zero on every nest

Noted in Attempt H and carried as a caveat through K: on a nest, `state.p` had a horizontal spread of
**exactly 0**, so every pressure-based number — the vortex pressure deficit above all — read zero.
The cause was not physics. The **low-memory anelastic projection computed the pressure potential φ
and discarded it**, returning only the velocity correction and the divergence residual. Any grid
above `_LOWMEM_N = 64 000` cells routes through that solver, which is *every fine nest in the entire
study*. `state.p` simply kept whatever stale value it had.

What this cost, stated precisely (an earlier draft of this section overstated it, and the code
disagrees). The tornado-like test is a *disjunction*:

```python
tlv = (v_theta_max >= 15.0) or (circulation >= 4.0e4 and pressure_deficit <= -200.0)
```

So the tier was **not** unreachable — the `V_θ` branch never depended on pressure, and that is the
branch Attempts A–K failed on its own merits (V_θ ≈ 7.8 against a 15 m s⁻¹ threshold). What was dead
is the **second** branch: a vortex that was broad and deep — large circulation, a real pressure
deficit — but whose peak tangential wind fell below 15 m s⁻¹ could not be recognised on a nest. That
is a plausible signature for an under-resolved tornado, which is exactly the regime these runs are
in, so the loss was not academic; but it does *not* overturn any A–K verdict.

The fix (`5c3c022`) returns φ on request and stores `p = φ/dt`; the factor follows from the two
paths' conventions (direct corrects `u -= (dt/ρ₀)∇p`, low-memory `u -= ∇φ/ρ₀`). Verified against
the direct solver: pressure and velocity agree to **2.5 × 10⁻⁷** relative, and the low-memory solve
is in fact the *more* exact of the two (divergence residual 10⁻¹⁷ against the direct solver's CG
tolerance 8 × 10⁻¹⁰). Two related traps were closed at the same time, both instances of the same
mistake as the `V_rot 37.9` artifact — **referencing a nest diagnostic against its own boundary**:
`vortex_report` still used the default `1/6` search margin (now settable), and `pressure_deficit`
measured "ambient" at the domain *edge*, which on a nest is the relaxation zone (an interior
reference ring is now available via `ambient_frac`).

Method note worth keeping: building the direct-vs-low-memory comparison surfaced a trap that looks
exactly like a solver bug and is not. The low-memory path **forces the boundary faces
BC-consistent before solving** (wrap faces synced when periodic, wall faces zeroed when not); the
direct path does not. Handed the same inconsistent random field, the two solvers therefore project
*different inputs* and disagree at O(1). The comparison is only meaningful once the input is made
BC-consistent up front.

---

### Attempt L — cascade to 22 m: two findings, and a confounded experiment (stopped early)

`scratchpad/tornado_intensity_L_gpu.py`. A third cascade level (600 → 200 → 67 → **22 m**), a longer
window with persistence tracking through the new `sample(sims, t)` hook, and — for the first time —
a measurable pressure deficit on a nest. Every level was diagnosed at the *same instants*, intended
as a controlled resolution experiment. **It is not one, and the run was stopped at t = 40 s of its
180 s window once that was clear.** Two solid findings came out of it anyway.

**Finding 1 — the low-level circulation sits 3.84 km from the mid-level mesocyclone.** Attempts J
and K centred their nests on the mid-level meso (|ζ| peak at z ≈ 500 m). At z ≈ 100 m the peak is
almost 4 km away. A 2.0 km fine nest centred on the meso therefore sat over *quiet air*: |ζ| = 0.002
against 0.238 on its own 67 m parent. Centring every level on the **low-level** peak instead
immediately gave |ζ| = 0.196 s⁻¹ (observed: 0.205) with a 25 m core. The same reasoning applies to
the nest tracker, whose column-max |ζ| follows the mid-level meso — hence the new
`follow_z_lo`/`follow_z_hi`. *Any* small fine nest in this problem must be placed and steered on the
low levels, or it resolves the wrong feature at great expense.

**Finding 2 — the nest pressure deficit is not trustworthy, and now says so.** A nest's projection
must also absorb the imbalance of its *imposed* lateral inflow, and that artifact scales like 1/dt,
so it worsens exactly as the mesh refines. Measured: Δp = −2860 Pa where the cyclostrophic scale
−ρv_θ² for the same vortex is ≈ −120 Pa — a factor of 24. Every Δp is now reported alongside that
scale rather than at face value. (Related: a *constant* pressure field now reports NaN instead of a
deficit of 0.0, which had been passing through as though it were a measurement.)

**Why the experiment is confounded.** A nest must sit inside its parent, so each finer level was
also a *smaller box*: 13.2 → 4.4 → 2.0 km. `dx` and domain size vary together, and the two metrics
disagree in a way that proves the point — peak |ζ| *rises* with refinement while mean V_rot *falls*:

| dx | domain | mean V_rot | peak V_rot | peak \|ζ\| | final ratio |
|---|---|---|---|---|---|
| 200 m | 13.2 km | 10.93 | 11.66 | 0.075 | 0.83 |
| 67 m | 4.4 km | 8.54 | 10.93 | 0.127 | 1.00 |
| 22 m | 2.0 km | 4.62 | 7.75 | **0.196** | 0.59 |

V_rot is the peak deviation within a fixed radius, so it depends on how much of the circulation fits
inside the box; |ζ| is local and does not. With a 0.2 border margin the 22 m level had only ~1.2 km
of trusted interior, so its vortex was boundary-controlled — its surface/aloft ratio decayed 1.00 →
0.59 (elevated again) while the 67 m level held 1.00. **No resolution conclusion can be drawn from
this table**, which is why the matched-domain test below exists.

### The matched-domain resolution test — `dx` alone, at a fixed 4.0 km domain — VERDICT: refinement does NOT intensify

`scratchpad/tornado_matched_domain_gpu.py`. Two cascades from the same cached parent at the same
instant, centred identically on the low-level circulation, scored with the same fixed 400 m radius
and 0.2 interior margin, whose **finest domain is 4.0 km in both**: `coarse` = 600/200/**67 m**
(60 cells), `fine` = 600/200/67/**22 m** (180 cells). Only `dx` differs; both use the moving nest
(`follow_interval=8`, low-level window 0–1500 m), stress-divergence + log-law drag, `z1 = 5.1 m`.
4.0 km is the smallest box that keeps both branches ≥ 60 cells across. Coarse completed in 15 min,
fine in 174 min; 45 matched samples each. *A first attempt at 2.0 km NaN-ed the coarse branch on
its first step (30 cells wide = no interior), which is what forced the retarget.*

**The answer to the narrow question is NO.** At a fixed physical domain, refining 67 → 22 m does
not intensify the surface vortex — it is equal-to-weaker:

| | peak V_sfc | mean V_sfc | peak \|ζ\| | mean \|ζ\| | mean ratio | peak v_θ |
|---|---|---|---|---|---|---|
| coarse 67 m | **12.28** | **6.35** | 0.134 | **0.075** | 0.78 | **8.29** |
| fine 22 m | 9.05 | 5.65 | **0.158** | 0.043 | **0.89** | 3.47 |

The coarse branch's intensification event (t ≈ 29–33 s → V_sfc 12.28, holding ~9 to t = 60) is the
whole result; the fine branch *passed through that window* showing only a bump (5.6 at t = 31.6)
and peaked later, at t = 42.2 (9.05, ζ 0.158, core 125 m — a genuine interior tight-core event),
then decayed monotonically to 6.65. The field snapshot at t = 60 confirms it: interior |ζ| 0.140
(coarse) vs 0.065 (fine), max wind in a fixed 400 m disc 10.6 vs 9.4 m s⁻¹.

So the tentative reading that the coarse *domain* alone had made Attempt L's 22 m level weak is
incomplete — at the same 4.0 km box the 22 m level is still not stronger than 67 m. **Horizontal
mesh resolution joins the measured-and-eliminated list** (mesh D, environment E, updraft F,
cold-pool source H): ~67 m is already sufficient to carry the rotation the storm supplies, and the
remaining limiter is the storm-scale angular-momentum supply feeding the corner flow — the one
lever not yet measured. Marginal upsides of 22 m worth noting: better connection (ratio 0.89 vs
0.78) and one brief finer-core event (peak ζ 0.158, the largest of either run); if a tornado-grade
core (Δv ≫ 50 m s⁻¹) ever forms, 22 m-class mesh will be needed to *resolve* it — it just does not
*produce* one here.

**Three new measurement traps found while closing this out** (all reproducible from
`outputs/matched_domain/`):

1. **The relaxation zone manufactures tornado-grade vorticity.** With the border *unmasked*, both
   final fields peak at the box corners/edges: |ζ| = 0.176 with 16 m s⁻¹ winds (67 m), **0.271 with
   23 m s⁻¹ (22 m)** — larger on the finer mesh, and 2–4× the interior maxima (0.140 / 0.065).
   Boundary-blending shear, not the vortex. `border_frac = 0.2` interior-only scoring is not
   cosmetic; an unmasked report here would have announced a violent tornado sitting in a nest corner.
2. **`dP` is NaN on every nest-move step** (period ~10.5 s = `follow_interval` in both runs): the
   first diagnostic instant after a regrid has no valid projection pressure. Skip move-step samples
   rather than averaging around them.
3. **The two reports can disagree violently on the same state** — at the fine peak (t = 42.2):
   V_sfc 9.05, ζ 0.158, core 125 m, but `vortex_report` v_θ = 0.94 m s⁻¹. The center-finders hop
   between neighbouring interior maxima (the series alternates two regimes, ζ ≈ 0.01 / r = 1.00 vs
   ζ ≈ 0.045 / r ≈ 0.85, at the re-centring period). Cross-check any single-instant claim against
   the saved fields (`fields_{coarse,fine}.npz`) before believing it.

#### RETRACTION — the "refinement does not intensify" verdict above is NOT SAFE

Following up trap 1 turned it from a *scoring* caveat into a *simulation* defect. Locating every
coherent vortex in the two final fields (a |ζ| maximum that also sits on a real `p_dyn`
depression, rather than any |ζ| maximum) gives:

| | 67 m | 22 m |
|---|---|---|
| coherent vortices found | **4** (at 333, 667, 1000, 1133 m from the edge; `dp_local` −99, −39, −73, −21 Pa) | **1** (at **111 m**, `dp_local` −30 Pa) |
| top-5 \|ζ\| maxima that are boundary junk | 1 of 5 | **4 of 5** (`dp_local` **+193**, +26, −10, −7 Pa) |
| mean \|ζ\|, 0–4 cells from edge | 0.0083 | **0.0343** |
| mean \|ζ\|, 4–12 cells | 0.0084 | 0.0218 |
| mean \|ζ\|, 12 cells–centre | 0.0095 | **0.0069** |

The 67 m nest is healthy: |ζ| is flat with distance from the boundary, and its edge feature is a
real anticyclone (`p_dyn` −118 Pa), the couplet's other member. The 22 m nest is not. It carries
**5× more vorticity at its boundary than in its interior**, its global maximum sits on a `p_dyn`
*maximum* (+368 Pa — not a vortex at all), and its **one** coherent vortex ended up 111 m from the
wall, i.e. *inside its own 89 m sponge*, where it was being nudged toward the parent state. It was
simultaneously damped by the relaxation and excluded from the interior scoring window.

So the fine branch's 9.05 is not "22 m resolves a weaker vortex"; it is "22 m measured the flank of
a vortex it had pushed into its own wall." **Horizontal mesh resolution is therefore NOT eliminated
— that row is withdrawn from the measured-and-eliminated list pending the re-test below.** What
does survive the pair unaffected: circulation agrees to 3 % (7.93e3 vs 8.16e3), and Attempt K's
flat-to-surface connected profile reproduces independently at both meshes.

Two root causes, both now fixed in `nesting.py` (opt-in, defaults byte-identical, 7 tests in
`tests/test_nest_sponge_and_follow.py`):

1. **`relax_width` is a cell COUNT, so the sponge shrinks physically with every refinement level**
   — 4 cells = 267 m at 67 m but only 89 m at 22 m. Damping crossing the band goes like
   `relax_rate * width / U`, so the fine nest damped incoming error ~3× less while having more of
   it to absorb. Fix: `NestSpec.relax_width_m` pins the band to a physical width
   (`effective_relax_width()`); `relax_width_m=None` reproduces the old behaviour exactly.
2. **The moving-nest tracker was steered by the artifact.** `tag_cells` thresholds at
   `frac * domain max`, so when the tracked grid is itself a nest, its sponge vorticity both
   inflates the normaliser and drags the tagged cluster into the wall — a feedback that pushes the
   real vortex toward the opposite boundary. Fix: `tag_cells(border=)` / `follow_spec(border=)`,
   plumbed as `run_multilevel_nest(follow_border=)` with `"auto"` = exclude a parent-nest's own
   band. This also better explains the fine run's two alternating regimes than the
   "center-finder hopping" reading in trap 3 above.

Re-test: `scratchpad/tornado_matched_domain_v2_gpu.py` (`sponge=267 m`, `follow_border=auto`,
`follow_filter=0.7`). Note the sponge fix is a **no-op at 67 m** (267 m *is* 4 cells there), so the
coarse branch re-run isolates the tracker fix alone — and it does change that branch, so the
previously reported coarse numbers (peak V_sfc 12.28) are themselves provisional.

#### The re-test came back: the fixes were NOT sufficient

Both v2 branches completed (67 m in 33 min, 22 m in 317 min). The two fixes are real and the
67 m branch is now demonstrably healthier, but **they did not make the 22 m branch measurable,
and they cost the headline number.**

**v2 at 67 m (only the tracker fix bites here -- 267 m *is* 4 cells at this mesh):**

| | v1 | v2 |
|---|---|---|
| peak V_sfc | **12.28** | **8.47** |
| peak \|zeta\| | 0.134 | 0.067 |
| peak v_theta | 8.29 | 5.21 |
| final ratio / connected | 0.94 / True | **1.00** / True |
| final v_theta | 0.18 | **2.25** |
| mean \|zeta\| edge / 4-12 / interior | 0.0083 / 0.0084 / 0.0095 | **0.0051 / 0.0060 / 0.0066** |
| edge/interior | 0.87 | **0.77** |

The field is healthier by every independent structural measure -- edge vorticity now *below* the
interior and rising monotonically inward, final v_theta 2.25 against v1's near-zero 0.18 -- while
peak V_sfc falls 31 % and peak |zeta| halves. **So 12.28 m/s is retracted**: it came from a run
whose tracker was steered by sponge vorticity, and the t = 29-33 s intensification event that
produced it does not survive tracking the true low-level vortex. The defensible best surface
rotation is now **8.47 m/s (33 % of the observed 26)**, and Attempt K's surface-intensified
profile becomes the more robust of the two headline claims.

**v2 at 22 m -- still not a valid measurement:**

| | v1 fine | v2 fine |
|---|---|---|
| mean \|zeta\| edge 0-4 | 0.0343 | **0.0148** |
| mean \|zeta\| 4-12 | 0.0218 | 0.0161 |
| mean \|zeta\| 12-centre | 0.0069 | **0.0026** |
| **edge/interior** | 4.98 | **5.61 (WORSE)** |
| strongest coherent vortex | 111 m from edge | **178 m from edge** |

Absolute edge vorticity did drop by more than half (the wider sponge damps better, as designed),
but the **interior emptied faster**, so contamination got worse in ratio; and the one coherent
vortex moved from 111 m to 178 m from the wall -- i.e. *deeper inside the now-267 m sponge*.
Widening the sponge to the physically correct width buried the vortex in it. The 67 m branch by
the same measure stays healthy (edge/interior 0.77, vortex 933 m in).

**Conclusion: the dominant defect is neither the sponge width nor the tracker.** The 3-level fine
cascade cannot keep the low-level vortex in trusted interior, and no fix so far addresses that.
The resolution question stays OPEN; mesh resolution stays OFF the eliminated list.

**A diagnostic caveat on the section above.** The "coherent vortex" test used here keys on a local
`p_dyn` depression. On the v2 *coarse* field those local values run -1418, -3250, -5161 and
**+7412 Pa** -- the known nest-projection artifact (`p_dyn ~ 1/dt`), not physics. A -20 Pa
threshold is meaningless against that noise, so the "6 coherent vortices" figure quoted for v2
coarse is NOT trustworthy. The v1 values were in a sane range (-99, -39, -73, -21 Pa), which is
why the test looked reliable when it was built. Trust the edge/interior |zeta| ratios (which use
no pressure at all) over the vortex counts.

**Next, before another full window** (`scratchpad/nest_tracking_diagnosis_gpu.py`): a SHORT (12 s)
instrumented run that reports, for *every* cascade level each sample, where the low-level |zeta|
peak sits inside that level's own box, whether it has entered that level's sponge, and the level's
edge/interior ratio. The first level whose peak crosses into its sponge is the one to fix. A full
fine window costs 317 min; this costs ~1 h and says what to change.

---

## 6. Reproduce

```
# observed target (WSL2: pip install --user arm_pyart nexradaws ; numpy<2.3)
python3 deploy/wsl2_nexrad_moore.py        # read KTLX Level II, reflectivity + couplet
python3 deploy/wsl2_extract_velocity.py    # extract the observed V_rot / Δv

# simulated attempts (this machine, GPU)
python scratchpad/moore_real_funnel_gpu.py # Attempt A: 3-level AMR to 46 m
python scratchpad/moore_sustained_gpu.py   # Attempt B: sustained maturation + follow nest
python scratchpad/moore_forced_gpu.py      # Attempt C: + sustained ascent forcing (§4)
python scratchpad/moore_cascade_gpu.py     # Attempt D: storm-relative deep cascade to 28 m
python scratchpad/moore_koun_cascade_gpu.py    # Attempt E: real KOUN sounding (SRH 254)
python scratchpad/moore_fineparent_gpu.py      # Attempt F: resolve the updraft (fine 250 m parent)
python scratchpad/moore_twoway_ab.py           # Attempt G: two-way coupling A/B (the lever)
python scratchpad/moore_twoway_deep_gpu.py     # Attempt H: two-way deep cascade + vorticity budget
python scratchpad/supercell_alignment_evolution.py  # Attempt I: freely-evolving supercell, alignment(t)
python scratchpad/tornado_occlusion_gpu.py     # Attempt J: resolve the vortex to 22 m (elevated)
python scratchpad/tornado_intensity_gpu.py     # Attempt K: all combined -> surface-INTENSIFIED
DEV=gpu WIN=240 python scratchpad/tornado_intensity_L_gpu.py   # Attempt L: 3-level to 22 m,
                                               #   persistence tracking + real Delta p on the nest

# feature tests
python -m pytest tests/test_forcing.py tests/test_vorticity_budget.py tests/test_vortex_diagnostics.py -q
python -m pytest tests/test_lowmem_pressure.py -q     # the nest pressure fix (p = phi/dt)
```
