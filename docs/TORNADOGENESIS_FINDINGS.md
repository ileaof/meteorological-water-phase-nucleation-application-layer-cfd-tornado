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

| Quantity | Observed value |
|---|---|
| Reflectivity core | ~70 dBZ |
| Velocity couplet Δv | **52 m s⁻¹** (inbound −26 / outbound +26) |
| Rotational velocity V_rot | **26 m s⁻¹** |
| Couplet separation | ~253 m (a tight **TVS**) |
| Azimuthal shear ≈ vertical vorticity | **~0.205 s⁻¹** |

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

---

## 5. Honest bottom line and remaining levers

The package reproduces the storm, the rotating updraft, and the mid-level mesocyclone **from real
downloaded data**; with the sustained-ascent forcing (§4) a **self-sustaining supercell**; and with
the storm-relative deep cascade (§4, Attempt D) a **correctly-structured, surface-connected
low-level mesocyclone** — the right feature in the right place. It does **not** yet reproduce the
tornado-scale intensity: the best low-level rotation (Attempt D, V_rot ~6.0 m s⁻¹ at 28 m) is
~23 % of the observed TVS (V_rot 26 m s⁻¹). This is exactly where operational tornadogenesis
science sits — not a bug, the frontier.

**What the resolution study taught us:** refining 9× (250 → 28 m) moved V_rot only +22 %, so the
ceiling is the **source circulation**, not the mesh. The levers, re-ranked by the evidence:

1. **Real low-level SRH** — use the KOUN radiosonde (SRH ~247) rather than ERA5 (~152), which
   smooths the very low-level shear that seeds the rotation. *(Now the top lever.)*
2. **Cold-pool microphysics** — the low-level (baroclinic) vorticity source is set by rain
   evaporation; the closure has to land in the narrow not-too-cold/not-too-warm window. A stronger,
   correctly-placed forward-flank baroclinic zone gives the low-level vortex a stronger parent
   circulation to stretch.
3. **A stronger, longer-lived updraft** — more vertical stretching of the available vorticity
   (longer nest window; a stronger/deeper sustained forcing while the storm organises).
4. **Resolution < 30 m at the storm base** — still needed to *resolve* the ~0.3 km TVS once the
   source circulation is strong enough, but no longer the first lever.

*Achieved so far:* sustainment (§4), the storm-relative deep cascade, and a surface-connected
low-level mesocyclone. *Next:* swap the ERA5 environment for the real KOUN hodograph (lever 1) and
re-run the D cascade.

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

# feature tests
python -m pytest tests/test_forcing.py -q
```
