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
| **C** | **sustained forcing** | 250 m | **4.9 m s⁻¹** | 6.6×10⁻³ s⁻¹ | **sustains** (self-holds after forcing off) |
| *observed* | *real* | ~250 m (0.5° beam) | **26 m s⁻¹** | 0.205 s⁻¹ | *tornado* |

*(ζ is Δx-dependent — A's 0.026 at 46 m is not comparable to C's at 250 m; **V_rot**, a velocity,
is the resolution-fair metric. C is the strongest and, uniquely, comes from a storm that **lives
without the crutch**.)*

**Honest read:** V_rot ~4.9 m s⁻¹ is the best low-level rotation we have obtained from real data,
but it is still only **~19 % of the observed 26 m s⁻¹** — *not* a tornado. The lever is validated
(sustained supercell + organising low-level meso), but the remaining gap is now squarely the
things in §5: the 250 m nest still far under-resolves the 0.3 km TVS, the 600 s nest window is
short for full stretching intensification, and ERA5's SRH (152) is well below the real KOUN (247).
The path forward is to **cascade the sustained 250 m nest down to < 50 m at the storm base with a
longer window** — now that we finally have a living parent to nest into.

---

## 5. Honest bottom line and remaining levers

The package reproduces the storm, the rotating updraft, and the mid-level mesocyclone **from real
downloaded data**, and — with the sustained-ascent forcing (§4) — a **self-sustaining supercell
with a monotonically organising low-level mesocyclone**. It does **not** yet reproduce the
tornado-scale low-level vortex: the best low-level rotation so far (Attempt C, V_rot ~4.9 m s⁻¹) is
~19 % of the observed TVS (V_rot 26 m s⁻¹). This is exactly where operational tornadogenesis
science sits — not a bug, the frontier.

Levers, in order of expected impact:

1. **Sustainment** — sustained ascent forcing (§4) so the parent supercell lives long enough for a
   low-level mesocyclone to organise. *(Being tested.)*
2. **Resolution < 30 m at the storm base**, not just aloft, to resolve the ~0.3 km TVS — a finer,
   storm-base-centred nest in the cascade.
3. **Cold-pool microphysics** — the low-level (baroclinic) vorticity source is set by rain
   evaporation; the closure has to land in the narrow not-too-cold/not-too-warm window.
4. **Real low-level SRH** — use the KOUN radiosonde (SRH ~247) rather than ERA5 (~146), which
   smooths the very low-level shear that seeds the rotation.
5. **Storm-following cascade + time-dependent lateral BCs** so the developing low-level meso is
   tracked and fed, not advected out.

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

# feature tests
python -m pytest tests/test_forcing.py -q
```
