# NEXRAD radial‑velocity validation (`atmospheric_data.radial` + `.validation`)

A WSR‑88D Doppler radar measures **radial velocity** `V_r = V · r̂` — the projection of the
wind onto the beam — and reflectivity, not the full 3‑D wind vector. Recovering 3‑D wind from a
**single** radar is under‑determined, so we validate the model **in radial space** by applying
the same projection to the CFD.

## The observation operator

For each radar gate at position `g` (model projected frame) with radial unit vector
`r̂ = (g − radar)/|g − radar|`, interpolate the model `(u,v,w)` to `g` and compute

```
V_r_sim = u r̂_x + v r̂_y + w r̂_z          (atmospheric_data.radial.cfd_radial_velocity)
```

This matches the observed `V_r` grid gate‑for‑gate. (Do **not** treat `V_r` as the 3‑D vector.)

## Reading Level II

`sources.nexrad.load` reads a Level II volume with **Py‑ART** (`pip install arm_pyart`) into a
volume dict (`azimuth, range, elevation, gate lat/lon/alt & x/y, reflectivity, radial_velocity,
spectrum_width, ZDR, ρhv, φdp`). Download is from AWS `noaa-nexrad-level2` (no credentials). If
Py‑ART is absent, `compare-radar` uses a labelled **synthetic** radar so the operator and
metrics still run.

## Metrics (`atmospheric_data.validation`)

`radar_metrics(obs, vr_sim, refl_sim)` returns:

* radial velocity: **RMSE, MAE, bias, correlation**;
* reflectivity: **CSI**, **FSS** (neighbourhood);
* **mesocyclone displacement** from the `|V_r|` couplet centroids (spatial), and the analysis
  time offset (temporal) when comparing across times.

## Reflectivity forward operator (diagnostic)

The idealised microphysics does not emit a calibrated reflectivity, so `atmospheric_data.
reflectivity` computes an **approximate Rayleigh (10 cm)** `Z_e` from the hydrometeor mixing
ratios (`qr,qs,qg`) — `Z_e = a_x (rho q_x)^{1.75}` summed over species (Marshall–Palmer, fixed
intercepts, ice dielectric for snow/graupel), giving ~41 dBZ for 1 g m⁻³ of rain. `compare-radar`
computes it on the model grid, interpolates to the radar gates, and scores CSI/FSS vs the
observed reflectivity. **Assumptions:** dry snow/graupel (no brightband), no attenuation, no
melting layer — indicative, not a validated forward operator.

## Not assimilation

This is an **observation operator + scores**, not data assimilation. Radar DA (3D/4D‑Var,
EnKF) is a separate effort (limitation 7).
