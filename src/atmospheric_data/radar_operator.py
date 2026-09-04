"""Radar FORWARD OPERATOR -- turn a model wind field into what a WSR-88D would actually report.

**Why this exists.**  The Moore 2013 validation target, ``V_rot = 26 m/s``, is not a wind speed.
It is half the peak-to-peak Doppler velocity difference across a couplet whose members are
separated by ~253 m, measured by KTLX at a range of ~30 km with a 0.5 deg sweep.  At that range
the 0.925 deg beam is ~480 m across and the sample volume is ~460 m above ground.  A tornado
with a 125 m radius of maximum wind is therefore SMALLER THAN THE BEAM: the radar cannot resolve
it and necessarily under-reports its peak wind.

Comparing a model's grid-point wind against that number is not a validation -- the two are
different quantities.  This module makes the comparison honest in the only direction that is
well posed: push the MODEL through the radar's sampling, and compare observable with observable.

A useful corollary falls out of the same machinery: for an assumed core radius, the operator
gives the UNDER-READING FACTOR, i.e. what true peak wind an observed 26 m/s implies.

**What is modelled**

* 4/3-earth beam height  h(r, theta) (standard atmospheric refraction);
* Gaussian two-way beam weighting across the 3 dB beamwidth in azimuth and elevation;
* range-gate averaging over the pulse length;
* projection onto the beam, ``V_r = V . r_hat`` (radial component only -- a single radar cannot
  see the full vector, which is why retrieval is under-determined and comparison is done here);
* the SAME couplet definition as the observation: ``V_rot = (max V_r - min V_r) / 2`` over a
  search window, with the couplet separation reported alongside.

**What is NOT modelled** (stated so results are not over-read): reflectivity weighting within
the volume (a real Doppler estimate is power-weighted, so debris/rain distribution biases it),
Nyquist folding, ground clutter, beam blockage, and the finite dwell/pulse statistics.
"""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

EARTH_RADIUS_M = 6371000.0
FOUR_THIRDS_EARTH_M = (4.0 / 3.0) * EARTH_RADIUS_M


@dataclass
class RadarSpec:
    """A WSR-88D-like radar.  Defaults are KTLX/88D nominal values.

    ``x_m, y_m`` are the radar position in the MODEL frame (metres, same origin as ``grid``).
    """
    x_m: float = 0.0
    y_m: float = 0.0
    alt_m: float = 0.0
    beamwidth_deg: float = 0.925      # WSR-88D 3 dB beamwidth
    gate_length_m: float = 250.0      # super-resolution range gate
    elevation_deg: float = 0.5        # the sweep the Moore V_rot came from
    name: str = "KTLX-like"

    def beam_diameter_m(self, range_m):
        """3 dB beam diameter at a given slant range (scalar in, scalar out; array in, array
        out -- the sweep needs it per gate)."""
        d = np.asarray(range_m, float) * np.radians(self.beamwidth_deg)
        return float(d) if np.ndim(range_m) == 0 else d


def beam_height_m(range_m, elevation_deg, radar_alt_m: float = 0.0,
                  effective_earth_m: float = FOUR_THIRDS_EARTH_M):
    """Height of the beam centre above ground, 4/3-earth refraction model.

        h = sqrt(r^2 + a^2 + 2 r a sin(elev)) - a + h_radar
    """
    r = np.asarray(range_m, float)
    e = np.radians(float(elevation_deg))
    a = float(effective_earth_m)
    return np.sqrt(r * r + a * a + 2.0 * r * a * np.sin(e)) - a + float(radar_alt_m)


def range_for_beam_height(target_height_m, elevation_deg, radar_alt_m: float = 0.0,
                          effective_earth_m: float = FOUR_THIRDS_EARTH_M):
    """Invert :func:`beam_height_m` -- the slant range at which the beam centre sits at a
    given height.  Used to place a model domain at the observation's sampling height."""
    h = float(target_height_m) - float(radar_alt_m)
    e = np.radians(float(elevation_deg))
    a = float(effective_earth_m)
    # solve r^2 + 2 a sin(e) r - (h^2 + 2 a h) = 0
    disc = (a * np.sin(e)) ** 2 + (h * h + 2.0 * a * h)
    return float(-a * np.sin(e) + np.sqrt(disc))


def _gauss_weights(n: int, beamwidth_deg: float):
    """Sub-beam offsets [deg] and two-way Gaussian weights across the 3 dB beamwidth.

    A Gaussian beam has one-way power ``exp(-4 ln2 (phi/theta_3dB)^2)``; the offsets span
    +-1 beamwidth, which captures the great majority of the weight."""
    if n <= 1:
        return np.array([0.0]), np.array([1.0])
    off = np.linspace(-1.0, 1.0, n) * float(beamwidth_deg)
    w = np.exp(-4.0 * np.log(2.0) * (off / float(beamwidth_deg)) ** 2)
    return off, w / w.sum()


def synthetic_sweep(u, v, w, x_m, y_m, z_m, radar: RadarSpec,
                    azimuths_deg=None, ranges_m=None,
                    n_beam: int = 5, n_gate: int = 3):
    """Sample a model wind field the way the radar would, on a polar (azimuth, range) grid.

    ``u, v, w`` are 3-D arrays on the model axes ``x_m, y_m, z_m`` (1-D, increasing, metres,
    same frame as ``radar.x_m/​y_m``).  Returns a dict with the gate geometry and ``v_r``.

    Each gate is a BEAM-VOLUME AVERAGE: ``n_beam`` x ``n_beam`` Gaussian-weighted sub-rays in
    azimuth and elevation, and ``n_gate`` samples along the pulse.  Setting ``n_beam=1,
    n_gate=1`` recovers naive point sampling, which is what makes the smoothing effect
    measurable rather than assumed.
    """
    from scipy.interpolate import RegularGridInterpolator

    x_m = np.asarray(x_m, float); y_m = np.asarray(y_m, float); z_m = np.asarray(z_m, float)
    interp = {nm: RegularGridInterpolator((x_m, y_m, z_m), np.asarray(F, float),
                                          bounds_error=False, fill_value=None)
              for nm, F in (("u", u), ("v", v), ("w", w))}

    if azimuths_deg is None or ranges_m is None:
        raise ValueError("azimuths_deg and ranges_m are required (use sweep_grid_for_domain)")
    az = np.asarray(azimuths_deg, float)
    rg = np.asarray(ranges_m, float)
    AZ, RG = np.meshgrid(az, rg, indexing="ij")

    d_off, d_wt = _gauss_weights(n_beam, radar.beamwidth_deg)      # azimuth + elevation offsets
    g_off = (np.linspace(-0.5, 0.5, n_gate) * radar.gate_length_m if n_gate > 1
             else np.array([0.0]))
    g_wt = np.full(g_off.shape, 1.0 / g_off.size)

    num = np.zeros(AZ.shape); den = 0.0
    for ia, a_off in enumerate(d_off):
        for ie, e_off in enumerate(d_off):
            wgt_ae = d_wt[ia] * d_wt[ie]
            elev = radar.elevation_deg + e_off
            azi = AZ + a_off
            for ig, r_off in enumerate(g_off):
                rr = RG + r_off
                hh = beam_height_m(rr, elev, radar.alt_m)
                # slant range -> horizontal ground distance (small-angle over these ranges)
                gd = rr * np.cos(np.radians(elev))
                gx = radar.x_m + gd * np.sin(np.radians(azi))     # azimuth measured from +y (N)
                gy = radar.y_m + gd * np.cos(np.radians(azi))
                pts = np.stack([np.clip(gx, x_m.min(), x_m.max()).ravel(),
                                np.clip(gy, y_m.min(), y_m.max()).ravel(),
                                np.clip(hh, z_m.min(), z_m.max()).ravel()], -1)
                uu = interp["u"](pts).reshape(AZ.shape)
                vv = interp["v"](pts).reshape(AZ.shape)
                ww = interp["w"](pts).reshape(AZ.shape)
                # radial unit vector from the radar to this sample
                rx = gx - radar.x_m; ry = gy - radar.y_m; rz = hh - radar.alt_m
                rn = np.sqrt(rx * rx + ry * ry + rz * rz) + 1e-9
                vr = (uu * rx + vv * ry + ww * rz) / rn
                wq = wgt_ae * g_wt[ig]
                num += wq * vr
                den += wq
    vr = num / max(den, 1e-30)

    hh0 = beam_height_m(RG, radar.elevation_deg, radar.alt_m)
    gd0 = RG * np.cos(np.radians(radar.elevation_deg))
    return {"azimuth_deg": AZ, "range_m": RG, "v_r": vr,
            "height_m": hh0,
            "x_m": radar.x_m + gd0 * np.sin(np.radians(AZ)),
            "y_m": radar.y_m + gd0 * np.cos(np.radians(AZ)),
            "beam_diameter_m": radar.beam_diameter_m(RG),
            "radar": {"name": radar.name, "elevation_deg": radar.elevation_deg,
                      "beamwidth_deg": radar.beamwidth_deg,
                      "gate_length_m": radar.gate_length_m,
                      "n_beam": n_beam, "n_gate": n_gate},
            "volume_averaged": bool(n_beam > 1 or n_gate > 1)}


def sweep_grid_for_domain(radar: RadarSpec, x_center_m, y_center_m, half_width_m,
                          az_resolution_deg: float = 0.5):
    """Azimuth/range vectors covering a box of half-width ``half_width_m`` about a model point,
    at the radar's native gate spacing and a given azimuthal sampling."""
    dx = float(x_center_m) - radar.x_m
    dy = float(y_center_m) - radar.y_m
    gd = float(np.hypot(dx, dy))
    az0 = float(np.degrees(np.arctan2(dx, dy)))                 # from +y (north)
    daz = float(np.degrees(np.arctan2(half_width_m, max(gd, 1.0))))
    e = np.radians(radar.elevation_deg)
    r0 = gd / max(np.cos(e), 1e-6)
    az = np.arange(az0 - daz, az0 + daz + 1e-9, az_resolution_deg)
    rg = np.arange(r0 - half_width_m, r0 + half_width_m + 1e-9, radar.gate_length_m)
    return az, rg


def vrot_from_sweep(sweep, max_separation_m: float = 4000.0):
    """``V_rot`` from a synthetic sweep, using the OBSERVATION's definition.

    The reported Moore value is half the peak-to-peak Doppler difference across a couplet
    (inbound -26 / outbound +26 -> ``V_rot = 26``).  We take the maximum and minimum radial
    velocity within the sweep, require them to be within ``max_separation_m`` of each other
    (a couplet, not two unrelated features), and return ``(max - min)/2`` together with the
    separation -- which is directly comparable to the observed ~253 m.
    """
    vr = np.asarray(sweep["v_r"], float)
    if not np.isfinite(vr).any():
        return {"v_rot_m_s": float("nan"), "delta_v_m_s": float("nan"),
                "couplet_separation_m": float("nan"), "valid": False,
                "reason": "no finite radial velocities"}
    i_hi = np.unravel_index(int(np.nanargmax(vr)), vr.shape)
    i_lo = np.unravel_index(int(np.nanargmin(vr)), vr.shape)
    x = np.asarray(sweep["x_m"], float); y = np.asarray(sweep["y_m"], float)
    sep = float(np.hypot(x[i_hi] - x[i_lo], y[i_hi] - y[i_lo]))
    dv = float(vr[i_hi] - vr[i_lo])
    ok = sep <= float(max_separation_m)
    return {"v_rot_m_s": 0.5 * dv if ok else float("nan"),
            "delta_v_m_s": dv if ok else float("nan"),
            "couplet_separation_m": sep,
            "inbound_m_s": float(vr[i_lo]), "outbound_m_s": float(vr[i_hi]),
            "beam_diameter_m": float(np.asarray(sweep["beam_diameter_m"])[i_hi]),
            "sample_height_m": float(np.asarray(sweep["height_m"])[i_hi]),
            "valid": bool(ok),
            "reason": "" if ok else ("extrema %.0f m apart exceed max_separation_m=%.0f -- "
                                     "not a couplet" % (sep, max_separation_m))}


# ------------------------------------------------------------------ analytic validation vortex
def rankine_vortex(x_m, y_m, z_m, center_xy, v_max_m_s: float, core_radius_m: float):
    """A Rankine combined vortex: solid-body inside ``core_radius_m``, 1/r outside.

    Its TRUE V_rot is exactly ``v_max_m_s`` (the peak tangential wind, at r = core radius), and
    its true couplet separation is ``2 * core_radius_m``.  That makes it the reference against
    which the operator's under-reading is measured."""
    X, Y, _ = np.meshgrid(np.asarray(x_m, float), np.asarray(y_m, float),
                          np.asarray(z_m, float), indexing="ij")
    dx = X - center_xy[0]; dy = Y - center_xy[1]
    r = np.sqrt(dx * dx + dy * dy) + 1e-9
    R = float(core_radius_m)
    vth = np.where(r <= R, v_max_m_s * r / R, v_max_m_s * R / r)
    u = -vth * dy / r
    v = vth * dx / r
    return u, v, np.zeros_like(u)


__all__ = ["RadarSpec", "beam_height_m", "range_for_beam_height", "synthetic_sweep",
           "sweep_grid_for_domain", "vrot_from_sweep", "rankine_vortex",
           "EARTH_RADIUS_M", "FOUR_THIRDS_EARTH_M"]
