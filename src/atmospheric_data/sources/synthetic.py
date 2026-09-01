"""Synthetic real-data samples (ROADMAP §3a) -- the backbone for offline tests and demos.

Produces small, physically-plausible stand-ins for the real sources so the WHOLE pipeline
(preprocess -> IC/BC -> QC -> radial validation) runs with **no downloads and no heavy
dependencies** (task: "testes ... não devem depender de grandes downloads").  These are clearly
labelled ``source='synthetic'`` -- never presented as real observations.
"""
from __future__ import annotations

import numpy as np

from .. import thermo
from ..internal import AtmosphericState
from ..project import Projection


def synthetic_atmosphere(cfg, nz=30, ny=24, nx=24, ntime=2):
    """A supercell-like environment on the case domain: veering winds, a moist unstable low
    level, a mid-level updraft anomaly, gentle terrain.  Returns an :class:`AtmosphericState`."""
    dom = cfg.domain
    proj = Projection(dom.center_lat, dom.center_lon, dom.projection)
    half_x = 0.5 * dom.width_km * 1000.0; half_z = dom.height_km * 1000.0
    x = np.linspace(-half_x, half_x, nx); y = np.linspace(-half_x, half_x, ny)
    z = np.linspace(0.0, half_z, nz)
    t0 = np.datetime64("%sT%s" % (cfg.case.date, cfg.case.start_time_utc.zfill(5)))
    times = t0 + (np.arange(ntime) * np.timedelta64(30, "m"))

    # environment profiles (1-D in z), then add a moving mesoscale anomaly in (x,y,t)
    theta0 = 300.0 + 0.0038 * z + 0.012 * np.maximum(z - 12000.0, 0.0)     # troposphere + stable top
    qv0 = 0.014 * np.exp(-z / 3200.0)
    beta = 0.5 * np.pi * np.clip(z / 3000.0, 0, 1)                          # quarter-circle veer
    u0 = 22.0 * np.sin(beta); v0 = -22.0 * (1.0 - np.cos(beta))
    p0, T0, _ = thermo.hydrostatic_base_pressure(z, theta0, qv0, 1.0e5)

    st = AtmosphericState.new(times, z, y, x, projection=dom.projection, source="synthetic")
    X, Y, Z = np.meshgrid(x, y, z, indexing="ij")                          # (nx,ny,nz)
    # assemble full 4-D arrays (time, z, y, x); the storm anomaly drifts east with time
    T4 = np.empty((ntime, nz, ny, nx)); th4 = np.empty_like(T4); qv4 = np.empty_like(T4)
    u4 = np.empty_like(T4); v4 = np.empty_like(T4); w4 = np.empty_like(T4); p4 = np.empty_like(T4)
    for it, tt in enumerate(times):
        cx = -60000.0 + it * 30 * 60 * 12.0
        r2 = ((X - cx) / 20000.0) ** 2 + (Y / 20000.0) ** 2
        bump = np.exp(-r2)
        th = theta0[None, None, :] + 2.5 * bump * np.exp(-Z / 6000.0)
        q = np.clip(qv0[None, None, :] + 0.003 * bump, 0.0, None)
        ww = 12.0 * bump * np.sin(np.pi * np.clip(Z / half_z, 0, 1))
        pp = np.broadcast_to(p0[None, None, :], th.shape)
        TT = thermo.temperature_from_theta(th, pp)
        uu = np.broadcast_to(u0[None, None, :], th.shape) - 6.0 * bump * (Y / 20000.0)
        vv = np.broadcast_to(v0[None, None, :], th.shape) + 6.0 * bump * ((X - cx) / 20000.0)
        tr = lambda a: np.transpose(a, (2, 1, 0))                          # (nx,ny,nz)->(nz,ny,nx)
        T4[it] = tr(TT); th4[it] = tr(th); qv4[it] = tr(q); p4[it] = tr(pp)
        u4[it] = tr(uu); v4[it] = tr(vv); w4[it] = tr(ww)
    for nm, arr in (("T", T4), ("theta", th4), ("qv", qv4), ("p", p4),
                    ("u", u4), ("v", v4), ("w", w4)):
        st.add(nm, arr, source="synthetic", original_name="synthetic_%s" % nm,
               valid_time=str(times[0]))
    terrain = 300.0 + 50.0 * np.exp(-((X[:, :, 0]) ** 2 + (Y[:, :, 0]) ** 2) / (150000.0 ** 2))
    st.add("terrain", np.transpose(terrain, (1, 0)), source="synthetic")   # (y,x), time-invariant
    st.provenance["projection_method"] = proj.method
    return st


def synthetic_pressure_level_state(cfg, nz=11, ny=10, nx=10):
    """A small PRESSURE-LEVEL state (as the HRRR/ERA5 readers produce before height conversion):
    the z coordinate is a pressure proxy and ``p,T,qv`` fields are present -- for testing
    ``sources._common.to_height_levels`` (pressure -> geometric height)."""
    dom = cfg.domain
    plev = np.array([1000, 925, 850, 700, 600, 500, 400, 300, 250, 200, 150.0]) * 100.0  # Pa
    z_proxy = -np.log(plev)                                  # monotone proxy (not height)
    x = np.linspace(-5e4, 5e4, nx); y = np.linspace(-5e4, 5e4, ny)
    t0 = np.datetime64("%sT%s" % (cfg.case.date, cfg.case.start_time_utc.zfill(5)))
    st = AtmosphericState.new(np.atleast_1d(t0), z_proxy, y, x, projection=dom.projection, source="synthetic")
    # a standard-atmosphere-ish column, broadcast over x,y
    theta = 300.0 + 0.03 * (np.arange(nz))                  # increasing potential temperature
    T = thermo.temperature_from_theta(theta, plev)
    qv = 0.012 * (plev / plev[0]) ** 3
    b = lambda a: np.broadcast_to(a[None, :, None, None], (1, nz, ny, nx)).copy()
    st.add("p", b(plev), source="synthetic", units="Pa")
    st.add("T", b(T), source="synthetic")
    st.add("qv", b(qv), source="synthetic")
    st.add("theta", b(theta), source="synthetic")
    return st


def synthetic_sounding():
    """A supercell proximity sounding as radiosonde columns (SI-friendly kwargs)."""
    z = np.array([100, 800, 1500, 3000, 4500, 6000, 8000, 10000, 12000, 14000.0])
    theta = 300.0 + 0.0038 * z
    p, T, _ = thermo.hydrostatic_base_pressure(z, theta, 0.012 * np.exp(-z / 3200.0), 1.0e5)
    beta = 0.5 * np.pi * np.clip(z / 3000.0, 0, 1)
    return {"height_m": z, "pressure_Pa": p, "temperature_K": T,
            "specific_humidity": 0.013 * np.exp(-z / 3200.0),
            "u_ms": 22.0 * np.sin(beta), "v_ms": -22.0 * (1.0 - np.cos(beta))}


def synthetic_radar(cfg, n_az=72, n_gate=60, gate_m=1000.0, elevation_deg=0.5):
    """A tiny NEXRAD-like single-sweep volume about the case domain centre: gate lat/lon/alt,
    azimuth/elevation, reflectivity [dBZ] and radial velocity [m/s] (a rotating couplet)."""
    dom = cfg.domain
    lat0, lon0 = dom.center_lat, dom.center_lon
    az = np.linspace(0, 360, n_az, endpoint=False)
    rng = (np.arange(n_gate) + 1) * gate_m
    AZ, RNG = np.meshgrid(az, rng, indexing="ij")
    el = np.radians(elevation_deg)
    ground = RNG * np.cos(el)
    alt = RNG * np.sin(el) + 380.0
    from ..project import Projection
    proj = Projection(lat0, lon0, dom.projection)
    xg = ground * np.sin(np.radians(AZ)); yg = ground * np.cos(np.radians(AZ))
    lat, lon = proj.to_lonlat(xg, yg)
    # a Rankine-like vortex couplet at ~25 km east: inbound/outbound radial velocity
    xc, yc = 25000.0, 0.0
    dx = xg - xc; dy = yg - yc; r = np.hypot(dx, dy) + 1.0
    vt = 30.0 * np.where(r < 3000.0, r / 3000.0, 3000.0 / r)                # tangential
    ux = -vt * dy / r; uy = vt * dx / r                                    # solid-body-ish rotation
    rhat_x = np.sin(np.radians(AZ)); rhat_y = np.cos(np.radians(AZ))
    vr = ux * rhat_x + uy * rhat_y                                         # radial projection
    refl = np.clip(45.0 * np.exp(-r / 15000.0) + 10.0, 0.0, 70.0)
    return {"azimuth_deg": az, "range_m": rng, "elevation_deg": elevation_deg,
            "lat": lat, "lon": lon, "alt_m": alt, "x_m": xg, "y_m": yg,
            "reflectivity": refl, "radial_velocity": vr,
            "radar_lat": lat0, "radar_lon": lon0, "radar_alt_m": 380.0,
            "station": cfg.data.radar_station, "source": "synthetic"}
