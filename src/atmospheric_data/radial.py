"""Radar radial-velocity observation operator (ROADMAP §3a, source C validation).

A Doppler radar measures ``V_r = V . r_hat`` -- the projection of the wind onto the beam, NOT
the 3-D vector.  To compare the CFD honestly we do the SAME projection: interpolate the model
``(u,v,w)`` to the radar gate positions and dot with the gate's radial unit vector from the
radar.  (Going the other way -- recovering 3-D wind from one radar -- is under-determined; that
is why we compare in radial space, per the task's explicit caution.)
"""
from __future__ import annotations

import numpy as np


def radial_unit_vectors(gate_x, gate_y, gate_z, radar_xyz):
    """Unit vectors from the radar to each gate (model projected frame, metres)."""
    rx = np.asarray(gate_x, float) - radar_xyz[0]
    ry = np.asarray(gate_y, float) - radar_xyz[1]
    rz = np.asarray(gate_z, float) - radar_xyz[2]
    r = np.sqrt(rx * rx + ry * ry + rz * rz) + 1e-9
    return rx / r, ry / r, rz / r


def project_to_radial(u, v, w, gate_x, gate_y, gate_z, radar_xyz):
    """``V_r = u r_x + v r_y + w r_z`` at each gate (all args already at gate positions)."""
    rx, ry, rz = radial_unit_vectors(gate_x, gate_y, gate_z, radar_xyz)
    return (np.asarray(u, float) * rx + np.asarray(v, float) * ry + np.asarray(w, float) * rz)


def cfd_radial_velocity(fields, x_model, y_model, z_model, radar, it=0):
    """Synthetic radial velocity from the model fields at the radar's gate positions.

    ``fields`` holds ``u,v,w`` as ``(time,nz,ny,nx)`` on model axes ``(x_model,y_model,
    z_model)`` in the projected frame centred at the domain centre; ``radar`` is a volume dict
    (``x_m,y_m,alt_m`` gate positions relative to the radar, plus radar position).  Returns
    ``V_r`` on the radar gate grid, matching ``radar['radial_velocity']`` in shape."""
    from scipy.interpolate import RegularGridInterpolator
    gx = np.asarray(radar["x_m"], float) + _radar_xy(radar)[0]
    gy = np.asarray(radar["y_m"], float) + _radar_xy(radar)[1]
    gz = np.asarray(radar["alt_m"], float)
    pts = np.stack([np.clip(gz, z_model.min(), z_model.max()).ravel(),
                    np.clip(gy, y_model.min(), y_model.max()).ravel(),
                    np.clip(gx, x_model.min(), x_model.max()).ravel()], -1)
    interp = lambda F: RegularGridInterpolator((z_model, y_model, x_model), F[it],
                                               bounds_error=False, fill_value=None)(pts).reshape(gx.shape)
    u = interp(fields["u"]); v = interp(fields["v"]); w = interp(fields["w"])
    radar_xyz = (_radar_xy(radar)[0], _radar_xy(radar)[1], radar.get("radar_alt_m", 0.0))
    return project_to_radial(u, v, w, gx, gy, gz, radar_xyz)


def _radar_xy(radar):
    """Radar position in the model (domain-centred projected) frame.  The synthetic radar sits
    at the domain centre (0,0); a real radar's offset is its projected (x,y)."""
    return radar.get("radar_x_m", 0.0), radar.get("radar_y_m", 0.0)
