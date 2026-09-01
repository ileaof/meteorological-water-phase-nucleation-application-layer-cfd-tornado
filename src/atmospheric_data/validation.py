"""Validation metrics against radar / surface obs (ROADMAP §3a).

Standard scores (task list): RMSE, MAE, bias, correlation, Fractions Skill Score (FSS),
Critical Success Index (CSI), and spatial/temporal displacement of a feature (e.g. the
mesocyclone).  All operate on matched arrays (obs vs simulated) with NaN-aware masking.
"""
from __future__ import annotations

import numpy as np


def _pair(obs, sim):
    o = np.asarray(obs, float).ravel(); s = np.asarray(sim, float).ravel()
    m = np.isfinite(o) & np.isfinite(s)
    return o[m], s[m]


def rmse(obs, sim):
    o, s = _pair(obs, sim); return float(np.sqrt(np.mean((s - o) ** 2))) if o.size else np.nan


def mae(obs, sim):
    o, s = _pair(obs, sim); return float(np.mean(np.abs(s - o))) if o.size else np.nan


def bias(obs, sim):
    o, s = _pair(obs, sim); return float(np.mean(s - o)) if o.size else np.nan


def correlation(obs, sim):
    o, s = _pair(obs, sim)
    if o.size < 2 or o.std() == 0 or s.std() == 0:
        return np.nan
    return float(np.corrcoef(o, s)[0, 1])


def critical_success_index(obs, sim, threshold):
    """CSI = hits / (hits + misses + false alarms) for exceedance of ``threshold``."""
    o, s = _pair(obs, sim)
    O = o >= threshold; S = s >= threshold
    hits = np.sum(O & S); miss = np.sum(O & ~S); fa = np.sum(~O & S)
    den = hits + miss + fa
    return float(hits / den) if den else np.nan


def fractions_skill_score(obs, sim, threshold, window=3):
    """FSS for a 2-D field at exceedance ``threshold``, neighbourhood ``window`` (odd) cells."""
    from scipy.ndimage import uniform_filter
    o = (np.asarray(obs, float) >= threshold).astype(float)
    s = (np.asarray(sim, float) >= threshold).astype(float)
    of = uniform_filter(o, window, mode="constant"); sf = uniform_filter(s, window, mode="constant")
    mse = np.mean((of - sf) ** 2); ref = np.mean(of ** 2) + np.mean(sf ** 2)
    return float(1.0 - mse / ref) if ref > 0 else np.nan


def feature_centroid(field, threshold):
    """(i,j) intensity-weighted centroid of the region above ``threshold`` (or None)."""
    f = np.asarray(field, float); m = f >= threshold
    if not m.any():
        return None
    idx = np.argwhere(m); wts = f[m]
    return tuple((idx * wts[:, None]).sum(0) / wts.sum())


def displacement_error(obs, sim, threshold, spacing=1.0):
    """Distance between the obs and sim feature centroids (× ``spacing`` for physical units)."""
    co = feature_centroid(obs, threshold); cs = feature_centroid(sim, threshold)
    if co is None or cs is None:
        return np.nan
    return float(np.hypot(co[0] - cs[0], co[1] - cs[1]) * spacing)


def radar_metrics(radar_obs, vr_sim, refl_sim=None, refl_thresh=20.0, vr_couplet_thresh=15.0):
    """Bundle the standard scores comparing observed radar to the synthetic-radial CFD."""
    out = {"radial_velocity": {"rmse": rmse(radar_obs["radial_velocity"], vr_sim),
                               "mae": mae(radar_obs["radial_velocity"], vr_sim),
                               "bias": bias(radar_obs["radial_velocity"], vr_sim),
                               "correlation": correlation(radar_obs["radial_velocity"], vr_sim)}}
    if refl_sim is not None and radar_obs.get("reflectivity") is not None:
        out["reflectivity"] = {"rmse": rmse(radar_obs["reflectivity"], refl_sim),
                               "csi": critical_success_index(radar_obs["reflectivity"], refl_sim, refl_thresh),
                               "fss": fractions_skill_score(radar_obs["reflectivity"], refl_sim, refl_thresh)}
    # mesocyclone displacement from the |Vr| couplet
    out["mesocyclone_displacement_gates"] = displacement_error(
        np.abs(radar_obs["radial_velocity"]), np.abs(vr_sim), vr_couplet_thresh)
    return out
