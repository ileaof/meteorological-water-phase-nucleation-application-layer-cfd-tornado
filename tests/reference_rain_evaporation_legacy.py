"""Frozen pre-factor routine; source SHA256: ea5959837de6da9d6136444f5fa1b5a8f22b016fe884864d3e09ff97382fdf7c"""
import math
from precip_microphysics import constants as C, thermo as th, size_distributions as sd
from precip_microphysics.processes import _arr, _cap, _VA, _VB, _ventilated_capacitance, _diffusional_denominator, Transfer

def rain_evaporation(st, cfg, dt):
    """q_r -> q_v in subsaturated air (ventilated, Rutledge & Hobbs 1983)."""
    if not cfg.processes.rain_evaporation:
        return []
    xp = st.xp
    T, P = st.T, st.P
    Sw = th.saturation_ratio_water(st.qv, T, P, xp=xp)
    sub = Sw < 1.0
    if not xp.any(sub & (_arr(st.qr, xp) > C.QSMALL)):
        return []
    lam = sd.lambda_slope(st.qr, st.rho, "rain", xp)
    vent = _ventilated_capacitance(C.N0_r, lam, _VA["rain"], _VB["rain"])
    denom = _diffusional_denominator(T, P, "water", xp)
    # dq_r/dt = (2 pi / rho) (S_w - 1) / (A_K+A_D) * vent   (<0 when subsaturated)
    rate = (2.0 * math.pi / xp.maximum(_arr(st.rho, xp), C.TINY)) * (Sw - 1.0) / denom * vent
    evap = xp.where(xp.isfinite(rate) & (Sw < 1.0), -rate, 0.0)   # positive magnitude
    # do not evaporate past saturation
    to_sat = xp.maximum(th.qsat_water(T, P, xp=xp) - _arr(st.qv, xp), 0.0)
    dq = _cap(xp.minimum(evap * dt, to_sat), st.qr, xp)
    return [Transfer("qr", "qv", dq, "rain_evaporation")] if xp.any(dq > 0) else []
