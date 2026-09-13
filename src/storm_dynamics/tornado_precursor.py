"""Tornado-precursor scoring: the necessary environment, and the geometry that discriminates.

The usual verbal precursor -- surface wind one way, the cloud layer shearing the other,
moisture rising from the ground -- is a statement about HELICITY: a turning hodograph
supplies horizontal vorticity, and low moisture supplies a low cloud base.  Those are the
classical environmental ingredients, and :func:`environment_precursor` measures them.

**This project's own measurements say they are necessary and not sufficient.**  Attempt E
raised storm-relative helicity from 254 to 648 m^2/s^2 and the low-level rotation changed
by about nothing; attempt H then factorised the tilting term and found the horizontal
vorticity abundant (|omega_h| ~ 1.1e-2 at 51 m) and the updraft gradient present, while
the ALIGNMENT cos(theta) sat at +0.04..+0.09 and reached -0.146 at the surface -- tilting
into ANTIcyclonic vorticity.  In attempt I, when the storm was left to occlude freely, the
alignment rose about tenfold and the streamwise fraction went 0.40 -> 0.64, and only then
did |zeta| grow.

So the two halves of this module are deliberately not merged into one number:

* :func:`environment_precursor` -- a 1-D sounding.  Supply side only.  It can say "this
  environment cannot produce a tornado"; it cannot say that it will.
* :func:`state_precursor` -- a 3-D model or analysis state.  It adds the part the sounding
  physically cannot contain, because it needs an updraft: how much of that vorticity is
  streamwise, and whether the updraft gradient is oriented to lift the vortex lines.

Reads only; no solver state is modified.  Thresholds are calibrated on this model's own
runs and are NOT a climatological skill claim -- see :data:`GATES`.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .vorticity_budget import horizontal_vorticity, streamwise_crosswise, tilting_efficiency
from .soundings import bunkers_storm_motion, bulk_shear, storm_relative_helicity

# Gate values taken from the measured A-K sequence in docs/TORNADOGENESIS_FINDINGS.md,
# not from an observational climatology.  They separate this model's own non-tornadic
# runs from its one surface-connected case; the sample is small and single-model.
GATES = {
    'omega_h_0_1km_s': 5.0e-3,      # abundance; attempt H measured 1.1e-2 and still failed
    'streamwise_fraction': 0.50,    # attempt I: 0.40 (no vortex) -> 0.64 (vortex)
    'alignment_updraft': 0.15,      # forced runs stuck at 0.09; the free case peaked at 0.20
    'lcl_m': 1200.0,                # low cloud base; classical, kept loose
    'srh_0_1km_m2_s2': 100.0,       # classical supercell floor
}


@dataclass
class PrecursorScore:
    """Components, gate outcomes, and an explicit statement of what is missing."""
    values: dict = field(default_factory=dict)
    gates: dict = field(default_factory=dict)
    kind: str = 'environment'
    note: str = ''

    @property
    def passed(self):
        return [k for k, v in self.gates.items() if v is True]

    @property
    def failed(self):
        return [k for k, v in self.gates.items() if v is False]

    def summary(self):
        lines = ['%s precursor score' % self.kind,
                 '  gates passed: %s' % (', '.join(self.passed) or 'none'),
                 '  gates failed: %s' % (', '.join(self.failed) or 'none')]
        for k, v in self.values.items():
            if isinstance(v, float):
                lines.append('  %-28s %+.5g' % (k, v))
        if self.note:
            lines.append('  NOTE: %s' % self.note)
        return '\n'.join(lines)


def _layer(z, lo, hi):
    sel = (z >= lo) & (z <= hi)
    return sel if sel.sum() >= 2 else (z <= max(hi, z[1]))


def environment_precursor(base, storm_motion=None, layers=((0.0, 1000.0), (0.0, 3000.0))):
    """Score a 1-D sounding: the supply side of the verbal precursor.

    ``base`` is a :class:`meteorological_flow.base_state.BaseState` (``zc``, ``u0``, ``v0``,
    ``theta0``, ``qv0``, ``p0``, ``T0``).  ``storm_motion`` defaults to the Bunkers
    right-mover, which is what makes the helicity storm-relative rather than ground-relative.

    The streamwise fraction here uses the ENVIRONMENTAL vorticity
    ``omega_h = (-dv/dz, du/dz)`` -- the dw terms do not exist in a sounding.  It therefore
    measures how the hodograph is oriented relative to the storm-relative flow, which is the
    quantity the verbal rule is really about, and it is an upper bound on what the storm can
    ingest, not what it actually tilts.
    """
    from meteorological_flow.base_state import sounding_diagnostics

    z = np.asarray(base.zc, dtype=float)
    u = np.asarray(base.u0, dtype=float)
    v = np.asarray(base.v0, dtype=float)
    cx, cy = storm_motion if storm_motion is not None else bunkers_storm_motion(base)

    dudz = np.gradient(u, z)
    dvdz = np.gradient(v, z)
    xi, eta = -dvdz, dudz                      # environmental (shear-only) horizontal vorticity
    omh = np.hypot(xi, eta)

    ur, vr = u - cx, v - cy
    speed = np.hypot(ur, vr) + 1e-12
    streamwise = (xi * ur + eta * vr) / speed
    crosswise = (-xi * vr + eta * ur) / speed

    values = {'storm_motion_x_m_s': float(cx), 'storm_motion_y_m_s': float(cy)}
    for lo, hi in layers:
        sel = _layer(z, lo, hi)
        tag = '%g_%gkm' % (lo / 1000.0, hi / 1000.0)
        mag = float(np.mean(omh[sel]))
        values['omega_h_' + tag + '_s'] = mag
        values['streamwise_mean_' + tag + '_s'] = float(np.mean(streamwise[sel]))
        values['crosswise_mean_' + tag + '_s'] = float(np.mean(crosswise[sel]))
        values['streamwise_fraction_' + tag] = float(
            np.mean(np.abs(streamwise[sel])) / (np.mean(omh[sel]) + 1e-20))
        values['srh_' + tag + '_m2_s2'] = float(
            storm_relative_helicity(base, z_top=hi, storm_motion=(cx, cy)))
    values['bulk_shear_0_1km_m_s'] = float(bulk_shear(base, 0.0, 1000.0))
    values['bulk_shear_0_6km_m_s'] = float(bulk_shear(base, 0.0, 6000.0))

    diag = sounding_diagnostics(base)
    for key in ('CAPE_J_kg', 'CIN_J_kg', 'LCL_m', 'LFC_m'):
        if diag.get(key) is not None:
            values[key] = float(diag[key])

    gates = {
        'omega_h_0_1km': values.get('omega_h_0_1km_s', 0.0) >= GATES['omega_h_0_1km_s'],
        'streamwise_fraction': values.get('streamwise_fraction_0_1km', 0.0)
                               >= GATES['streamwise_fraction'],
        'srh_0_1km': values.get('srh_0_1km_m2_s2', 0.0) >= GATES['srh_0_1km_m2_s2'],
    }
    if 'LCL_m' in values:
        gates['lcl'] = values['LCL_m'] <= GATES['lcl_m']

    return PrecursorScore(values=values, gates=gates, kind='environment',
                          note='Supply side only. A sounding contains no updraft, so the '
                               'alignment cos(theta) that discriminated in this model cannot '
                               'be evaluated here. Passing every gate is necessary, not '
                               'sufficient: attempt E raised SRH 254 -> 648 for ~0% change.')


def state_precursor(uc, vc, wc, grid, storm_motion=(0.0, 0.0), z_lo=0.0, z_hi=1000.0,
                    w_updraft=1.0):
    """Score a 3-D state: the geometry the sounding cannot contain.

    Adds the discriminating quantity -- ``alignment = cos(theta)`` between the horizontal
    vorticity and the horizontal updraft gradient, from
    :func:`vorticity_budget.tilting_efficiency`.  Because a domain mean washes it out (the
    tilting that matters happens inside the updraft), the alignment and streamwise fraction
    are reported BOTH unconditionally over the layer and conditioned on ``w >= w_updraft``.

    ``uc, vc, wc`` are cell-centred velocity components.  Set ``storm_motion`` to the storm
    motion (Bunkers, or the tracked vortex translation) or the streamwise split is measured
    in the ground-relative frame, which is not the frame the physics cares about.
    """
    xp = grid.xp
    z = np.asarray(_to_cpu(grid.zc), dtype=float)
    sel = (z >= z_lo) & (z <= z_hi)
    if not sel.any():
        raise ValueError('no model levels inside %g-%g m' % (z_lo, z_hi))

    eff = tilting_efficiency(uc, vc, wc, grid)
    xi, eta = horizontal_vorticity(uc, vc, wc, grid)
    sw, cw = streamwise_crosswise(uc, vc, xi, eta, grid, storm_motion=storm_motion)

    layer = lambda a: _to_cpu(a)[:, :, sel]
    omh = layer(eff['omega_h'])
    gw = layer(eff['grad_h_w'])
    tilt = layer(eff['tilting'])
    align = layer(eff['alignment'])
    sww = layer(sw)
    w = layer(wc)
    up = w >= w_updraft

    values = {
        'omega_h_mean_s': float(np.mean(omh)),
        'grad_h_w_mean_s': float(np.mean(gw)),
        'tilting_mean_s2': float(np.mean(tilt)),
        'alignment_mean': float(np.mean(align)),
        'streamwise_fraction': float(np.mean(np.abs(sww)) / (np.mean(omh) + 1e-20)),
        'updraft_fraction': float(np.mean(up)),
    }
    unweighted = False
    if up.any():
        # Weight by the updraft so the number describes the air actually being tilted.
        # With a threshold at or below zero the selection can contain no rising air at
        # all, and the weights then sum to zero; fall back to a plain mean and say so
        # rather than dividing by it.
        wgt = np.clip(w[up], 0.0, None)
        if not np.any(wgt > 0):
            wgt = np.ones_like(wgt)
            unweighted = True
        values['alignment_updraft'] = float(np.average(align[up], weights=wgt))
        values['omega_h_updraft_s'] = float(np.average(omh[up], weights=wgt))
        values['tilting_updraft_s2'] = float(np.average(tilt[up], weights=wgt))
        values['streamwise_fraction_updraft'] = float(
            np.average(np.abs(sww[up]), weights=wgt) / (values['omega_h_updraft_s'] + 1e-20))
    else:
        values['alignment_updraft'] = float('nan')

    gates = {
        'omega_h': values['omega_h_mean_s'] >= GATES['omega_h_0_1km_s'],
        'streamwise_fraction': values.get('streamwise_fraction_updraft',
                                          values['streamwise_fraction'])
                               >= GATES['streamwise_fraction'],
    }
    a = values['alignment_updraft']
    gates['alignment'] = bool(a >= GATES['alignment_updraft']) if np.isfinite(a) else False

    note = ('Alignment is the quantity that separated this model\'s cases; a negative value '
            'means tilting produces ANTIcyclonic vorticity.')
    if np.isfinite(a) and a < 0:
        note += ' Alignment is NEGATIVE here.'
    if not up.any():
        note += ' No cell in the layer reaches w >= %g m/s, so the conditional scores are ' \
                'undefined and the unconditional ones describe air that is not rising.' % w_updraft
    elif unweighted:
        note += ' No rising air in the selection, so the conditional scores are plain means, ' \
                'not updraft-weighted.'
    return PrecursorScore(values=values, gates=gates, kind='state', note=note)


def _to_cpu(a):
    return a.get() if hasattr(a, 'get') else np.asarray(a)


__all__ = ['GATES', 'PrecursorScore', 'environment_precursor', 'state_precursor']
