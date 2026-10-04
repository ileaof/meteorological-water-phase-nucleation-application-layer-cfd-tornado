"""Passive per-stage momentum flux accounting, independent of endpoint residuals.

Fluxes are reevaluated by the native routines on their unmodified input state.
The reconstructed increments are checked against the actual solver/tracer changes.
No routine here advances the physical state.
"""
from types import SimpleNamespace

import numpy as np

from .momentum import _u_tendency, _v_tendency
from .vorticity_provenance import (
    VorticityProvenanceTracer, _frozen_u_tendency, _frozen_v_tendency,
)

TERMS = ('total_les', 'total_muscl', 'total_flux_x', 'total_flux_y', 'total_flux_z',
         'total_roundoff', 'total_other', 'les_injection', 'les_muscl',
         'les_flux_x', 'les_flux_y', 'les_flux_z', 'les_roundoff')


def horizontal_curl(uv, grid):
    u, v = uv
    uc, vc = (u[:-1]+u[1:])*0.5, (v[:, :-1]+v[:, 1:])*0.5
    return grid._central_x(vc)-grid._central_y(uc)


def transpose_derivative(mask, spacing, axis):
    m = np.moveaxis(mask, axis, 0)
    result = np.zeros_like(m, dtype=float)
    result[0] -= m[0]/spacing
    result[1] += m[0]/spacing
    result[-2] -= m[-1]/spacing
    result[-1] += m[-1]/spacing
    result[:-2] -= m[1:-1]/(2*spacing)
    result[2:] += m[1:-1]/(2*spacing)
    return np.moveaxis(result, 0, axis)


def flux_increments(flux, grid, dt):
    """Reconstruct directional native u/v increments, including native edge rules."""
    xp = grid.xp
    answer = {axis: [] for axis in 'xyz'}
    for component, shape, normal in (('u', grid.u_shape, 0), ('v', grid.v_shape, 1)):
        for axis, name, spacing in ((0, 'x', grid.dx), (1, 'y', grid.dy)):
            f = xp.moveaxis(flux[f'F{name}_{component}'], axis, 0)
            out = xp.zeros(shape)
            target = xp.moveaxis(out, axis, 0)
            if axis == normal:
                target[1:-1] = -(f[1:]-f[:-1])/spacing
                if grid.periodic:
                    target[0] = -(f[0]-f[-1])/spacing
                    target[-1] = target[0]
            elif grid.periodic:
                target[:] = -(xp.roll(f, -1, axis=0)-f)/spacing
            else:
                target[1:-1] = -(f[2:-1]-f[1:-2])/spacing
            answer[name].append(dt*out)
        f = flux[f'Fz_{component}']
        dz = grid.dz_c[None, None, :] if grid.stretched else grid.dz
        answer['z'].append(-dt*(f[:, :, 1:]-f[:, :, :-1])/dz)
    return answer


class DirectMusclCapture:
    def __init__(self, sim, masks, nk=17):
        self.grid = sim.grid
        self.xp = sim.grid.xp
        self.to = sim.grid.backend.to_cpu
        self.nk = nk
        self.tracer = VorticityProvenanceTracer(sim, record_stage_metrics=False)
        self.label_index = self.tracer.source_index['les']
        self.gates = {}
        self.stage_order = []
        self.set_masks(masks)

    def label_uv(self):
        return tuple(a[self.label_index] for a in self.tracer.velocity_sources[:2])

    def q(self):
        return tuple(horizontal_curl(uv, self.grid)[:, :, :self.nk].copy()
                     for uv in (self.tracer.previous_velocity[:2], self.label_uv()))

    def set_masks(self, masks):
        self.masks_cpu = np.asarray(masks, dtype=float)
        self.masks = self.xp.asarray(self.masks_cpu)
        self.wx = self.xp.asarray(np.stack([transpose_derivative(m, self.grid.dx, 0) for m in masks]))
        self.wy = self.xp.asarray(np.stack([transpose_derivative(m, self.grid.dy, 1) for m in masks]))
        self.accum = {}
        self.block_start = self.q()
        self.records = []

    def check(self, name, actual, predicted, floor, units):
        xp = self.xp
        delta = actual-predicted
        err = float(xp.sqrt(xp.mean(delta*delta)).item())
        scale = max(float(xp.sqrt(xp.mean(actual*actual)).item()),
                    float(xp.sqrt(xp.mean(predicted*predicted)).item()), floor)
        absolute = float(xp.max(xp.abs(delta)).item())
        row = self.gates.setdefault(name, dict(relative_max=0., absolute_max=0., floor=floor, units=units))
        row['relative_max'] = max(row['relative_max'], err/scale)
        row['absolute_max'] = max(row['absolute_max'], absolute)
        if not np.isfinite(err/scale) or err/scale > 1e-10:
            raise RuntimeError(f'{name}: relative={err/scale:.6g}, absolute={absolute:.6g} {units}')

    def project(self, field, weights=None):
        return self.grid.dx*self.grid.dy*self.xp.tensordot(
            self.masks if weights is None else weights, field, axes=([1, 2], [0, 1]))

    def stats(self, uv):
        u, v = [a[:, :, :self.nk] for a in uv]
        zeta = horizontal_curl((u, v), self.grid)
        signed = self.project(zeta)
        absolute = self.project(self.xp.abs(zeta))
        boundary = self.project((v[:, :-1]+v[:, 1:])*0.5, self.wx)
        boundary -= self.project((u[:-1]+u[1:])*0.5, self.wy)
        self.check('curl_adjoint', signed, boundary, 1., 'm2/s')
        return self.to(self.xp.stack((signed, absolute, boundary), axis=-1))

    def measured_flux(self, before, dt, label=None):
        flux = {}
        if label is None:
            state = SimpleNamespace(**dict(zip(('u', 'v', 'w'), before)))
            du = _u_tendency(state, self.grid, self.tracer.momentum_order, self.grid.periodic, flux)
            dv = _v_tendency(state, self.grid, self.tracer.momentum_order, self.grid.periodic, flux)
        else:
            du = _frozen_u_tendency(label, before, self.grid, self.tracer.momentum_order, self.grid.periodic, flux)
            dv = _frozen_v_tendency(label, before, self.grid, self.tracer.momentum_order, self.grid.periodic, flux)
        for key in ('Fz_u', 'Fz_v'):
            if bool(self.xp.any(flux[key][:, :, (0, -1)] != 0)):
                raise RuntimeError('Unexpected nonzero vertical wall flux')
        return flux_increments(flux, self.grid, dt), (dt*du, dt*dv)

    def mark(self, sim, name, dt):
        before = self.tracer.previous_velocity
        if name == 'begin':
            self.step_start = self.q()
            self.step_terms = {}
            self.others = [self.xp.zeros_like(a[:, :, :self.nk]) for a in before[:2]]
            self.stage_order = []
            self.tracer.mark(sim, name, dt)
            return
        self.stage_order.append(name)
        label_before = tuple(a.copy() for a in self.label_uv()) if name in ('les', 'advection') else None
        directions = {}
        if name == 'advection':
            directions['total'] = self.measured_flux(before, dt)
            label = tuple(a[self.label_index] for a in self.tracer.velocity_sources)
            directions['les'] = self.measured_flux(before, dt, label)
        self.tracer.mark(sim, name, dt)
        actual = tuple(a-b for a, b in zip(self.tracer.previous_velocity[:2], before[:2]))
        if name in ('les', 'advection'):
            labelled = tuple(a-b for a, b in zip(self.label_uv(), label_before))
            if name == 'les':
                self.step_terms['total_les'] = actual
                self.step_terms['les_injection'] = labelled
                for i in range(2):
                    self.check('les_injection_'+str(i), labelled[i], actual[i], 1e-12, 'm/s')
            else:
                for prefix, measured in (('total', actual), ('les', labelled)):
                    self.step_terms[prefix+'_muscl'] = measured
                    directional, native_increment = directions[prefix]
                    reconstructed = tuple(sum(directional[d][i] for d in 'xyz') for i in range(2))
                    before_update = before[:2] if prefix == 'total' else label_before
                    for i in range(2):
                        self.check(prefix+'_flux_reconstruction_'+str(i), native_increment[i], reconstructed[i], 1e-12, 'm/s')
                        rounded_update = (before_update[i]+native_increment[i])-before_update[i]
                        self.check(prefix+'_actual_update_'+str(i), measured[i], rounded_update, 1e-12, 'm/s')
                    for d in 'xyz':
                        self.step_terms[prefix+'_flux_'+d] = directional[d]
                    self.step_terms[prefix+'_roundoff'] = tuple(a-b for a,b in zip(measured,reconstructed))
        else:
            for dest, delta in zip(self.others, actual):
                dest += delta[:, :, :self.nk]
        if name == 'transport_microphysics_bcs':
            self.step_terms['total_other'] = self.others
            for quantity, names in enumerate((('total_les','total_muscl','total_other'),
                                               ('les_injection','les_muscl'))):
                predicted = sum(horizontal_curl(self.step_terms[n], self.grid)[:, :, :self.nk] for n in names)
                self.check(('total' if quantity==0 else 'les')+'_step_closure',
                           self.q()[quantity]-self.step_start[quantity], predicted, 1e-12, '1/s')
            records = []
            for term in TERMS:
                uv = tuple(a[:, :, :self.nk] for a in self.step_terms[term])
                records.append(self.stats(uv))
                if term not in self.accum:
                    self.accum[term] = [a.copy() for a in uv]
                else:
                    for dest, delta in zip(self.accum[term], uv):
                        dest += delta
            self.last_record = np.stack(records)

    def finish_block(self, mask_pairs):
        end = self.q()
        statistics = np.stack([self.stats(self.accum[t]) for t in TERMS])
        balances = []
        for m, (ma, mb, symmetric) in enumerate(mask_pairs):
            for quantity, names in enumerate((('total_les','total_muscl','total_other'),
                                               ('les_injection','les_muscl'))):
                qa, qb = self.block_start[quantity], end[quantity]
                a, b = self.xp.asarray(ma), self.xp.asarray(mb)
                first = self.grid.dx*self.grid.dy*self.xp.sum(a[:, :, None]*qa, axis=(0,1))
                last = self.grid.dx*self.grid.dy*self.xp.sum(b[:, :, None]*qb, axis=(0,1))
                qmotion = (qa+qb)*.5 if symmetric else qa
                motion = self.grid.dx*self.grid.dy*self.xp.sum((b-a)[:, :, None]*qmotion, axis=(0,1))
                sums = self.xp.asarray(sum(statistics[TERMS.index(n), m, :, 0] for n in names))
                self.check('mask_balance_'+str(quantity), last-first, sums+motion, 1., 'm2/s')
                balances.append(self.to(self.xp.stack((first,last,motion),axis=-1)))
        return statistics, np.stack(balances).reshape(len(mask_pairs),2,self.nk,3)
