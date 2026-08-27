"""StormSimulation -- the rotating deep-convection time loop.

A fork of :class:`meteorological_flow.simulation.Simulation` that swaps in the
rotational dynamical core while **reusing** the repo's grid, anelastic Chorin
projection, conservative scalar transport, moist buoyancy, bulk microphysics
(evaporative cold pool = baroclinic vorticity source) and the validated
nucleation kernel unchanged.

Per step (anelastic projection method, explicit predictor):

1. boundary conditions + diagnose thermodynamics;
2. momentum predictor -- the pieces the demonstration core lacks:
   * **conservative flux-form momentum advection** (item 1: tilting + stretching),
   * moist buoyancy on ``w`` (reused),
   * **f-plane Coriolis** on the perturbation wind (item 2),
   * **LES (Smagorinsky) subgrid diffusion** instead of Rayleigh drag + clip
     (item 3), applied to the perturbation so the environmental shear persists,
   * **surface bulk-drag** on the lowest level (item 4);
   (no Rayleigh drag, no velocity clip -- only an extreme documented guard);
3. anelastic projection -> div(rho0 u) ~ 0 (reused);
4. conservative scalar transport (theta, q_v, hydrometeors) + LES scalar diffusion;
5. two-way microphysics: growth + latent heat + **evaporative cold pool** +
   sedimentation (reused);
6. rotation diagnostics (zeta, updraft helicity, near-surface / mid-level trackers).
"""
from __future__ import annotations

import os
import time as _time

import numpy as np

from meteorological_flow import advection as adv
from meteorological_flow import boundary_conditions as bc
from meteorological_flow import buoyancy as buo
from meteorological_flow import diagnostics as diag
from meteorological_flow import thermodynamics as th
from meteorological_flow.base_state import warm_bubble
from meteorological_flow.grid import Grid
from meteorological_flow.pressure_solver import PressureSolver
from meteorological_flow.state import FlowState

from . import momentum as mom
from . import rotation as rot
from . import surface_drag as sfc
from . import turbulence as les
from .config import StormConfig
from .soundings import build_sounding


def _pressure_method(grid: Grid) -> str:
    n = grid.nx * grid.ny * grid.nz
    if getattr(grid, "stretched", False):
        return "direct"
    return "direct" if n <= 64_000 else "cg"


class StormSimulation:
    def __init__(self, scfg: StormConfig, base=None, backend=None):
        self.scfg = scfg
        self.cfg = scfg.sim                       # reused meteorological_flow config
        self.dyn = scfg.dyn
        cfg = self.cfg
        # compute backend.  An explicit ``backend=`` (e.g. from the example's
        # --device flag) always wins; otherwise honour cfg.performance.device
        # (build_storm_config defaults it to "cpu", preserving the prior
        # CPU-only behaviour).  "gpu" fails loudly; "auto" falls back to CPU.
        from meteorological_flow.backend import get_backend
        self.backend = (backend if backend is not None
                        else get_backend(getattr(cfg.performance, "device", "cpu")))
        periodic = getattr(cfg.boundaries, "x_west", None) == "periodic"
        self.grid = Grid(nx=cfg.grid.nx, ny=cfg.grid.ny, nz=cfg.grid.nz,
                         Lx=cfg.domain.Lx, Ly=cfg.domain.Ly, Lz=cfg.domain.Lz,
                         z_stretch=getattr(cfg.grid, "z_stretch", 1.0),
                         periodic=periodic, backend=self.backend)
        g = self.grid
        xp = g.xp
        # curved-hodograph base state (item 5); reuse WK thermodynamics.
        self.base = base if base is not None else build_sounding(g, self.dyn.hodograph)
        self.theta0_field = self.base.field(self.base.theta0, g.center_shape, xp=xp)
        self.qv0_field = self.base.field(self.base.qv0, g.center_shape, xp=xp)
        self.f = self.dyn.f_value()
        # environmental wind on the faces (for the initial condition + Coriolis)
        self._u0_face = xp.broadcast_to(xp.asarray(self.base.u0, float)[None, None, :],
                                        g.u_shape).copy()
        self._v0_face = xp.broadcast_to(xp.asarray(self.base.v0, float)[None, None, :],
                                        g.v_shape).copy()
        self.state = self._initial_state()
        self.state.diagnose(cfg)
        self.T_ref = float(self.state.T.mean()) if cfg.physics.T_ref is None else cfg.physics.T_ref
        self.qv_ref = float(self.state.qv.mean())
        self.rho0 = cfg.physics.P0 / (th.R_d * self.T_ref)
        # anelastic reference density profiles
        self.dynamics = getattr(cfg.physics, "dynamics", "anelastic")
        rho0_prof = np.asarray(self.base.rho0, dtype=float)
        rho0_wface_h = np.interp(g.backend.to_cpu(g.zf), g.backend.to_cpu(g.zc), rho0_prof)
        self.rho0_c = xp.asarray(rho0_prof)
        self.rho0_wface = xp.asarray(rho0_wface_h)
        if self.dynamics == "anelastic":
            self._transport_rho_c = self.rho0_c
            self._transport_rho_wf = self.rho0_wface
            self.rho_ref = self.rho0_c
        else:
            self._transport_rho_c = xp.ones(g.nz)
            self._transport_rho_wf = xp.ones(g.nz + 1)
            self.rho_ref = self.rho0
        self.pressure = PressureSolver(g, method=_pressure_method(g))
        # two-way microphysics (cold pool) -- reused unchanged
        from meteorological_flow.microphysics_coupling import MicrophysicsCoupler
        self.coupler = MicrophysicsCoupler()
        # optional: couple the validated nucleation kernel as the embryo source
        # (eq39 pathway), exactly as meteorological_flow.Simulation does.
        self.couple_nucleation = bool(getattr(cfg.nucleation, "couple_kernel", False))
        self.adapter = None
        self.lookup = None
        if self.couple_nucleation:
            from meteorological_flow.nucleation_adapter import NucleationAdapter
            self.adapter = NucleationAdapter(cfg.nucleation)
            if cfg.nucleation.method == "lookup":
                from meteorological_flow.simulation import _build_lookup
                self.lookup, _ = _build_lookup(cfg, cfg.output.outdir, self.adapter)
                self.adapter.set_lookup(self.lookup)
        self.tracker = rot.TornadogenesisTracker()
        self.history = []
        self.step = 0
        self.t = 0.0
        self._last_res = 0.0
        self._last_iters = 0
        self._t0 = _time.perf_counter()

    # ---- initial condition: stratified base + environmental wind + warm bubble ----
    def _initial_state(self) -> FlowState:
        g = self.grid
        xp = g.xp
        cfg = self.cfg
        st = FlowState.zeros(g)
        dth, _ = warm_bubble(g, dtheta=cfg.physics.bubble_dtheta)
        st.theta = self.theta0_field + dth
        st.qv = xp.maximum(self.qv0_field.copy(), 0.0)
        st.p0_field = self.base.field(self.base.p0, g.center_shape, xp=xp)
        st.u[:] = self._u0_face
        st.v[:] = self._v0_face
        bc.apply_velocity_bcs(st, g, cfg)
        bc.apply_scalar_bcs(st, g, cfg, theta0=self.theta0_field, qv0=self.qv0_field)
        st.diagnose(cfg)
        # saturate the bubble core so it condenses on first ascent (honest trigger)
        if cfg.physics.bubble_dtheta > 0.0:
            qsat = th.q_v_from_p_v(th.psat_water(st.T, xp=xp), st.P_total, xp=xp)
            core = dth > 0.5 * cfg.physics.bubble_dtheta
            st.qv = xp.where(core, xp.maximum(st.qv, 0.97 * qsat), st.qv)
            st.diagnose(cfg)
        return st

    # ---- adaptive CFL (advective + diffusive), now including advected momentum ----
    def _dt(self) -> float:
        g = self.grid
        xp = g.xp
        st = self.state
        uc = 0.5 * (st.u[:-1] + st.u[1:]); vc = 0.5 * (st.v[:, :-1] + st.v[:, 1:])
        wc = 0.5 * (st.w[:, :, :-1] + st.w[:, :, 1:])
        umax = float(xp.abs(uc).max()) if uc.size else 0.0
        vmax = float(xp.abs(vc).max()) if vc.size else 0.0
        wmax = float(xp.abs(wc).max()) if wc.size else 0.0
        dzmin = g.dz if not getattr(g, "stretched", False) else float(g.dz_c.min())
        inv_adv = umax / g.dx + vmax / g.dy + wmax / dzmin
        self._inv_adv = inv_adv
        adv_dt = self.cfg.time.cfl / max(1.25 * inv_adv, 1e-12)
        # diffusive limit from the peak LES viscosity (bounded below by the config nu)
        Kmax = float(self._Km.max()) if getattr(self, "_Km", None) is not None else self.cfg.flow.nu
        diff_dt = 0.5 / (max(Kmax, 1e-12) * (1.0/g.dx**2 + 1.0/g.dy**2 + 1.0/dzmin**2))
        dt = min(adv_dt, diff_dt, self.cfg.time.dt_max)
        return max(dt, 1e-4)

    # ---- one anelastic projection step ----
    def _step(self, dt: float) -> None:
        cfg = self.cfg
        g = self.grid
        xp = g.xp
        st = self.state
        # 1. BCs + diagnose
        bc.apply_velocity_bcs(st, g, cfg)
        bc.apply_scalar_bcs(st, g, cfg, theta0=self.theta0_field, qv0=self.qv0_field)
        st.diagnose(cfg)
        # 2. momentum predictor -----------------------------------------------
        # (a) LES subgrid closure: eddy viscosity from the resolved strain,
        #     momentum diffusion applied to the PERTURBATION so the environmental
        #     hodograph persists (analogous to the scalar perturbation-diffusion).
        Km = les.strain_and_viscosity(st, g, self.dyn.les, theta0=self.theta0_field)
        self._Km = Km
        st.u -= self._u0_face; st.v -= self._v0_face
        les.apply_les_momentum(st, g, Km, dt)
        st.u += self._u0_face; st.v += self._v0_face
        # (b) conservative flux-form momentum advection (the enabling term)
        if self.dyn.momentum_advection:
            mom.add_momentum_advection(st, g, dt, order=self.dyn.momentum_order)
        # (c) moist buoyancy on w (reused, perturbation vs base state)
        Bf = buo.buoyancy_w_tendency(st, g, cfg, self.T_ref, self.qv_ref,
                                     theta0=self.theta0_field, qv0=self.qv0_field)
        st.w += dt * Bf
        # (d) f-plane Coriolis on the perturbation wind
        if self.dyn.coriolis:
            from .coriolis import add_coriolis
            add_coriolis(st, g, dt, self.f, self._u0_face, self._v0_face)
        # (e) surface bulk drag on the lowest level
        sfc.apply_surface_drag(st, g, dt, self.dyn.drag)
        # extreme numerical guard ONLY (documented; not a physical cap)
        vg = self.dyn.v_guard
        xp.clip(st.u, -vg, vg, out=st.u); xp.clip(st.v, -vg, vg, out=st.v)
        xp.clip(st.w, -vg, vg, out=st.w)
        # 3. anelastic projection -> div(rho0 u) ~ 0
        bc.apply_velocity_bcs(st, g, cfg)
        if self.dynamics == "anelastic":
            res, it = self.pressure.project_anelastic(st, dt, self.rho0_c, self.rho0_wface)
        else:
            res, it = self.pressure.project(st, dt, self.rho0)
        self._last_res, self._last_iters = res, it
        bc.apply_velocity_bcs(st, g, cfg)
        # 4. conservative scalar transport with the divergence-free face velocity
        trc, trwf = self._transport_rho_c, self._transport_rho_wf
        _adv = lambda fld: adv.advect_center_massflux(fld, st.u, st.v, st.w, g, dt,
                                                      trc, trwf, order=2)
        st.theta = _adv(st.theta)
        st.qv = xp.maximum(_adv(st.qv), 0.0)
        st.ensure_hydrometeors()
        for nm in ("ql", "qi", "qr", "qs", "qg", "qh"):
            setattr(st, nm, xp.maximum(_adv(getattr(st, nm)), 0.0))
        # LES scalar diffusion (perturbation only for the stratified reference)
        st.theta = les.les_scalar_diffusion(st.theta, Km, g, self.dyn.les, dt,
                                             base=self.theta0_field)
        st.qv = xp.maximum(les.les_scalar_diffusion(st.qv, Km, g, self.dyn.les, dt,
                                                    base=self.qv0_field), 0.0)
        bc.apply_scalar_bcs(st, g, cfg, theta0=self.theta0_field, qv0=self.qv0_field)
        bc.apply_velocity_bcs(st, g, cfg)
        st.diagnose(cfg)
        # 5. two-way microphysics: growth + latent heat + cold pool + sedimentation.
        #    when kernel coupling is on, evaluate the validated 2nd-order rate on the
        #    current (post-transport) state and pass it as the embryo source.
        nf = None
        if self.couple_nucleation and self.adapter is not None:
            nf = self.adapter.evaluate_field(st, dt, g.cell_vol)
            self.last_nf = nf
        self.coupler.apply(st, g, dt, nf=nf)
        self.coupler.sediment(st, g, dt, rho_ref=self._transport_rho_c)
        self.coupler.zero_inflow_hydrometeors(st)
        bc.apply_scalar_bcs(st, g, cfg, theta0=self.theta0_field, qv0=self.qv0_field)
        st.diagnose(cfg)
        # deep-column T guard (documented; only bites at extreme cells)
        Tc = xp.clip(st.T, 180.0, 335.0)
        if not bool(xp.array_equal(Tc, st.T)):
            st.theta = th.theta_from_T(Tc, st.P_total, th.P0_REF, xp=xp)
            st.diagnose(cfg)
        st.t = self.t + dt

    # ---- main loop ----
    def run(self, progress=None, record_interval=None) -> dict:
        cfg = self.cfg
        self._Km = None
        initial = diag.initial_budgets(self.state, self.rho_ref)
        duration = cfg.time.duration
        interval = record_interval or max(1, cfg.output.interval_steps)
        self.tracker.update(self.t, self.state, self.grid)
        self._record(initial)
        while self.t < duration - 1e-9:
            dt = self._dt()
            if self.t + dt > duration:
                dt = duration - self.t
            self._step(dt)
            self.step += 1
            self.t = float(self.state.t)
            if self.step % interval == 0 or self.t >= duration - 1e-9:
                self.tracker.update(self.t, self.state, self.grid)
                self._record(initial)
            if progress and (self.step % max(1, min(interval, 10)) == 0):
                progress(self.t, duration, self.step)
        return self._finalise(initial)

    def _record(self, initial: dict) -> None:
        st = self.state
        bud = diag.conservation_budgets(st, initial, self.rho_ref)
        rep = rot.rotation_report(st, self.grid, base=self.base)
        row = {"time": self.t, "step": self.step,
               "total_water_rel_err": bud["total_water_rel_err"],
               "total_energy_rel_err": bud["total_energy_rel_err"],
               "solver_residual": self._last_res,
               **{k: rep[k] for k in ("zeta_abs_max", "midlevel_mesocyclone",
                                      "near_surface_zeta_max", "updraft_helicity_max",
                                      "w_max", "w_min")}}
        self.history.append(row)

    def _finalise(self, initial: dict) -> dict:
        cfg = self.cfg
        bud = diag.conservation_budgets(self.state, initial, self.rho_ref)
        mres = diag.mass_continuity_residual(
            self.state,
            rho0_c=self.rho0_c if self.dynamics == "anelastic" else None,
            rho0_wface=self.rho0_wface if self.dynamics == "anelastic" else None)
        rep = rot.rotation_report(self.state, self.grid, base=self.base)
        report = {
            "n_steps": self.step, "final_time": self.t,
            "wall_clock_s": _time.perf_counter() - self._t0,
            "dynamics": self.dynamics,
            "f": self.f,
            "backend": {"name": self.backend.name,
                        "fallback_reason": self.backend.fallback_reason},
            "kernel_coupled": self.couple_nucleation,
            "rotation": rep,
            "rotation_peak": self.tracker.peak(),
            "conservation": {
                "total_water_rel_err": bud["total_water_rel_err"],
                "total_energy_rel_err": bud["total_energy_rel_err"],
                "mass_continuity_residual_abs": mres["abs_max"],
                "mass_continuity_residual_norm": mres["normalised"],
            },
            "limitations": [
                "Idealised simulation, NOT operational forecasting: no data "
                "assimilation, no real-event initial/boundary conditions, no "
                "observational verification.",
                "Demonstration-scale grid: the vortex is under-resolved; zeta / "
                "updraft-helicity magnitudes are indicative, not quantitative, "
                "until the grid resolves the ~100 m vortex (M3, not delivered).",
                "f-plane Coriolis (constant f); anelastic deep-column core.",
            ],
        }
        return report


__all__ = ["StormSimulation"]
