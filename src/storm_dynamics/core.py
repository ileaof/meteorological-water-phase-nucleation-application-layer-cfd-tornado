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

from . import forcing as frc
from . import momentum as mom
from . import rotation as rot
from . import surface_drag as sfc
from . import turbulence as les
from .config import StormConfig
from .soundings import build_sounding


def _pressure_method(grid: Grid) -> str:
    """Pick the pressure Poisson solver.  Stretched grids (all storm grids) use the exact
    direct (splu) solve; large *uniform* grids fall back to CG.

    KNOWN LIMITATION (memory / fine AMR nests): the direct solve stores a full sparse LU
    factorisation, which OOMs at ~48³.  Switching to the existing CG is **not** a fix —
    Jacobi-preconditioned CG does *not* converge on the stretched anelastic operator
    (verified: NaN on the periodic pure-Neumann system, and residual ~8 after 5000 iters on
    the wall system; `test_pressure_cg_does_not_converge_on_stretched_operator`).  The real
    low-memory fix is a proper solver — an **FFT/DST-in-(x,y) + tridiagonal-in-z** direct
    solve (exploiting the stretched-z-only structure: low memory, exact) or an
    **ILU/multigrid-preconditioned CG** — see docs/ROADMAP.md §2b/§3f.  Until then keep fine
    nests small (`render_tornado_3d.py --sub-half-frac`)."""
    n = grid.nx * grid.ny * grid.nz
    if getattr(grid, "stretched", False):
        return "direct"
    return "direct" if n <= 64_000 else "cg"


# Above this many cells the direct sparse-LU factorisation OOMs (verified ~48³); route such
# grids through the low-memory solvers (ROADMAP §3f).  Smaller grids keep the exact direct
# solve, so every existing small-grid run/test is byte-for-byte unchanged.
_LOWMEM_N = 64_000


def separable(grid: Grid) -> bool:
    """Whether the pressure Poisson **separates** — uniform horizontal spacing, coefficients
    homogeneous in x,y, separable homogeneous BCs — so the exact FFT/DCT-in-(x,y) +
    tridiagonal-in-z solver applies (``pressure_fft``).  The current storm and *all* its AMR
    nests are separable (uniform dx,dy; anelastic ρ0=ρ0(z)); terrain-following coordinates or
    x,y-varying reference states (not yet in the model) would make this ``False`` and select
    the general Jacobi-CG fallback (``pressure_iterative``).  Never a silent default: the
    caller routes on this predicate."""
    if getattr(grid, "terrain", None) is not None:
        return False                          # orography / immersed boundary -> non-separable
    if getattr(grid, "horizontal_stretch", False):
        return False                          # x,y grid stretching -> non-separable
    return True


def _project_anelastic_lowmem(st, grid: Grid, rho0_c, rho0_wface):
    """Low-memory anelastic projection routed by grid structure: **FFT+tridiag** when the grid
    is :func:`separable` (the fast path — the case today), else the general **Jacobi-CG**
    fallback.  Operates on host arrays (device round-trip on GPU) and writes the
    divergence-free velocity back into ``st``.  Returns ``max|div(ρ0 u)|`` after projection.

    Physically equivalent to :meth:`PressureSolver.project_anelastic` (the divergence-free
    projection is unique for the given BCs), but with no stored LU factorisation, so a fine
    nest that OOMs the direct solve fits.  On a **GPU** backend the separable FFT path runs
    entirely on the device (cuFFT / cupyx DCT + batched Thomas) with **no host round-trip**;
    the non-separable Jacobi-CG fallback is scipy/host-only, so it round-trips."""
    to_cpu = grid.backend.to_cpu
    xp = grid.xp
    zc = np.asarray(to_cpu(grid.zc), float); zf = np.asarray(to_cpu(grid.zf), float)
    dzc = zf[1:] - zf[:-1]; dzf = np.diff(zc)                 # tiny host coefficients
    dx, dy = float(grid.dx), float(grid.dy)
    periodic_h = bool(getattr(grid, "periodic", True))
    sep = separable(grid)
    is_gpu = type(st.u).__module__.split(".", 1)[0] == "cupy"
    if sep and is_gpu:                                        # GPU fast path: solve on-device
        from .pressure_fft import project_anelastic_fft as _proj
        return float(_proj(st.u, st.v, st.w, rho0_c, rho0_wface, dx, dy, dzc, dzf,
                           periodic_h=periodic_h))            # modifies st.{u,v,w} in place
    u = np.ascontiguousarray(to_cpu(st.u), float)             # host path (CPU, or the
    v = np.ascontiguousarray(to_cpu(st.v), float)             # scipy Jacobi-CG fallback)
    w = np.ascontiguousarray(to_cpu(st.w), float)
    rc = np.asarray(to_cpu(rho0_c), float); rw = np.asarray(to_cpu(rho0_wface), float)
    if sep:
        from .pressure_fft import project_anelastic_fft as _proj
        res = _proj(u, v, w, rc, rw, dx, dy, dzc, dzf, periodic_h=periodic_h)
    else:
        from .pressure_iterative import project_anelastic_iterative as _proj
        res = _proj(u, v, w, rc, rw, dx, dy, dzc, dzf, periodic_h=periodic_h)
    st.u = xp.asarray(u); st.v = xp.asarray(v); st.w = xp.asarray(w)
    return float(res)


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
        # large anelastic grids: route the projection through the low-memory solver
        # (ROADMAP §3f) -- the direct LU OOMs past ~48³.  Small grids keep the exact
        # direct solve (unchanged).  Opt-out/force via scfg if ever needed.
        self._lowmem_pressure = (self.dynamics == "anelastic"
                                 and g.nx * g.ny * g.nz > _LOWMEM_N)
        # When the low-memory FFT/tridiag path handles projection, the direct
        # PressureSolver (host splu of an NxN Laplacian, assembled in a Python index
        # loop) is never used -- skip building it: it costs minutes and GBs at ~1e6
        # cells (a stretched grid forces the direct solve on the GPU).  Small grids
        # keep the exact direct solve, byte-identical to before.
        self.pressure = None if self._lowmem_pressure else PressureSolver(g, method=_pressure_method(g))
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
    # The step is split into phases so the composite (parent+nest) projection can
    # run once over BOTH levels (storm_dynamics.nesting.run_concurrent_nest) in
    # place of the two per-level solves: _predictor -> _project -> _transport.
    def _step(self, dt: float) -> None:
        Km = self._predictor(dt)
        self._project(dt)
        self._transport(dt, Km)

    def _predictor(self, dt: float) -> "object":
        "BCs + diagnose, then the momentum predictor (LES, advection, forces)."
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
        if self.dyn.les.model == "tke15":       # prognostic Deardorff TKE-1.5 (evolve e(t))
            Km, self._tke = les.deardorff_tke_step(st, g, self.dyn.les,
                                                   getattr(self, "_tke", None), dt, self.theta0_field)
        else:
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
        # (f) sustained low-level mesoscale-ascent forcing (dryline/convergence proxy):
        #     keeps lifting parcels through a real CIN cap so a supercell can establish
        #     instead of a one-shot bubble that decays.  Opt-in; off by default.
        fcfg = getattr(self.dyn, "forcing", None)
        if fcfg is not None and getattr(fcfg, "enabled", False):
            frc.apply_meso_forcing(st, g, fcfg, float(self.t), dt,
                                   center=getattr(fcfg, "center", None), xp=xp)
        # extreme numerical guard ONLY (documented; not a physical cap)
        vg = self.dyn.v_guard
        xp.clip(st.u, -vg, vg, out=st.u); xp.clip(st.v, -vg, vg, out=st.v)
        xp.clip(st.w, -vg, vg, out=st.w)
        bc.apply_velocity_bcs(st, g, cfg)
        return Km

    def _project(self, dt: float) -> None:
        "Step 3: this level's own anelastic projection (skipped in composite mode)."
        st = self.state
        if self.dynamics == "anelastic" and getattr(self, "_lowmem_pressure", False):
            res, it = _project_anelastic_lowmem(st, self.grid, self.rho0_c, self.rho0_wface), 0
        elif self.dynamics == "anelastic":
            res, it = self.pressure.project_anelastic(st, dt, self.rho0_c, self.rho0_wface)
        else:
            res, it = self.pressure.project(st, dt, self.rho0)
        self._last_res, self._last_iters = res, it
        bc.apply_velocity_bcs(st, self.grid, self.cfg)

    def _transport(self, dt: float, Km: "object") -> None:
        cfg = self.cfg
        g = self.grid
        xp = g.xp
        st = self.state
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
    def run(self, progress=None, record_interval=None, capture_frames=False) -> dict:
        """Run the storm.  ``capture_frames`` stores lightweight 2-D vorticity /
        vertical-velocity slices at each record interval into ``self.frames`` for
        animation (see :func:`storm_dynamics.plotting.animate_rotation`)."""
        cfg = self.cfg
        self._Km = None
        self._capture_frames = bool(capture_frames)
        self.frames = []
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

    def _capture_frame(self) -> None:
        """Append a 2-D snapshot (mid-level ζ, w, near-surface ζ) for animation."""
        g = self.grid
        to_cpu = g.backend.to_cpu
        _, _, zeta = rot.vorticity_3d(self.state, g)
        _, _, wc = rot._centered_velocity(self.state, g)
        z = to_cpu(g.zc)
        import numpy as _np
        kmid = int(_np.argmin(_np.abs(z - 4000.0)))
        knear = int(_np.argmin(_np.abs(z - 500.0)))
        self.frames.append({
            "t": self.t,
            "zeta_mid": to_cpu(zeta[:, :, kmid]),
            "w_mid": to_cpu(wc[:, :, kmid]),
            "zeta_near": to_cpu(zeta[:, :, knear]),
        })

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
        if getattr(self, "_capture_frames", False):
            self._capture_frame()

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
