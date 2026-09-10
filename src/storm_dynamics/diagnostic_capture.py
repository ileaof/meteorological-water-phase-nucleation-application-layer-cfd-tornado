"""Opt-in, passive fixed-grid capture. Requires the optional h5py package.

Native momentum increments are accumulated over *all* steps in each output
interval. Independent continuous kinematics are sampled before momentum
advection, dt-weighted. Their difference from the actual advective curl is
retained, never silently called diffusion. No arrays owned by the solver are
modified. Nest/composite synchronization needs a separate hierarchy observer.
"""
from __future__ import annotations
import json
import numpy as np


def centered(velocity):
    u,v,w=velocity
    return (u[:-1]+u[1:])/2, (v[:,:-1]+v[:,1:])/2, (w[:,:,:-1]+w[:,:,1:])/2


def grad(a, coordinates):
    return tuple(np.gradient(a,c,axis=i,edge_order=2) for i,c in enumerate(coordinates))


def curl(velocity, coordinates):
    u,v,w=centered(velocity)
    ux,uy,uz=grad(u,coordinates); vx,vy,vz=grad(v,coordinates); wx,wy,wz=grad(w,coordinates)
    return np.array((wy-vz,uz-wx,vx-uy))


def kinematics(velocity, coordinates):
    u,v,w=centered(velocity)
    derivatives=[grad(a,coordinates) for a in (u,v,w)]
    ux,uy,uz=derivatives[0]; vx,vy,vz=derivatives[1]; wx,wy,wz=derivatives[2]
    omega=np.array((wy-vz,uz-wx,vx-uy)); divergence=ux+vy+wz
    stretching=np.array([sum(omega[j]*derivatives[i][j] for j in range(3)) for i in range(3)])
    dilatation=-omega*divergence
    advection=np.array([-sum(vel*d for vel,d in zip((u,v,w),grad(om,coordinates))) for om in omega])
    # Separate the vertical component's axial stretching and tilting.
    return dict(omega=omega,vector_stretching=stretching,vector_dilatation=dilatation,
                vector_advection=advection,stretching=omega[2]*wz,
                tilting=omega[0]*wx+omega[1]*wy,convergence=-(ux+vy))


class DiagnosticCapture:
    def __init__(self, sim, path, metadata, interval=30., top=2000., full_column=True):
        import h5py
        self.file=h5py.File(path,'x')  # Never overwrite a prior run.
        self.file.attrs['schema']='storm-diagnostic-sequence-v1'
        self.file.attrs['metadata']=json.dumps(metadata,default=str)
        self.file.attrs['status']='running'
        self.file.attrs['pressure_semantics']='p: thermodynamic perturbation; p_dyn: projection phi/dt, not coupled to saturation'
        self.file.attrs['increments_semantics']='native velocity changes summed across every solver step; divide by interval duration for mean tendency'
        self.file.attrs['kinematics_semantics']='relative-vorticity terms, independently differenced before advection, dt-weighted integrals'
        self.file.attrs['microphysics_semantics']='anelastic rho0-weighted transfer mass and latent energy, summed over every native step'
        self.interval=interval; self.next_output=sim.t+interval; self.count=0
        self.to=sim.grid.backend.to_cpu
        g=sim.grid; z=np.asarray(self.to(g.zc))
        self.nk=min(len(z),int(np.searchsorted(z,top,side='right'))+3)
        self.coords=(np.asarray(self.to(g.xc)),np.asarray(self.to(g.yc)),z[:self.nk])
        self.state_nk=g.nz if full_column else self.nk
        self.file.attrs['state_nz_saved']=self.state_nk
        self.file.attrs['full_column_state']=bool(full_column)
        self.file.attrs['budget_top_m']=top
        self.file.attrs['budget_nz_with_halo']=self.nk
        gg=self.file.create_group('grid')
        for name in ('xc','yc','zc','xf','yf','zf'):
            gg.create_dataset(name,data=self.to(getattr(g,name)))
        gg.attrs['origin_m']=(0.,0.,0.); gg.attrs['periodic_xy']=g.periodic
        bg=self.file.create_group('base')
        for name,value in vars(sim.base).items():
            if isinstance(value,np.ndarray): bg.create_dataset(name,data=value)
        bg.create_dataset('rho0_wface',data=self.to(sim.rho0_wface))
        bg.attrs['f_s-1']=sim.f
        self.snapshots=self.file.create_group('snapshots')
        self.steps=self.file.create_dataset('steps',shape=(0,3),maxshape=(None,3),dtype='f8')
        self.steps.attrs['columns']='t_start,dt,step_index'
        dz=np.asarray(self.to(g.dz_c if getattr(g,'stretched',False) else np.full(g.nz,g.dz)))
        weights=np.asarray(self.to(sim.rho0_c))*float(g.dx)*float(g.dy)*dz
        self.mass_weights=g.xp.asarray(weights[None,None,:])
        self.micro_mass={}; self.micro_latent={}
        sim.coupler.scheme.process_observer=self.microphysics_process
        self.previous_surface={name:np.array(self.to(value),copy=True) for name,value in sim.state.surface_precip.items()}
        self.previous=self.native(sim); self.start_t=sim.t; self.increments={}; self.integrals={}
        self.save(sim,initial=True)

    def microphysics_process(self,name,dq,latent_per_dq):
        xp=dq.__class__.__module__.split('.')[0]
        weighted=(dq*self.mass_weights).sum()
        self.micro_mass[name]=self.micro_mass.get(name,0.0)+weighted
        self.micro_latent[name]=self.micro_latent.get(name,0.0)+weighted*float(latent_per_dq)*1005.0

    def native(self,sim):
        return tuple(np.array(self.to(getattr(sim.state,n)[:,:,:self.nk+(n=='w')]),copy=True) for n in ('u','v','w'))

    def mark(self,sim,name,dt):
        current=self.native(sim)
        if name=='begin':
            if any(not np.array_equal(a,b) for a,b in zip(current,self.previous)):
                raise RuntimeError('Unobserved velocity mutation / unsupported hierarchy coupling')
            n=len(self.steps); self.steps.resize((n+1,3)); self.steps[n]=(sim.t,dt,sim.step)
            return
        change=tuple(a-b for a,b in zip(current,self.previous))
        if name not in self.increments: self.increments[name]=[np.zeros_like(a) for a in change]
        for accum,delta in zip(self.increments[name],change): accum+=delta
        self.previous=current
        if name=='les':
            terms=kinematics(current,self.coords)
            for key in ('vector_stretching','vector_dilatation','vector_advection','stretching','tilting'):
                if key not in self.integrals: self.integrals[key]=np.zeros_like(terms[key])
                self.integrals[key]+=dt*terms[key]
        if name=='transport_microphysics_bcs' and sim.state.t>=self.next_output:
            self.save(sim, completed_step=sim.step+1)
            self.next_output+=self.interval

    def dataset(self,group,name,array,units):
        array=np.asarray(array)
        if array.ndim:
            d=group.create_dataset(name,data=array,compression='gzip',compression_opts=1,shuffle=True)
        else:
            d=group.create_dataset(name,data=array)
        d.attrs['units']=units
        return d

    def save(self,sim,initial=False,completed_step=None):
        st=sim.state; time=float(st.t)
        group=self.snapshots.create_group(f'{self.count:05d}')
        group.attrs['time_s']=time; group.attrs['interval_start_s']=self.start_t
        # Stage callbacks precede the driver's increment; close() follows it.
        group.attrs['step']=sim.step if completed_step is None else completed_step
        # Full column native state retained for 3D trajectories and saturation.
        for name in ('u','v','w','p','p_dyn','theta','qv','ql','qi','qr','qs','qg','qh','T','rho','RH_w'):
            value=getattr(st,name,None)
            if value is not None:
                value=value[:,:,:self.state_nk+(name=='w')]
                units='m s-1' if name in ('u','v','w') else ('Pa' if name in ('p','p_dyn') else ('K' if name in ('theta','T') else ('kg m-3' if name=='rho' else ('1' if name=='RH_w' else 'kg kg-1'))))
                self.dataset(group,name,self.to(value),units)
        if not initial:
            inc=group.create_group('increments')
            for name,arrays in self.increments.items():
                stage=inc.create_group(name)
                for component,array in zip(('u','v','w'),arrays): self.dataset(stage,component,array,'m s-1')
            ki=group.create_group('kinematic_integrals')
            for name,array in self.integrals.items(): self.dataset(ki,name,array,'s-1')
            micro=group.create_group('microphysics_integrals')
            for name,value in self.micro_mass.items():
                self.dataset(micro,name+'_mass',self.to(value),'kg')
                self.dataset(micro,name+'_latent_energy',self.to(self.micro_latent[name]),'J')
            precip=group.create_group('surface_precipitation_increment')
            area=float(sim.grid.dx)*float(sim.grid.dy)
            for name,value in st.surface_precip.items():
                current=np.asarray(self.to(value)); delta=(current-self.previous_surface[name])*area
                self.dataset(precip,name,delta,'kg')
                self.previous_surface[name]=current.copy()
        group.attrs['complete']=True; self.file.flush()
        self.count+=1; self.start_t=time; self.increments={}; self.integrals={}; self.micro_mass={}; self.micro_latent={}

    def close(self,sim,completed=True):
        if float(sim.state.t)>self.start_t: self.save(sim)
        self.file.attrs['status']='complete' if completed else 'interrupted'
        self.file.flush(); self.file.close()
        if getattr(sim.coupler.scheme,'process_observer',None)==self.microphysics_process:
            sim.coupler.scheme.process_observer=None
