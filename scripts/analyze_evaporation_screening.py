"""Compare the fixed, predeclared symmetric screening amplitudes."""
from pathlib import Path
import csv,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'outputs/causal_evap_screening_summary';OUT.mkdir(exist_ok=True)
names={.8:'m20',.9:'m10',.95:'m05',1.:'control',1.05:'p05',1.1:'p10',1.2:'p20'}
data={f:[{k:float(v) for k,v in r.items()} for r in csv.DictReader((ROOT/f'outputs/causal_evap_screen_{n}/screening_metrics.csv').open())] for f,n in names.items()}
control=data[1.]; metrics=('cumulative_rain_evaporation_mass_kg','cold_area_m1_km2','theta_v_min_K','theta_v_cold_mean_K','B_min_ms2','B_cold_mean_ms2','B_gradient_max_s2','w_max_ms','total_condensate_kg')
rows=[]
for delta in (.05,.1,.2):
    weak,strong=data[1-delta],data[1+delta]
    row={'delta':delta}
    for key in metrics:
        w,c,s=weak[-1][key],control[-1][key],strong[-1][key]
        row.update({key+'_weak':w,key+'_control':c,key+'_strong':s,key+'_pair_difference':s-w})
    row['area_pair_cells']=row['cold_area_m1_km2_pair_difference']/.36
    row['wmax_largest_relative_difference']=max(abs(row['w_max_ms_weak']/row['w_max_ms_control']-1),abs(row['w_max_ms_strong']/row['w_max_ms_control']-1))
    row['condensate_largest_relative_difference']=max(abs(row['total_condensate_kg_weak']/row['total_condensate_kg_control']-1),abs(row['total_condensate_kg_strong']/row['total_condensate_kg_control']-1))
    rows.append(row)
with (OUT/'screening_comparison.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
selection={'selected_delta':.05,'cases':{'WEAK-EVAP':.95,'CONTROL':1.,'STRONG-EVAP':1.05},'selection_rule':'smallest tested symmetric amplitude with pair separation > one 0.36 km2 mask cell and continuous cold metrics separated, while wmax and condensate remain within 1% of control','comparison_at_2790_s':rows[0],'larger_candidates_retained_for_screening_only':[.1,.2],'new_amplitude_after_results_forbidden':True}
(OUT/'selection.json').write_text(json.dumps(selection,indent=2))
fig,axs=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
for ax,key,title in zip(axs.flat,metrics[:6],('Evaporação acumulada','Área θv′<−1 K','θv′ mínimo','θv′ médio frio','B mínimo','B médio frio')):
    for factor in sorted(data):ax.plot([r['time_s'] for r in data[factor]],[r[key] for r in data[factor]],label=f'{factor:.2f}')
    ax.set(title=title,xlabel='Tempo (s)',ylabel=key);ax.grid(alpha=.2)
axs[0,0].legend(ncol=2,fontsize=8,title='f')
fig.suptitle('Triagem simétrica a partir de 2370 s; amplitude escolhida antes do piloto: δ=0,05')
fig.savefig(OUT/'screening.png',dpi=150);plt.close(fig)
print(json.dumps(selection,indent=2))
