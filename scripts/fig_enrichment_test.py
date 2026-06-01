import csv, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

rows=[]
with open('../data/table5_enrichment_test_S8_dualmodel.csv') as f:
    for r in csv.DictReader(f):
        rows.append(r)

# use S8 (largest sample) as the main panel
shell='S8'
sub=[r for r in rows if r['shell_K']==shell and int(r['nodes'])>=1000]
omega=[int(r['omega_big']) for r in sub]
old=[float(r['obs/pred_old']) for r in sub]
new=[float(r['obs/pred_new']) for r in sub]
meanE_new=[float(r['meanE_new']) for r in sub]
nodes=[int(r['nodes']) for r in sub]

fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,4.6))

# left panel: obs/pred vs omega, comparing the two models
ax1.axhline(1.0,color='gray',lw=1,ls='--',alpha=.7)
ax1.plot(omega,old,'o-',color='#c0392b',lw=2,ms=8,label=r'naive  $\prod p/(p-2)$')
ax1.plot(omega,new,'s-',color='#185FA5',lw=2,ms=8,label=r'corrected  $\prod (p-1)/(p-3)$')
ax1.set_xlabel(r'$\omega_{>3}(N)$  (distinct prime factors $>3$)',fontsize=11)
ax1.set_ylabel('observed / predicted twin rate',fontsize=11)
ax1.set_title(f'Factor-resolved test of the local model  (shell {shell}, $\\sim$10$^8$ nodes)',fontsize=11)
ax1.legend(fontsize=10,loc='upper left')
ax1.grid(alpha=.25); ax1.set_ylim(0.85,1.22)
ax1.set_xticks(omega)

# right panel: observed vs predicted (scatter, corrected model), by shell
ax2.axhline(1.0,color='gray',lw=1,ls='--',alpha=.7)
colors={'S5':'#7fb069','S6':'#e8a33d','S7':'#9b59b6','S8':'#185FA5'}
for sh,c in colors.items():
    ss=[r for r in rows if r['shell_K']==sh and int(r['nodes'])>=1000]
    if not ss: continue
    me=[float(r['meanE_new']) for r in ss]
    op=[float(r['obs/pred_new']) for r in ss]
    ax2.scatter(me,op,color=c,s=55,label=sh,zorder=3,edgecolor='white',lw=.6)
ax2.set_xlabel(r'mean enrichment  $\overline{E_{\rm new}(N)}=\prod(p-1)/(p-3)$',fontsize=11)
ax2.set_ylabel('observed / predicted (corrected model)',fontsize=11)
ax2.set_title('Corrected model holds across shells & enrichment range',fontsize=11)
ax2.legend(fontsize=9,title='shell',ncol=2)
ax2.grid(alpha=.25); ax2.set_ylim(0.85,1.22)

plt.tight_layout()
plt.savefig('fig_enrichment_test.pdf',bbox_inches='tight')
plt.savefig('fig_enrichment_test.png',dpi=150,bbox_inches='tight')
print("figure generated")
