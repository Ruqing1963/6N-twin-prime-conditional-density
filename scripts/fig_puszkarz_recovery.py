# Figure 2. Values R0/(a)/(b)/(c) are the 10^10-scale outputs of puszkarz_S10_correct.py.
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
R0=2.2899; Ra=2.2899; Rb=2.4269; Rc=2.4266   # 10^10 数字
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,4.5))
labels=['unbiased\n$R_0=\\pi^2/3-1$','model prediction\n$\\overline{E}$-weighted','observed twins\n(= Puszkarz $R_2$)']
vals=[Ra,Rb,Rc]; colors=['#999999','#185FA5','#c0392b']; y=np.arange(3)
ax1.barh(y,vals,color=colors,height=0.55,zorder=3)
ax1.set_yticks(y); ax1.set_yticklabels(labels,fontsize=10)
ax1.set_xlim(2.27,2.44); ax1.invert_yaxis()
ax1.axvline(R0,color='#999999',ls='--',lw=1,zorder=1)
for yi,v in zip(y,vals): ax1.text(v+0.002,yi,f'{v:.4f}',va='center',fontsize=10,fontweight='bold')
ax1.set_xlabel('nonsquarefree / squarefree ratio at twin centres',fontsize=10.5)
ax1.set_title('Model reproduces Puszkarz bias at $10^{10}$ scale',fontsize=11)
ax1.grid(axis='x',alpha=.25,zorder=0)
total=Rc-Ra; expl=Rb-Ra; resid=Rc-Rb
ax2.bar(['total\n$(c)-(a)$'],[total],color='#c0392b',width=0.5,zorder=3)
ax2.bar(['explained\n$(b)-(a)$'],[expl],color='#185FA5',width=0.5,zorder=3)
ax2.bar(['residual\n$(c)-(b)$'],[resid],color='#e8a33d',width=0.5,zorder=3)
for i,v in enumerate([total,expl,resid]):
    ax2.text(i,v+0.004 if v>=0 else v-0.004,f'{v:+.4f}\n({100*v/total:.0f}%)',
             ha='center',va='bottom' if v>=0 else 'top',fontsize=9.5,fontweight='bold')
ax2.axhline(0,color='black',lw=.8); ax2.set_ylim(-0.02,0.16)
ax2.set_ylabel('contribution to the bias ratio',fontsize=10.5)
ax2.set_title('Decomposition: model $\\sim$100%, residual $-0.2\\%$',fontsize=11)
ax2.grid(axis='y',alpha=.25,zorder=0)
plt.tight_layout()
plt.savefig('fig_puszkarz_recovery.pdf',bbox_inches='tight')
plt.savefig('fig_puszkarz_recovery.png',dpi=150,bbox_inches='tight')
print("updated fig done")
