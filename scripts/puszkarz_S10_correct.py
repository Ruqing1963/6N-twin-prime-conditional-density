# S10 (10^9..10^10) 上的 Puszkarz 还原检验——用已验证的完整分解+区间筛逻辑。
# 自己重筛孪生，不依赖外部"中心列表"。输出 (a)无权重比 (b)E加权比 (c)实测孪生比。
import numpy as np, math, time
LO = 10**9 // 6 + 1      # N 下界，使 6N >= 10^9
HI = 10**10 // 6          # N 上界，使 6N < 10^10
SEG = 3_000_000
PB = int(math.isqrt(10**10))+1
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
BP=primes_upto(PB)
print(f"[setup] N in [{LO:,},{HI:,}], 6N in [10^9,10^10), base primes={len(BP):,}")

cnt_sqf=cnt_nsf=0; sumE_sqf=sumE_nsf=0.0; tw_sqf=tw_nsf=0
t0=time.time(); n=LO; seg=0
while n<=HI:
    nh=min(n+SEG,HI+1); sz=nh-n
    rem=np.arange(n,nh,dtype=np.int64)
    E=np.ones(sz); has_sq=np.zeros(sz,bool)
    d2=np.zeros(sz,bool); d3=np.zeros(sz,bool)
    for p in BP:
        if p*p>nh-1: break
        f=((n+p-1)//p)*p
        if f>=nh: continue
        idx=np.arange(f-n,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; c=np.zeros(idx.size,np.int8); m=(sub%p)==0
        while m.any(): sub[m]//=p; c[m]+=1; m=(sub%p)==0
        rem[idx]=sub
        has_sq[idx[c>=2]]=True
        if p==2: d2[idx]=True
        elif p==3: d3[idx]=True
        else: E[idx]*=(p-1)/(p-3)
    lo=rem>1
    if lo.any():
        lp=rem[lo].astype(float); E[lo]*=(lp-1)/(lp-3)
    sixN_sqf=(~has_sq)&(~d2)&(~d3)
    # 实测孪生：区间筛 6N±1
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st=max(p*p,((vlo+p-1)//p)*p)
        if st>vhi: continue
        comp[st-vlo:span:p]=True
    N=np.arange(n,nh,dtype=np.int64)
    tw=(~comp[(6*N-1)-vlo])&(~comp[(6*N+1)-vlo])
    cnt_sqf+=int(sixN_sqf.sum()); cnt_nsf+=int((~sixN_sqf).sum())
    sumE_sqf+=float(E[sixN_sqf].sum()); sumE_nsf+=float(E[~sixN_sqf].sum())
    tw_sqf+=int((tw&sixN_sqf).sum()); tw_nsf+=int((tw&~sixN_sqf).sum())
    seg+=1
    if seg%20==0 or nh>HI:
        print(f"  N={nh-1:,} ({100*(nh-LO)/(HI-LO):.1f}%) twins={tw_sqf+tw_nsf:,} {time.time()-t0:.0f}s")
    n=nh
R0=math.pi**2/3-1
print("\n===== S10 (10^9..10^10) Puszkarz recovery =====")
print(f"twin pairs total: {tw_sqf+tw_nsf:,}  (sqfree centre {tw_sqf:,} / nonsqfree {tw_nsf:,})")
print(f"R0 unbiased        = {R0:.4f}")
print(f"(a) unweighted     = {cnt_nsf/cnt_sqf:.4f}")
print(f"(b) E-weighted     = {sumE_nsf/sumE_sqf:.4f}")
print(f"(c) observed twins = {tw_nsf/tw_sqf:.4f}")
tot=tw_nsf/tw_sqf-cnt_nsf/cnt_sqf; exp=sumE_nsf/sumE_sqf-cnt_nsf/cnt_sqf
print(f"explained: {exp:+.4f} ({100*exp/tot:.0f}%)   residual: {tot-exp:+.4f} ({100*(tot-exp)/tot:.0f}%)")
print(f"elapsed {time.time()-t0:.0f}s")
