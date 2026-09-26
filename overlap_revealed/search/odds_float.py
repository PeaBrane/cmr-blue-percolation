# Floating-point version of the odds certificate pipeline (to choose tangent weights); exact verification is separate.
import numpy as np
from math import comb
from scipy.optimize import minimize_scalar
def binom_pmf(n,p): return np.array([comb(n,i)*p**i*(1-p)**(n-i) for i in range(n+1)])
def num_tangent(beta,m,k,wmix):
    """tangent point c_i (i=#plus in N+ block, k tilted coords) from nu0=(1-w)aligned+w one-off (on N+) x fair (N-).
       returns min over j+ of Xtilde(j+), and the per-j values"""
    km=m-k; a=np.exp(2*beta)/(2*np.cosh(2*beta))
    i=np.arange(k+1); S=2*i-k
    if k>0:
        R=(1-wmix)*2*np.cosh(beta*S)+wmix*((i/k)*2*np.cosh(beta*(S-2))+((k-i)/k)*2*np.cosh(beta*(S+2)))
    else:
        R=np.array([2.0])
    c=np.cosh(beta)**km*R
    fair=binom_pmf(km,0.5); Fv=2*np.arange(km+1)-km
    vals=[]
    for jp in range(k//2+1):
        PA=binom_pmf(jp,a); PB=binom_pmf(k-jp,a)
        tot=0.0
        for al in range(jp+1):
            for b in range(k-jp+1):
                A=2*al-jp; B=2*b-(k-jp); z=A-B; ii=al+b
                EW=np.sum(fair*2*np.cosh(beta*(z+Fv)))
                tot+=PA[al]*PB[b]*(3/c[ii]**2-2*EW/c[ii]**3)
        vals.append(tot)
    return min(vals),vals
def den_max(beta,m,k):
    km=m-k; a=np.exp(2*beta)/(2*np.cosh(2*beta))
    best=0;arg=None
    for jm in range(km//2+1):
        P=np.array([1.0])
        for _ in range(jm): P=np.convolve(P,[1-a,a])
        for _ in range(km-jm): P=np.convolve(P,[a,1-a])
        for _ in range(k): P=np.convolve(P,[0.5,0.5])
        Z=2*np.arange(m+1)-m
        val=np.sum(P*0.25/np.cosh(beta*Z)**2)
        if val>best: best=val;arg=jm
    return best,arg
def best_w(beta,m,k):
    if k<2: return 0.0
    f=lambda w: -num_tangent(beta,m,k,w)[0]
    res=minimize_scalar(f,bounds=(0,0.9),method='bounded',options={'xatol':1e-10})
    w0=-f(0.0)
    return res.x if -res.fun>w0 else 0.0
