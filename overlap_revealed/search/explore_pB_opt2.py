# Optimize tangent points c(U)>0 in the signed residual (epigraph + SLSQP), floating point exploration.
import numpy as np
from math import comb
from scipy.optimize import minimize
def sdist(n,p):
    P=np.array([comb(n,i)*p**i*(1-p)**(n-i) for i in range(n+1)])
    return 2*np.arange(n+1)-n, P
def build(beta,m):
    a=np.exp(2*beta)/(2*np.cosh(2*beta)); n1=m-1
    Us=np.arange(-n1,n1+1,2)
    PUm=[];NUm=[];SUm=[]
    for j in range(m):
        vA,PA=sdist(j,a); vB,PB=sdist(n1-j,a)
        A=vA[:,None]; Bp=vB[None,:]; W=PA[:,None]*PB[None,:]
        U=A+Bp; dot=A-Bp
        PUm.append([np.sum(W[U==u]) for u in Us]); NUm.append([np.sum((W*np.cosh(beta*(1+dot)))[U==u]) for u in Us]); SUm.append([np.sum((W*0.25/np.cosh(beta*(1+dot))**2)[U==u]) for u in Us])
    return a,Us,np.array(PUm),np.array(NUm),np.array(SUm)
def solve(beta,lam,Us,PU,NU,SU):
    k=np.exp(-4*beta*Us)-lam
    neg=k<=0; pos=~neg
    base=(SU[:,pos]*k[pos]).sum(1)
    kn=k[neg]; P=PU[:,neg]; N=NU[:,neg]
    y0=P[-1]/(2*N[-1])  # aligned tangent
    def Gv(z):
        y=np.exp(z); return base+((3*P*y**2-4*N*y**3)*kn).sum(1)
    x0=np.concatenate([np.log(y0),[np.max(Gv(np.log(y0)))]])
    cons={'type':'ineq','fun':lambda x: x[-1]-Gv(x[:-1])}
    res=minimize(lambda x:x[-1],x0,constraints=[cons],method='SLSQP',options={'maxiter':2000,'ftol':1e-16})
    return np.max(Gv(res.x[:-1])),np.exp(res.x[:-1]),neg
