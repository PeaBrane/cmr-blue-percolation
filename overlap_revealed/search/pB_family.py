# Float: compare tangent-point families for the signed-residual p_B certificate (k'=m-1).
import numpy as np
from math import comb
from scipy.optimize import minimize
from explore_pB_opt2 import build, solve
def Gfam(beta,lam,Us,PU,NU,SU,kappa,lamexp,m):
    k=np.exp(-4*beta*Us)-lam; neg=k<=0; pos=~neg
    base=(SU[:,pos]*k[pos]).sum(1)
    c=kappa*(2*np.cosh(beta*(1+Us[neg])))**lamexp*(2*np.cosh(beta)**m)**(1-lamexp); y=1/c
    return base+((3*PU[:,neg]*y**2-4*NU[:,neg]*y**3)*k[neg]).sum(1)
