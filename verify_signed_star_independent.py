"""Independent exact enumeration of the five signed d22 star cases.

The complete certificate is verify_signed_star.py. This separately written
checker accumulates by absolute cavity field and checks only the new local part.
"""

from fractions import Fraction as F
from itertools import product
from math import comb
import json
import time
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

M=44
T=F(11,125)
ETA=(1-T)/(1+T)
C=F(119,100)
TARGETS=[F(n,100000) for n in (5860,5253,4650,4080,3610)]

def verify(h,r):
    lam=1/TARGETS[h]; s=M-h-r
    f={i:4*ETA**i/(1+ETA**i)**2 for i in range(0,M+1,2)}
    cosh={i:(ETA**(i//2)+ETA**(-i//2))/2 for i in range(0,M+1,2)}
    negkern={i:2*cosh[i]/C**3-3/C**2 for i in f}
    center={}
    for e in product((-1,1),repeat=r-1):
        act=1+sum(e)
        active=ETA**h*(ETA**act+ETA**(-act))/4
        a=ETA**h*(1-ETA**2)*ETA**(-act)/4
        center[e]=[]
        for n in range(s+1):
            S=2*n-s
            H=active+(ETA**S+ETA**(-S))/4-lam*a
            center[e].append((max(H,F(0)),max(-H,F(0))))
    maxima=[]
    for x in product((-1,1),repeat=r-1):
        for j in range(s+1):
            pos=[F(0)]*(M//2+1);neg=[F(0)]*(M//2+1)
            for e in product((-1,1),repeat=r-1):
                shift=1+sum(a*b for a,b in zip(x,e))
                for p in range(j+1):
                    for q in range(s-j+1):
                        hp,hm=center[e][p+q]
                        weight=comb(j,p)*comb(s-j,q)
                        orient=2*p-2*q+s-2*j
                        for k in range(h+1):
                            L=abs(shift+orient+2*k-h)
                            assert L%2==0
                            w=weight*comb(h,k)
                            pos[L//2]+=hp*w
                            neg[L//2]+=hm*w
            value=sum(pos[i//2]*f[i]+neg[i//2]*negkern[i] for i in f)/2**(M-1)
            maxima.append((value,x,j))
    value,x,j=max(maxima)
    assert value<0,(h,r,float(value),x,j)
    return {'h':h,'r':r,'target':str(TARGETS[h]),'lambda':str(lam),'c':str(C),'cases':len(maxima),'max_M':float(value),'max_x':x,'max_j':j,'exact_sign':'negative'}

if __name__=='__main__':
    for h,r in [(0,1),(0,2),(1,1),(1,2),(2,1)]:
        then=time.monotonic(); result=verify(h,r)
        result['seconds']=time.monotonic()-then
        print(json.dumps(result),flush=True)
