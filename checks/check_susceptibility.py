from fractions import Fraction as Q
from itertools import product
from collections import deque
import json
import sys

if sys.flags.optimize:
    raise SystemExit("Verification requires assertions: run Python without -O.")

checks = 0
for a,b,c in product([Q(i,4) for i in range(-4,5)], repeat=3):
    for t in [Q(0),Q(1,10),Q(1,2),Q(9,10)]:
        exact = sum(((a+u*b)/(1+u*c))**2 for u in [t,-t])/2
        stated = ((a*a+t*t*b*b)*(1+t*t*c*c)-4*t*t*a*b*c)/(1-t*t*c*c)**2
        aa = (1+3*t*t)/(1-t*t)**2
        bb = t*t*(3+t*t)/(1-t*t)**2
        assert exact == stated
        assert exact <= aa*a*a+bb*b*b
        checks += 1

summaries=[]
for name,n,edges in [
    ('path6',6,[(i,i+1) for i in range(5)]),
    ('cycle6',6,[(i,(i+1)%6) for i in range(6)]),
    ('complete4',4,[(i,j) for i in range(4) for j in range(i+1,4)]),
    ('triangle_with_tail',6,[(0,1),(1,2),(2,0),(2,3),(3,4),(4,5)]),
]:
    degree=max(sum(x in e for e in edges) for x in range(n))
    t=Q(1,10)
    aa=(1+3*t*t)/(1-t*t)**2
    bb=t*t*(3+t*t)/(1-t*t)**2
    kappa=degree*bb*aa**(degree-1)
    assert kappa<1
    masks=list(range(1<<n))
    corr_sq=[Q(0) for mask in masks]
    spins=list(product([-1,1], repeat=n))
    for signs in product([-1,1], repeat=len(edges)):
        weights=[]
        for spin in spins:
            weight=1
            for (i,j),sign in zip(edges,signs):
                weight*=10+sign*spin[i]*spin[j]
            weights.append(weight)
        z=sum(weights)
        for mask in masks:
            moment=0
            for spin,w in zip(spins,weights):
                monomial=1
                for i in range(n):
                    if mask&(1<<i): monomial*=spin[i]
                moment+=w*monomial
            corr_sq[mask]+=Q(moment,z)**2/(1<<len(edges))
    dist=[]
    for i in range(n):
        dd=[n+1]*n
        dd[i]=0
        queue=deque([i])
        while queue:
            x=queue.popleft()
            for u,v in edges:
                y=v if u==x else u if v==x else None
                if y is not None and dd[y]>dd[x]+1:
                    dd[y]=dd[x]+1
                    queue.append(y)
        dist.append(dd)
    max_sus=Q(0)
    for i in range(n):
        sus=sum(corr_sq[(1<<i)^(1<<j)] for j in range(n))
        assert sus <= 1/(1-kappa)
        max_sus=max(max_sus,sus)
        for j in range(n):
            assert corr_sq[(1<<i)^(1<<j)] <= kappa**dist[i][j]
    odds=[mask for mask in masks if mask.bit_count()%2]
    pair_checks=0
    for a,b in product(odds,repeat=2):
        distance=min(dist[i][j] for i in range(n) if a&(1<<i) for j in range(n) if b&(1<<j))
        assert corr_sq[a^b] <= kappa**distance
        pair_checks+=1
    summaries.append(dict(graph=name,sign_assignments=1<<len(edges),odd_pair_checks=pair_checks,max_susceptibility=float(max_sus),bound=float(1/(1-kappa)),kappa=float(kappa)))
print(json.dumps(dict(edge_identity_and_inequality_checks=checks,graphs=summaries),indent=2))
