"""General CPM/APM--PP seed extraction and deterministic extension.

Default rules of manuscript Algorithms 1/2. No distance selection.
The reported high-distance families use separately frozen specializations.
Requires only Python standard library and the bundled, self-certifying zsnf.
"""
from fractions import Fraction
from math import gcd,prod
from functools import lru_cache
import zsnf


def inverse(A):
    n=len(A); a=[[Fraction(x) for x in row]+[Fraction(i==j) for j in range(n)] for i,row in enumerate(A)]
    for j in range(n):
        k=next(k for k in range(j,n) if a[k][j]);a[k],a[j]=a[j],a[k]
        c=a[j][j];a[j]=[x/c for x in a[j]]
        for k in range(n):
            if k!=j:
                c=a[k][j];a[k]=[x-c*y for x,y in zip(a[k],a[j])]
    out=[row[n:] for row in a]
    assert all(x.denominator==1 for row in out for x in row)
    return [[int(x) for x in row] for row in out]


def smith(A):
    d,cols=zsnf.smith(A);n=len(A[0]);d += [0]*(n-len(d))
    return d,[[cols[j][i] for j in range(n)] for i in range(n)]


def rows(J,L,M,a=None,P=None):
    n=2*J*L; out=[]
    for i in range(J):
        for j in range(J):
            for u,v in M[i][j]:
                q=(i*L+u,i*L+v,J*L+j*L+u,J*L+j*L+v)
                rho=1 if a is None else a[q[2]]*pow(a[q[0]],-1,P)%P
                if a is not None:
                    assert rho==a[q[3]]*pow(a[q[1]],-1,P)%P
                row=[0]*n
                for k,c in zip(q,(-rho,rho,1,-1)):row[k]+=c
                out.append(row)
    return out


@lru_cache(None)
def roots(e,P):
    return tuple(r for r in range(1,P) if gcd(r,P)==1 and pow(r,e,P)==1)


def multiplicative(A,r,P):
    return [prod(pow(x,c,P) for x,c in zip(r,row))%P for row in A]


def additive(A,x,P):
    return [sum(c*y for c,y in zip(row,x))%P for row in A]


def validate(J,L,P,a,b,M=None):
    if not (P>=2 and J>=1 and L%2==0 and len(a)==len(b)==2*J*L):raise ValueError('Invalid dimensions')
    if not all(gcd(x,P)==1 for x in a):raise ValueError('Nonunit slope')
    recovered=[]
    for i in range(J):
        rr=[]
        for j in range(J):
            groups={}
            for l in range(L):
                x=i*L+l;z=J*L+j*L+l;rho=a[z]*pow(a[x],-1,P)%P
                groups.setdefault((rho,(b[z]-rho*b[x])%P),[]).append(l)
            if any(len(v)%2 for v in groups.values()):raise ValueError('Not a PP input')
            pairs=sorted((v[k],v[k+1]) for v in groups.values() for k in range(0,len(v),2));rr.append(pairs)
        recovered.append(rr)
    if M is None:M=recovered
    if len(M)!=J or any(len(row)!=J for row in M):raise ValueError('Pairing dimensions')
    for i in range(J):
        for j in range(J):
            if any(len(pair)!=2 for pair in M[i][j]) or sorted(q for pair in M[i][j] for q in pair)!=list(range(L)):raise ValueError('Not a partition')
    try:B=rows(J,L,M,a,P)
    except AssertionError:raise ValueError('Invalid slope pairing')
    if any(additive(B,b,P)):raise ValueError('Invalid shift pairing')
    return M


def extract(J,L,P,a,b,M=None):
    a=[x%P for x in a];b=[x%P for x in b];M=validate(J,L,P,a,b,M)
    e,Q=smith(rows(J,L,M));r=multiplicative(inverse(Q),a,P)
    lam=[roots(d,P).index(x) for d,x in zip(e,r)]
    d,V=smith(rows(J,L,M,a,P));c=[P//gcd(P,x) for x in d]
    y=additive(inverse(V),b,P);assert all(x%f==0 for x,f in zip(y,c))
    theta=[x//f for x,f in zip(y,c)]
    data=dict(J=J,L=L,P0=P,M=M,e=e,Q=Q,lambda_indices=lam,theta=theta)
    assert construct(data,P)==(a,b)
    return data


def construct(data,P):
    if P<2:raise ValueError('P must be at least 2')
    J,L=data['J'],data['L'];r=[roots(e,P)[k%len(roots(e,P))] for e,k in zip(data['e'],data['lambda_indices'])]
    a=multiplicative(data['Q'],r,P)
    d,V=smith(rows(J,L,data['M'],a,P));c=[P//gcd(P,x) for x in d]
    b=additive(V,[x*f for x,f in zip(data['theta'],c)],P)
    validate(J,L,P,a,b,data['M'])
    return a,b


def physical(J,L,P,a,b):
    sides=[]
    for side in (0,1):
        rows_out=[]
        for i in range(J):
            for v in range(P):
                support=[]
                for l in range(L):
                    q=side*J*L+i*L+l
                    u=(pow(a[q],-1,P)*(v-b[q]))%P
                    support.append(l*P+u)
                rows_out.append(sum(1<<q for q in support))
        sides.append(rows_out)
    return sides
