"""Reconstruct every published matrix and check ranks, CSS, weights and witnesses.

Python >=3.10, standard library only. This does not rerun distance exclusion.
"""
from pathlib import Path
import hashlib
import json
import math
import sys
import zsnf

ROOT = Path(__file__).resolve().parents[1]
if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

def read(path):
    return json.loads((ROOT/path).read_text())

def shifts(P, a, M, theta, J, L):
    N = 2*J*L
    B = []
    for i in range(J):
        for j in range(J):
            for u,v in M[i][j]:
                rho = a[J*L+j*L+u]*pow(a[i*L+u],-1,P) % P
                assert rho == a[J*L+j*L+v]*pow(a[i*L+v],-1,P) % P
                row = [0]*N
                for k,c in ((J*L+j*L+u,1),(J*L+j*L+v,-1),(i*L+u,-rho),(i*L+v,rho)):
                    row[k] += c
                B.append(row)
    diag,Q = zsnf.smith(B)
    diag += [0]*(N-len(diag))
    c = [P//math.gcd(x,P) if x else 1 for x in diag]
    b = [sum(theta[j]*c[j]*Q[j][i] for j in range(N)) % P for i in range(N)]
    assert all(sum(x*y for x,y in zip(r,b)) % P == 0 for r in B)
    return b

def construct(rule, P):
    J,L = rule['J'],rule['L']
    if rule['kind']=='square_zero':
        ratio=P//192
        assert P%192==0 and ratio>0 and ratio&(ratio-1)==0
        t=ratio.bit_length()-1;D=96*2**(t//2)
        a=[(1+D*s)%P for s in rule['S']]
        b=shifts(P,a,rule['M'],rule['theta'],J,L)
    elif rule['kind']=='power9':
        a=[pow(9,s,P) for s in rule['S']]
        b=shifts(P,a,rule['M'],rule['theta'],J,L)
    else:
        assert rule['kind']=='cpm'
        a=[1]*(2*J*L)
        b=[sum(c*v[i] for c,v in zip(rule['theta'],rule['integer_basis']))%P for i in range(2*J*L)]
    matrices=[]
    for start in (0,J):
        H=[]
        for f in range(J):
            for r in range(P):
                row=[]
                for l in range(L):
                    k=(start+f)*L+l
                    u=pow(a[k],-1,P)*(r-b[k])%P
                    row.append(L*u+l if rule['index']=='q=L*u+l' else l*P+u)
                H.append(sorted(row))
        matrices.append(H)
    for i in range(J):
        for j in range(J):
            for u,v in rule['M'][i][j]:
                ru=a[J*L+j*L+u]*pow(a[i*L+u],-1,P)%P
                rv=a[J*L+j*L+v]*pow(a[i*L+v],-1,P)%P
                assert ru==rv
                assert (b[J*L+j*L+u]-b[J*L+j*L+v]-ru*(b[i*L+u]-b[i*L+v]))%P==0
    return dict(a=a,b=b,hx=matrices[0],hz=matrices[1])

def basis(rows):
    piv={}
    for v in rows:
        while v:
            j=v.bit_length()-1
            if j not in piv:
                piv[j]=v;break
            v^=piv[j]
    return piv

def remainder(v,piv):
    while v and (j:=v.bit_length()-1) in piv:
        v^=piv[j]
    return v

def verify_code(code):
    rule=read(code['rule']);P=code['P'];n=code['n']
    rebuilt=construct(rule,P)
    assert rebuilt['a']==code['a'] and rebuilt['b']==code['b']
    assert rebuilt['hx']==code['hx'] and rebuilt['hz']==code['hz']
    H=[]
    for side in ('hx','hz'):
        rows=code[side];cols=[0]*n
        for row in rows:
            assert row==sorted(set(row)) and len(row)==code['check_weight']
            for q in row:
                assert 0<=q<n;cols[q]+=1
        assert set(cols)=={code['column_weight_per_sector']}
        H.append([sum(1<<q for q in row) for row in rows])
    assert all((x&z).bit_count()%2==0 for x in H[0] for z in H[1])
    bx,bz=map(basis,H)
    assert [len(bx),len(bz)]==code['ranks']
    assert n-len(bx)-len(bz)==code['k']
    assert not remainder((1<<n)-1,bx) and not remainder((1<<n)-1,bz)
    for side,w in code['witnesses'].items():
        support=w['support'];assert len(set(support))==w['weight']
        assert all(0<=q<n for q in support)
        v=sum(1<<q for q in support);opposite=H[1] if side=='X' else H[0]
        assert all((v&r).bit_count()%2==0 for r in opposite)
        assert remainder(v,bx if side=='X' else bz)
    assert min(w['weight'] for w in code['witnesses'].values())==code['distance']['upper']
    return dict(id=code['id'],n=n,k=code['k'],physical_checks='PASS',distance_exclusion_recomputed=False)

def main():
    for path,expected in read('SHA256SUMS.json').items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==expected,path
    results=[verify_code(read(row['file'])) for row in read('data/index.json')]
    print(json.dumps(dict(status='PASS',codes=results,scope='Reconstruction, ranks, CSS, row/column weights, parity and upper witnesses; archived lower-bound audits are not replayed'),indent=2))

if __name__=='__main__':
    main()
