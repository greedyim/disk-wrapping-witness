#!/usr/bin/env python3
"""
Independent exact verifier for the one-pleat decagon wrapping certificate.

This verifier intentionally uses a different proof strategy from
`exact_certificate.py`:

* source regular-decagon metric is represented by a hand-written quadratic
  field Q(sqrt(5)); no SymPy is used;
* triangle-triangle intersection is decided as an exact barycentric linear
  feasibility problem, not by plane-section interval overlap;
* orientability is established by a face-sign propagation on the dual graph;
* nonconvexity is certified by an exact face-plane witness.

All target coordinates are exact Fractions obtained from finite decimal strings.
R is normalized to 1.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from collections import defaultdict, deque, Counter
from itertools import combinations

def require(condition, message="certificate check failed"):
    """Raise explicitly if a certificate condition fails (active even under python -O)."""
    if not condition:
        raise AssertionError(message)

# ---------- Q(sqrt(5)) exact arithmetic ----------
@dataclass(frozen=True)
class Q5:
    a: F = F(0)  # a + b sqrt(5)
    b: F = F(0)
    def __add__(self, o):
        o=q5(o); return Q5(self.a+o.a,self.b+o.b)
    __radd__=__add__
    def __neg__(self): return Q5(-self.a,-self.b)
    def __sub__(self,o): return self+(-q5(o))
    def __rsub__(self,o): return q5(o)-self
    def __mul__(self,o):
        o=q5(o); return Q5(self.a*o.a+5*self.b*o.b,self.a*o.b+self.b*o.a)
    __rmul__=__mul__
    def __truediv__(self,o):
        o=q5(o); den=o.a*o.a-5*o.b*o.b
        if den==0: raise ZeroDivisionError
        return Q5((self.a*o.a-5*self.b*o.b)/den,(self.b*o.a-self.a*o.b)/den)
    def __eq__(self,o):
        o=q5(o); return self.a==o.a and self.b==o.b
    def sign(self):
        a,b=self.a,self.b
        if b==0: return (a>0)-(a<0)
        if a==0: return (b>0)-(b<0)
        if a>0 and b>0: return 1
        if a<0 and b<0: return -1
        cmp=(a*a > 5*b*b) - (a*a < 5*b*b)
        if cmp==0: return 0
        if a>0 and b<0: return cmp
        if a<0 and b>0: return -cmp
        raise AssertionError
    def __gt__(self,o): return (self-q5(o)).sign()>0
    def __ge__(self,o): return (self-q5(o)).sign()>=0
    def __repr__(self): return f"Q5({self.a!r},{self.b!r})"

def q5(x):
    return x if isinstance(x,Q5) else Q5(F(x),F(0))

SQ5=Q5(F(0),F(1))
COS = {
    0: Q5(F(1),0),
    1: Q5(F(1,4),F(1,4)),
    2: Q5(F(-1,4),F(1,4)),
    3: Q5(F(1,4),F(-1,4)),
    4: Q5(F(-1,4),F(-1,4)),
    5: Q5(F(-1),0),
}
def cos36k(k):
    k%=10
    if k>5: k=10-k
    return COS[k]

# ---------- target certificate coordinates ----------
BASE_STR = {
    ('J', 0): ('0.114996453431', '0.652042463578', '0.359274110539'),
    ('a', 1): ('0.713527595885', '0.556758865942', '0.238247176413'),
    ('a', 2): ('0.904510611077', '0.150919970833', '-0.186942534628'),
    ('b', 1): ('-0.269259146678', '0.427202794161', '0.787947375499'),
    ('b', 2): ('-0.41123294113', '-0.173302559313', '0.822628398965'),
    ('c', 1): ('-0.407447365206', '0.632905197995', '0.029649729336'),
    ('center', 0): ('0', '0', '0'),
    ('inner', 0): ('0.65', '0', '0'),
    ('inner', 1): ('0.525861046344', '0.38206041399', '0'),
    ('inner', 2): ('0.201227167335', '0.617682813279', '0.021807551587'),
    ('inner', 3): ('-0.138361314859', '0.605243587298', '-0.192448295874'),
    ('inner', 4): ('-0.192247060016', '0.586039489637', '0.205179883277'),
    ('inner', 5): ('-0.33215460465', '0.312020954927', '0.463482731389'),
    ('inner', 6): ('-0.42443756412', '-0.078307526511', '0.486025395896'),
    ('inner', 7): ('-0.126467920004', '0.156558885543', '0.618057586773'),
    ('inner', 8): ('0.208351084601', '0.348929212409', '0.507285156764'),
    ('inner', 9): ('0.525861046344', '0.263795278601', '0.276373318043'),
}
SCALE=F(99999,100000)
Q={u:tuple(SCALE*F(s) for s in xyz) for u,xyz in BASE_STR.items()}
N=10; RHO=F(13,20)

boundary=[('a',2),('a',1),('J',0),('c',1),('J',0),('b',1),('b',2),('b',1),('J',0),('a',1)]
C=('center',0); I=[('inner',i) for i in range(N)]

def vl(v):
    typ,k=v
    if typ=='c': return C
    if typ=='i': return I[k]
    return boundary[k]

faces_orig=[]
for i in range(N):
    j=(i+1)%N
    faces_orig.append((('c',0),('i',i),('i',j)))
    faces_orig.append((('i',i),('b',i),('b',j)))
    faces_orig.append((('i',i),('b',j),('i',j)))
faces=[tuple(vl(v) for v in f) for f in faces_orig]
nodes=sorted(set(v for f in faces for v in f),key=repr)

# ---------- vector arithmetic over rationals ----------
def vsub(a,b): return tuple(x-y for x,y in zip(a,b))
def vadd(a,b): return tuple(x+y for x,y in zip(a,b))
def vmul(s,a): return tuple(s*x for x in a)
def dot(a,b): return sum((x*y for x,y in zip(a,b)),F(0))
def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
def det3(a,b,c): return dot(a,cross(b,c))

# ---------- source metric and exact Lipschitz test ----------
def sradius(v):
    return F(0) if v[0]=='c' else (RHO if v[0]=='i' else F(1))
def sindex(v): return 0 if v[0]=='c' else v[1]
def sd2(a,b):
    ra,rb=sradius(a),sradius(b)
    if ra==0 or rb==0: return q5(ra*ra+rb*rb)
    return q5(ra*ra+rb*rb) - q5(2*ra*rb)*cos36k(sindex(a)-sindex(b))
def td2(a,b):
    d=vsub(Q[a],Q[b]); return dot(d,d)

def target_gram(vertices):
    p0,p1,p2=[Q[v] for v in vertices]
    u=vsub(p1,p0); v=vsub(p2,p0)
    return ((dot(u,u),dot(u,v)),(dot(u,v),dot(v,v)))

def source_gram(f):
    d01,d02,d12=sd2(f[0],f[1]),sd2(f[0],f[2]),sd2(f[1],f[2])
    off=(d01+d02-d12)/q5(2)
    return ((d01,off),(off,d02))

lip=[]
for k,(fo,ft) in enumerate(zip(faces_orig,faces)):
    S=source_gram(fo); T=target_gram(ft)
    m00=S[0][0]-q5(T[0][0]); m01=S[0][1]-q5(T[0][1]); m11=S[1][1]-q5(T[1][1])
    det=m00*m11-m01*m01
    require(m00>0 and m11>0 and det>0, f"Lipschitz failure face {k}")
    lip.append((m00,m11,det))

# ---------- combinatorics and orientability ----------
edge_faces=defaultdict(list)
for fi,f in enumerate(faces):
    require(len(set(f))==3)
    for i in range(3):
        e=frozenset((f[i],f[(i+1)%3])); edge_faces[e].append(fi)
require((len(nodes),len(edge_faces),len(faces))==(17,45,30))
require(all(len(x)==2 for x in edge_faces.values()))
require(len(nodes)-len(edge_faces)+len(faces)==2)
G=defaultdict(set)
for e in edge_faces:
    a,b=tuple(e); G[a].add(b);G[b].add(a)
seen={nodes[0]}; dq=deque(seen)
while dq:
    u=dq.popleft()
    for w in G[u]:
        if w not in seen: seen.add(w);dq.append(w)
require(len(seen)==len(nodes))
for u in nodes:
    L=defaultdict(set)
    for f in faces:
        if u in f:
            x,y=[z for z in f if z!=u]; L[x].add(y);L[y].add(x)
    require(L and all(len(L[x])==2 for x in L))
    s=next(iter(L)); ss={s}; dq=deque([s])
    while dq:
        x=dq.popleft()
        for y in L[x]:
            if y not in ss: ss.add(y);dq.append(y)
    require(len(ss)==len(L))
def edge_dir(face,a,b):
    for i in range(3):
        if face[i]==a and face[(i+1)%3]==b: return 1
        if face[i]==b and face[(i+1)%3]==a: return -1
    raise AssertionError
orient={0:1}; dq=deque([0])
while dq:
    fi=dq.popleft()
    for e in [frozenset((faces[fi][i],faces[fi][(i+1)%3])) for i in range(3)]:
        a,b=tuple(e); fs=edge_faces[e]; fj=fs[0] if fs[1]==fi else fs[1]
        need=-orient[fi]*edge_dir(faces[fi],a,b)*edge_dir(faces[fj],a,b)
        if fj in orient: require(orient[fj]==need, "orientation propagation inconsistency")
        else: orient[fj]=need;dq.append(fj)
require(len(orient)==len(faces))
require(set(orient.values())=={1}, "listed face order is not globally consistent")

# ---------- exact triangle intersection via barycentric feasibility ----------
def rref_affine_line(M,b):
    A=[list(row)+[rhs] for row,rhs in zip(M,b)]
    rows=len(A); cols=4; pivot_cols=[]; r=0
    for c in range(cols):
        piv=next((i for i in range(r,rows) if A[i][c]!=0),None)
        if piv is None: continue
        A[r],A[piv]=A[piv],A[r]
        p=A[r][c]; A[r]=[x/p for x in A[r]]
        for i in range(rows):
            if i!=r and A[i][c]!=0:
                q=A[i][c]; A[i]=[A[i][j]-q*A[r][j] for j in range(cols+1)]
        pivot_cols.append(c); r+=1
        if r==rows: break
    require(r==3, "unexpected parallel/coplanar face planes")
    free=next(c for c in range(cols) if c not in pivot_cols)
    x0=[F(0)]*cols; d=[F(0)]*cols; d[free]=F(1)
    for i,c in enumerate(pivot_cols):
        x0[c]=A[i][cols]; d[c]=-A[i][free]
    return x0,d,free

def intersect_ge0_interval(lo,hi,c0,c1):
    if c1==0:
        return (lo,hi) if c0>=0 else None
    bound=-c0/c1
    if c1>0:
        lo=bound if lo is None or bound>lo else lo
    else:
        hi=bound if hi is None or bound<hi else hi
    if lo is not None and hi is not None and lo>hi: return None
    return lo,hi

def triangle_intersection_interval(Aface,Bface):
    A0,A1,A2=[Q[x] for x in Aface]; B0,B1,B2=[Q[x] for x in Bface]
    u=vsub(A1,A0);v=vsub(A2,A0);p=vsub(B1,B0);q=vsub(B2,B0)
    M=[[u[r],v[r],-p[r],-q[r]] for r in range(3)]
    rhs=list(vsub(B0,A0))
    x0,d,free=rref_affine_line(M,rhs)
    forms=[(x0[0],d[0]),(x0[1],d[1]),(1-x0[0]-x0[1],-d[0]-d[1]),
           (x0[2],d[2]),(x0[3],d[3]),(1-x0[2]-x0[3],-d[2]-d[3])]
    lo=hi=None
    for c0,c1 in forms:
        z=intersect_ge0_interval(lo,hi,c0,c1)
        if z is None: return None
        lo,hi=z
    require(lo is not None and hi is not None)
    return lo,hi,x0,d

def point_on_A(sol,t,Aface):
    lo,hi,x0,d=sol
    a=x0[0]+d[0]*t; b=x0[1]+d[1]*t
    A0,A1,A2=[Q[x] for x in Aface]
    return vadd(A0,vadd(vmul(a,vsub(A1,A0)),vmul(b,vsub(A2,A0))))

pair_counts=Counter()
for i,j in combinations(range(30),2):
    A,B=faces[i],faces[j]; shared=set(A)&set(B); pair_counts[len(shared)]+=1
    sol=triangle_intersection_interval(A,B)
    if not shared:
        require(sol is None, f"unexpected intersection {i},{j}")
    elif len(shared)==1:
        require(sol is not None)
        lo,hi,_,_=sol; require(lo==hi, "vertex-only intersection is not a single point")
        require(point_on_A(sol,lo,A)==Q[next(iter(shared))])
    elif len(shared)==2:
        require(sol is not None)
        lo,hi,_,_=sol
        got={point_on_A(sol,lo,A),point_on_A(sol,hi,A)}
        expect={Q[x] for x in shared}
        require(got==expect)
    else: raise AssertionError
require(pair_counts==Counter({0:260,1:130,2:45}))

# ---------- exact nonconvexity witness ----------
f0=faces[0]; A,B,C0=[Q[x] for x in f0]; n=cross(vsub(B,A),vsub(C0,A))
pos=dot(n,vsub(Q[('J',0)],A)); neg=dot(n,vsub(Q[('a',2)],A))
require(pos>0 and neg<0, 'nonconvexity witness failed')

# ---------- exact volume ----------
V=sum((det3(Q[f[0]],Q[f[1]],Q[f[2]]) for f in faces),F(0))/6
V=abs(V)
EXPECTED=F(359686988325832312742347877865122384142690317017149,
           2000000000000000000000000000000000000000000000000000)
require(V==EXPECTED, 'exact volume does not match the paper')
require(V>F(179843,1000000), 'certified volume lower bound failed')

print('INDEPENDENT CERTIFICATE VERIFIED')
print('V,E,F =',len(nodes),len(edge_faces),len(faces),'Euler =',len(nodes)-len(edge_faces)+len(faces))
print('triangle pairs =',dict(pair_counts),'(barycentric exact-feasibility method)')
print('nonconvex witness: face 0 has J and a2 on opposite sides')
print('volume exact =',V)
print('volume decimal =',float(V))
print('all 30 triangle maps strictly 1-Lipschitz in Q(sqrt(5))')
