#!/usr/bin/env python3
"""
Exact certificate for the 0.99999-scaled regular-decagon one-pleat wrapping candidate.

All target coordinates are finite decimals parsed as fractions, so every target-side
calculation is exact rational arithmetic. Source-side triangle metric calculations
are exact SymPy expressions in sqrt(5).

Checks:
  1. quotient mesh is a connected orientable triangulated 2-sphere;
  2. each of 30 affine triangle maps is strictly 1-Lipschitz;
  3. every pair of target triangles intersects exactly in its prescribed shared
     simplex (empty / vertex / edge), hence the PL sphere is embedded;
  4. an exact face-plane witness certifies nonconvexity;
  5. exact enclosed volume equals the rational value stated in the paper and
     exceeds 0.179843 R^3.

R is normalized to 1.
"""

from fractions import Fraction
import itertools, collections
import sympy as sp

def require(condition, message="certificate check failed"):
    """Raise explicitly if a certificate condition fails (active even under python -O)."""
    if not condition:
        raise AssertionError(message)

def F(s):
    return Fraction(s)

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
SCALE = Fraction(99999,100000)
BASE = {u: tuple(F(x) for x in xyz) for u, xyz in BASE_STR.items()}
Q = {u: tuple(SCALE*x for x in xyz) for u, xyz in BASE.items()}

N = 10
rho = sp.Rational(13,20)

def traversal_labels(ma,mb,mc):
    J=('J',0)
    def nd(a,k): return J if k==0 else (a,k)
    s=[nd('a',ma)]
    for k in range(ma-1,-1,-1): s.append(nd('a',k))
    for k in range(1,mc+1): s.append(nd('c',k))
    for k in range(mc-1,-1,-1): s.append(nd('c',k))
    for k in range(1,mb+1): s.append(nd('b',k))
    for k in range(mb-1,-1,-1): s.append(nd('b',k))
    for k in range(1,ma+1): s.append(nd('a',k))
    require(s[0] == s[-1])
    return s[:-1]

boundary = traversal_labels(2,2,1)
C=('center',0)
I=[('inner',i) for i in range(N)]

faces_orig=[]
for i in range(N):
    j=(i+1)%N
    faces_orig.append((('c',0),('i',i),('i',j)))
    faces_orig += [
        (('i',i),('b',i),('b',j)),
        (('i',i),('b',j),('i',j)),
    ]

def vl(v):
    typ,k=v
    if typ=='c': return C
    if typ=='i': return I[k]
    return boundary[k]

faces=[tuple(vl(v) for v in f) for f in faces_orig]
nodes=sorted(set(x for f in faces for x in f), key=repr)

# ---------- exact source triangle metric ----------
def src_radius(v):
    typ,k=v
    if typ=='c': return sp.Rational(0)
    if typ=='i': return rho
    return sp.Rational(1)

def src_index(v):
    return 0 if v[0]=='c' else v[1]

def src_d2(a,b):
    ra,rb=src_radius(a),src_radius(b)
    if ra==0 or rb==0:
        return sp.simplify(ra**2+rb**2)
    dk=(src_index(a)-src_index(b)) % N
    return sp.simplify(ra**2+rb**2
                       -2*ra*rb*sp.cos(sp.pi*sp.Rational(dk,5)))

def sub(a,b): return tuple(x-y for x,y in zip(a,b))
def add(a,b): return tuple(x+y for x,y in zip(a,b))
def mul(s,a): return tuple(s*x for x in a)
def dot(a,b): return sum(x*y for x,y in zip(a,b))
def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],
            a[2]*b[0]-a[0]*b[2],
            a[0]*b[1]-a[1]*b[0])
def zero(v): return all(x==0 for x in v)

def tgt_d2(u,v):
    d=sub(Q[u],Q[v])
    return dot(d,d)

# Strict 1-Lipschitz: for each triangle, source Gram - target Gram is positive definite.
lip_records=[]
for fi,f in enumerate(faces_orig):
    u=[vl(x) for x in f]
    s01,s02,s12=src_d2(f[0],f[1]),src_d2(f[0],f[2]),src_d2(f[1],f[2])
    S00,S11=s01,s02
    S01=sp.simplify((s01+s02-s12)/2)
    t01,t02,t12=(sp.Rational(tgt_d2(u[0],u[1]).numerator,tgt_d2(u[0],u[1]).denominator),
                  sp.Rational(tgt_d2(u[0],u[2]).numerator,tgt_d2(u[0],u[2]).denominator),
                  sp.Rational(tgt_d2(u[1],u[2]).numerator,tgt_d2(u[1],u[2]).denominator))
    T00,T11=t01,t02
    T01=(t01+t02-t12)/2
    m00=sp.simplify(S00-T00)
    m11=sp.simplify(S11-T11)
    det=sp.simplify(m00*m11-(S01-T01)**2)
    def positive_in_Qsqrt5(expr):
        # Every source-metric expression here lies in Q(sqrt(5)).
        # Decide a+b*sqrt(5)>0 using rational arithmetic only.
        e=sp.expand(sp.radsimp(sp.simplify(expr)))
        rt=sp.sqrt(5)
        b=sp.simplify(e.coeff(rt))
        a=sp.simplify(e-b*rt)
        require(a.is_Rational and b.is_Rational)
        a=Fraction(int(sp.numer(a)),int(sp.denom(a)))
        b=Fraction(int(sp.numer(b)),int(sp.denom(b)))
        if b==0:
            return a>0
        if a>=0 and b>0:
            return True
        if a<=0 and b<0:
            return False
        if b>0 and a<0:
            return 5*b*b > a*a
        if b<0 and a>0:
            return a*a > 5*b*b
        return False

    require(positive_in_Qsqrt5(m00))
    require(positive_in_Qsqrt5(m11))
    require(positive_in_Qsqrt5(det))
    lip_records.append((m00,m11,det))

# ---------- combinatorial closed sphere ----------
edge_faces=collections.defaultdict(list)
adj=collections.defaultdict(set)
for fi,f in enumerate(faces):
    require(len(set(f))==3)
    for a,b in ((f[0],f[1]),(f[1],f[2]),(f[2],f[0])):
        e=tuple(sorted((a,b),key=repr))
        edge_faces[e].append(fi)
        adj[a].add(b); adj[b].add(a)
require(len(nodes)==17 and len(edge_faces)==45 and len(faces)==30)
require(all(len(fs)==2 for fs in edge_faces.values()))
require(len(nodes)-len(edge_faces)+len(faces)==2)

# connected
seen={nodes[0]}; stack=[nodes[0]]
while stack:
    u=stack.pop()
    for v in adj[u]:
        if v not in seen:
            seen.add(v); stack.append(v)
require(len(seen)==len(nodes))

# every vertex link is one cycle
for u in nodes:
    L=collections.defaultdict(set)
    for f in faces:
        if u in f:
            a,b=[x for x in f if x!=u]
            L[a].add(b); L[b].add(a)
    require(all(len(L[x])==2 for x in L))
    s=next(iter(L)); seenL={s}; st=[s]
    while st:
        x=st.pop()
        for y in L[x]:
            if y not in seenL:
                seenL.add(y); st.append(y)
    require(len(seenL)==len(L))

# face orientations agree along every shared edge
def oriented_edges(f):
    return ((f[0],f[1]),(f[1],f[2]),(f[2],f[0]))
for e,fs in edge_faces.items():
    u,v=e
    dirs=[]
    for fi in fs:
        oe=oriented_edges(faces[fi])
        dirs.append(1 if (u,v) in oe else -1 if (v,u) in oe else 0)
    require(sorted(dirs)==[-1,1])

# ---------- exact pairwise triangle intersection ----------
def normal(face):
    A,B,C0=[Q[u] for u in face]
    n=cross(sub(B,A),sub(C0,A))
    require(not zero(n))
    return n

def plane_section(face,n,p0):
    pts=[Q[u] for u in face]
    ds=[dot(n,sub(x,p0)) for x in pts]
    if all(d>0 for d in ds) or all(d<0 for d in ds):
        return []
    out=[]
    for i,d in enumerate(ds):
        if d==0:
            out.append(pts[i])
    for i,j in ((0,1),(1,2),(2,0)):
        di,dj=ds[i],ds[j]
        if di*dj<0:
            t=di/(di-dj)
            out.append(add(pts[i],mul(t,sub(pts[j],pts[i]))))
    uniq=[]
    for x in out:
        if x not in uniq: uniq.append(x)
    return uniq

def interval(sec,d):
    if not sec: return None
    p=[dot(x,d) for x in sec]
    return min(p),max(p)

pair_counts=collections.Counter()
for i,j in itertools.combinations(range(len(faces)),2):
    A,B=faces[i],faces[j]
    nA,nB=normal(A),normal(B)
    line=cross(nA,nB)
    # In this explicit certificate no two face planes are parallel.
    require(not zero(line))
    IA=interval(plane_section(A,nB,Q[B[0]]),line)
    IB=interval(plane_section(B,nA,Q[A[0]]),line)
    if IA is None or IB is None:
        inter=None
    else:
        lo=max(IA[0],IB[0]); hi=min(IA[1],IB[1])
        inter=(lo,hi) if lo<=hi else None

    shared=set(A)&set(B)
    pair_counts[len(shared)]+=1
    if len(shared)==0:
        require(inter is None)
    elif len(shared)==1:
        s=next(iter(shared)); t=dot(Q[s],line)
        require(inter==(t,t))
    elif len(shared)==2:
        vals=[dot(Q[s],line) for s in shared]
        require(inter==(min(vals),max(vals)))
    else:
        raise AssertionError("duplicate face")

require(pair_counts==collections.Counter({0:260,1:130,2:45}), "unexpected triangle-pair incidence counts")

# ---------- exact nonconvexity witness ----------
f0=faces[0]
A0,B0,C0=[Q[x] for x in f0]
n0=cross(sub(B0,A0),sub(C0,A0))
pos=dot(n0,sub(Q[('J',0)],A0))
neg=dot(n0,sub(Q[('a',2)],A0))
require(pos>0 and neg<0, "nonconvexity witness failed")

# ---------- exact enclosed volume ----------
def det3(a,b,c): return dot(a,cross(b,c))
V=sum((det3(Q[f[0]],Q[f[1]],Q[f[2]]) for f in faces),Fraction(0))/6
Vabs=abs(V)
EXPECTED=Fraction(
    359686988325832312742347877865122384142690317017149,
    2000000000000000000000000000000000000000000000000000,
)
require(Vabs == EXPECTED, "exact volume does not match the paper")
require(Vabs > Fraction("0.179843"), "certified volume lower bound failed")

print("CERTIFIED")
print("V,E,F =",len(nodes),len(edge_faces),len(faces),"Euler =",len(nodes)-len(edge_faces)+len(faces))
print("triangle pairs =",dict(pair_counts),"(all intersections exactly prescribed)")
print("nonconvex witness: face 0 has J and a2 on opposite sides")
print("volume exact =",Vabs)
print("volume decimal =",float(Vabs))
print("volume > 0.179843 =",Vabs > Fraction("0.179843"))
print("all 30 triangle maps are strictly 1-Lipschitz")
