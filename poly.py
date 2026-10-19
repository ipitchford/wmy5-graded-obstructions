"""Small exact sparse polynomial implementation; Python standard library only.

This verifier does not search for a Groebner basis. It differentiates the stated
maps, expands determinants, and checks supplied polynomial identities.
"""
from fractions import Fraction
from functools import lru_cache

class Ring:
    def __init__(self, names):
        self.names=tuple(names); self.n=len(self.names)
        if len(set(self.names))!=self.n: raise ValueError('Duplicate variable')
        self.index={name:i for i,name in enumerate(self.names)}
        self.zero_monom=(0,)*self.n
    def const(self,c): return Poly(self,{self.zero_monom:c} if c else {})
    def var(self,name):
        m=[0]*self.n;m[self.index[name]]=1;return Poly(self,{tuple(m):1})
    def vars(self): return tuple(self.var(n) for n in self.names)

class Poly:
    def __init__(self,R,terms):
        self.R=R;self.d={tuple(m):c for m,c in terms.items() if c}
        if any(len(m)!=R.n or min(m,default=0)<0 for m in self.d): raise ValueError('Bad monomial')
    def coerce(self,p):
        if isinstance(p,Poly):
            if p.R.names!=self.R.names: raise ValueError('Different polynomial rings')
            return p
        return self.R.const(p)
    def __add__(self,p):
        p=self.coerce(p);d=self.d.copy()
        for m,c in p.d.items():
            v=d.get(m,0)+c
            if v:d[m]=v
            elif m in d:del d[m]
        return Poly(self.R,d)
    __radd__=__add__
    def __neg__(self):return Poly(self.R,{m:-c for m,c in self.d.items()})
    def __sub__(self,p):return self+-self.coerce(p)
    def __rsub__(self,p):return self.coerce(p)+-self
    def __mul__(self,p):
        if not isinstance(p,Poly):return Poly(self.R,{m:c*p for m,c in self.d.items() if c*p})
        p=self.coerce(p);d={}
        for m,c in self.d.items():
            for n,e in p.d.items():
                k=tuple(a+b for a,b in zip(m,n));v=d.get(k,0)+c*e
                if v:d[k]=v
                elif k in d:del d[k]
        return Poly(self.R,d)
    __rmul__=__mul__
    def __pow__(self,k):
        if not isinstance(k,int) or k<0:raise ValueError('Nonnegative integer powers only')
        p=self.R.const(1);b=self
        while k:
            if k&1:p=p*b
            k//=2
            if k:b=b*b
        return p
    def __eq__(self,p):
        try:p=self.coerce(p)
        except ValueError:return False
        return self.d==p.d
    def __bool__(self):return bool(self.d)
    def key(self):return tuple(sorted(self.d.items()))
    def diff(self,name):
        i=self.R.index[name];d={}
        for m,c in self.d.items():
            if m[i]:
                n=list(m);n[i]-=1;d[tuple(n)]=c*m[i]
        return Poly(self.R,d)
    def set(self,name,value):
        i=self.R.index[name];d={}
        for m,c in self.d.items():
            n=list(m);n[i]=0;n=tuple(n);v=d.get(n,0)+c*value**m[i]
            if v:d[n]=v
            elif n in d:del d[n]
        return Poly(self.R,d)
    def value(self,values):
        out=0
        for m,c in self.d.items():
            for name,k in zip(self.R.names,m):
                if k:c*=values[name]**k
            out+=c
        return out
    def mul_term(self,m,c):
        return Poly(self.R,{tuple(a+b for a,b in zip(n,m)):v*c for n,v in self.d.items()})
    def reindex(self,R):
        d={}
        for m,c in self.d.items():
            n=[0]*R.n
            for name,k in zip(self.R.names,m):
                if k:
                    if name not in R.index:raise ValueError('Cannot drop used variable')
                    n[R.index[name]]=k
            n=tuple(n);d[n]=d.get(n,0)+c
        return Poly(R,d)
    def subst(self,R,images):
        out=R.const(0);cache={}
        for m,c in self.d.items():
            p=R.const(c)
            for name,k in zip(self.R.names,m):
                if k:
                    if (name,k) not in cache:
                        image=images[name] if name in images else R.var(name)
                        if not isinstance(image,Poly):image=R.const(image)
                        cache[name,k]=image**k
                    p=p*cache[name,k]
                    if not p:break
            out=out+p
        return out

def determinant(M):
    n=len(M);R=M[0][0].R
    @lru_cache(None)
    def rec(row,cols):
        if not cols:return R.const(1)
        out=R.const(0)
        for k,j in enumerate(cols):
            if M[row][j]:
                out=out+((-1)**k)*M[row][j]*rec(row+1,cols[:k]+cols[k+1:])
        return out
    return rec(0,tuple(range(n)))

def charpoly(H,variables,lam='L',specialise_x=False):
    R=H[0].R;l=R.var(lam)
    M=[[ (l if i==j else 0)-H[i].diff(v) for j,v in enumerate(variables)] for i in range(len(H))]
    if specialise_x:M=[[p.set('x',1) for p in row] for row in M]
    return determinant(M)

def coefficients(p,physical,parameters):
    R=Ring(parameters);out={};ip=[p.R.index[n] for n in physical];ir=[p.R.index[n] for n in parameters]
    covered=set(ip+ir)
    for m,c in p.d.items():
        if any(k and i not in covered for i,k in enumerate(m)):raise ValueError('Unaccounted variable')
        key=tuple(m[i] for i in ip);mon=tuple(m[i] for i in ir)
        if key not in out:out[key]={}
        out[key][mon]=out[key].get(mon,0)+c
    ps=[Poly(R,d) for _,d in sorted(out.items(),reverse=True)]
    return unique([p for p in ps if p])

def unique(ps):
    out=[];seen=set()
    for p in ps:
        if p.key() not in seen:out.append(p);seen.add(p.key())
    return out

def decode(R,terms):
    d={}
    for m,n,den in terms:
        if len(m)!=R.n or any(not isinstance(x,int) or x<0 for x in m):raise ValueError('Invalid exponent')
        if not isinstance(n,int) or not isinstance(den,int) or den==0:raise ValueError('Invalid rational')
        key=tuple(m)
        if key in d:raise ValueError('Duplicate monomial in encoding')
        d[key]=Fraction(n,den)
    return Poly(R,d)

def bareiss(M):
    a=[list(map(int,row)) for row in M];n=len(a);sgn=1;last=1
    for k in range(n-1):
        if a[k][k]==0:
            ix=next((i for i in range(k+1,n) if a[i][k]),None)
            if ix is None:return 0
            a[k],a[ix]=a[ix],a[k];sgn=-sgn
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                num=a[i][j]*pivot-a[i][k]*a[k][j]
                if num%last:raise ArithmeticError('Inexact Bareiss division')
                a[i][j]=num//last
        for i in range(k+1,n):a[i][k]=0
        last=pivot
    return sgn*a[-1][-1] if n else 1
