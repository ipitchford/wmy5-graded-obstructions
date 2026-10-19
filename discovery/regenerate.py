#!/usr/bin/env python3
"""Recompute the three proof certificates with SymPy. Not needed to verify."""
from pathlib import Path
import sympy as s,json
from cert_gb import discover
ROOT=Path(__file__).resolve().parents[1]
x,y,z,a,b,L=s.symbols('x y z a b L');al,be,de,ep,t=s.symbols('al be de ep t');p=s.symbols('p0:4');r=s.symbols('r0:4')
q=y-6*x*z;f=x*x*(3*q+2*x*z);g2=24*x*q*(q+x*z);g3=2*q*(2*q+3*x*z)
P=sum(p[k]*x**(2+k)*y**(3-k)*z**k for k in range(4));Q=sum(r[k]*x**(1+k)*y**(3-k)*z**k for k in range(4))
A=al*x*y+be*x*x*z;B=(be+t*al-3*t)*x*x*y+de*x**3*z;C=ep*y-(al+t*ep)*x*z
H=s.Matrix([f,a+t*x*b-g2,b-g3,-P+A*a+B*b,-Q+C*a-A*b])
def equations(H):
 cp=s.expand(H.jacobian([x,y,z,a,b]).subs(x,1).charpoly(L).as_expr())
 es=[]
 for i in range(3,-1,-1):es.extend(s.Poly(cp.coeff(L,i),y,z,a,b).coeffs())
 return list(dict.fromkeys(s.factor(e) for e in es if e!=0))
def save(name,c):
 (ROOT/'certificates'/name).write_text(json.dumps(c,separators=(',',':')))
eq=equations(H);pure=[e for e in eq if e.free_symbols.isdisjoint(p+r)]
save('cross_pure_certificate.json',discover(pure,(al,be,de,ep,t),(be**2,ep**2),120))
es=list(dict.fromkeys(s.factor(e.subs({al:0,be:0,ep:0,de:16*t})) for e in eq));es=[e for e in es if e!=0]
save('cross_final_certificate.json',discover(es,p+r+(t,),(s.Integer(1),),120))
r0,ss,j=s.symbols('r ss j');v=s.symbols('v0:5');w=s.symbols('w0:3')
R=r0*x*y+ss*x*x*z
P=sum(p[k]*x**(1+k)*y**(3-k)*z**k for k in range(4));Q=sum(v[k]*x**(2+k)*y**(4-k)*z**k for k in range(5));S=sum(w[k]*x**(2+k)*y**(2-k)*z**k for k in range(3))
H=s.Matrix([f,t*x*a-g2,a-g3,b+P+R*a,Q+S*a+j*x*x*a*a-R*b])
es=list(dict.fromkeys(s.factor(e.subs({t:0,r0:0,ss:0,j:0})) for e in equations(H)));es=[e for e in es if e!=0]
save('chain_final_certificate.json',discover(es,p+v+w,(s.Integer(1),),120))
print('Certificates regenerated. Run the standalone verifier; release hashes no longer apply after edits.')
