"""Exact Buchberger discovery with polynomial-identity proof steps."""
from sympy.polys.rings import ring
from sympy.polys.domains import QQ
from sympy.polys.orderings import grevlex
import sympy as sp, heapq, time, json

def encoded(poly):
 return [[list(m),int(c.numerator),int(c.denominator)] for m,c in sorted(poly.items())]
def discover(expressions,variables,targets,limit=35):
 R,*_=ring(','.join(map(str,variables)), QQ, grevlex)
 polys=[R.from_expr(sp.expand(e)) for e in expressions];steps=[];ninputs=len(polys);basis=[];pairs=[];done=set();start=time.time();zero=(0,)*len(variables)
 def monsub(u,v):
  w=tuple(a-b for a,b in zip(u,v));return w if min(w)>=0 else None
 def reducep(p,combo):
  rem=R.zero
  while p:
   m,c=p.LT
   for ix in basis:
    g=polys[ix];d=monsub(m,g.LM)
    if d is not None:
     a=c/g.LC;p-=g.mul_term((d,a));combo.append([ix,list(d),-int(a.numerator),int(a.denominator)]);break
   else:
    q=R.term_new(m,c);rem+=q;p-=q
  return rem,combo
 def appendp(p,combo):
  lc=p.LC;p=p.quo_ground(lc) if hasattr(p,'quo_ground') else p/lc
  combo=[[i,m,n*int(lc.denominator),d*int(lc.numerator)] for i,m,n,d in combo]
  ix=len(polys);polys.append(p);steps.append({'terms':combo,'poly':encoded(p)})
  for j in basis:
   u,v=p.LM,polys[j].LM
   if all(min(a,b)==0 for a,b in zip(u,v)): continue
   l=tuple(max(a,b) for a,b in zip(u,v));heapq.heappush(pairs,(sum(l),grevlex(l),j,ix))
  basis.append(ix)
  return ix
 for i in range(ninputs):
  p,comb=reducep(polys[i],[[i,list(zero),1,1]])
  if p:appendp(p,comb)
 while pairs:
  if time.time()-start>limit:raise TimeoutError((len(polys),len(pairs)))
  _,_,i,j=heapq.heappop(pairs);key=(i,j)
  if key in done:continue
  done.add(key);u,v=polys[i],polys[j];l=tuple(max(a,b) for a,b in zip(u.LM,v.LM));di=monsub(l,u.LM);dj=monsub(l,v.LM)
  p=u.mul_term((di,1/u.LC))-v.mul_term((dj,1/v.LC))
  combo=[[i,list(di),int((1/u.LC).numerator),int((1/u.LC).denominator)],[j,list(dj),-int((1/v.LC).numerator),int((1/v.LC).denominator)]]
  p,combo=reducep(p,combo)
  if p:
   ix=appendp(p,combo)
   if polys[ix]==R.one:break
 targetproof=[]
 for tar in targets:
  p=R.from_expr(tar);rem,comb=reducep(p,[])
  if rem:raise ValueError(('target not proved',str(tar),str(rem)))
  targetproof.append({'target':encoded(p),'terms':[[i,m,-n,d] for i,m,n,d in comb]})
 print('certificate',len(polys),'polys',len(steps),'steps',len(pairs),'pairs left','secs',time.time()-start,flush=True)
 return {'variables':list(map(str,variables)),'inputs':[encoded(p) for p in polys[:ninputs]],'steps':steps,'targets':targetproof}

