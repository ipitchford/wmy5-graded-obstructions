#!/usr/bin/env python3
"""Replay exact algebraic certificates without SymPy or a CAS.

Run: python3 -I verify.py --report verification_report.json
The -I option requires adding this script's directory explicitly below.
"""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import argparse,json,time,hashlib
from fractions import Fraction
from poly import Ring,Poly,charpoly,coefficients,unique,decode,bareiss
import models
ROOT=Path(__file__).resolve().parent

def require(ok,msg):
    if not ok:raise AssertionError(msg)

def load(name):return json.loads((ROOT/'certificates'/name).read_text())

def check_certificate(cert,source,expected_targets):
    R=Ring(cert['variables']);source=[p.reindex(R) for p in source]
    inp=[decode(R,t) for t in cert['inputs']]
    require({p.key() for p in inp}=={p.key() for p in source},'Certificate inputs do not match derived coefficient equations')
    ps=inp[:]
    def combination(terms):
        out=R.const(0)
        for ix,m,n,d in terms:
            require(isinstance(ix,int) and 0<=ix<len(ps),'Invalid proof reference')
            require(len(m)==R.n and all(isinstance(k,int) and k>=0 for k in m),'Invalid multiplier')
            require(isinstance(n,int) and isinstance(d,int) and d!=0,'Invalid proof rational')
            out=out+ps[ix].mul_term(tuple(m),Fraction(n,d))
        return out
    for k,step in enumerate(cert['steps']):
        got=combination(step['terms']);want=decode(R,step['poly'])
        require(got==want,'Invalid polynomial identity at proof step '+str(k));ps.append(want)
    require(len(cert['targets'])==len(expected_targets),'Wrong number of targets')
    for item,expected in zip(cert['targets'],expected_targets):
        expected=expected.reindex(R);want=decode(R,item['target'])
        require(want==expected,'Target was changed')
        require(combination(item['terms'])==expected,'Target not derived')
    return {'inputs':len(inp),'identity_steps':len(cert['steps']),'targets':len(expected_targets)}

def check_six():
    R,H=models.six();L=R.var('L');vs=('x','y','z','a','b','c')
    residual=charpoly(H,vs)-L**6;require(not residual,'Six-dimensional characteristic polynomial')
    pts=[(0,1,-4,0,0,0),(-1,-11,2,12,-6,4),(1,7,Fraction(5,4),-3,Fraction(-3,2),Fraction(1,2))]
    for pt in pts:
        val=dict(zip(vs,pt))
        require(all(H[i].value(val)-pt[i]+(1 if i==1 else 0)==0 for i in range(6)),'Six-dimensional zero')
    return {'characteristic_residual_terms':len(residual.d),'checked_zeros':len(pts)}

def check_naive():
    R,H=models.naive_five();L=R.var('L');vs=('x','y','z','a','b')
    cp=charpoly(H,vs)
    for name,v in zip(vs,(1,0,0,-1,0)):cp=cp.set(name,v)
    require(cp==L**3*(L**2-6),'Naive compression witness polynomial')
    pts=[(0,1,-4,0,0),(-1,-11,2,4,-6),(1,7,Fraction(5,4),-1,Fraction(-3,2))]
    for pt in pts:
        val=dict(zip(vs,pt));require(all(H[i].value(val)-pt[i]+(1 if i==1 else 0)==0 for i in range(5)),'Naive compressed zero')
    return {'preserved_zeros':3,'witness_point':[1,0,0,-1,0],'DH_characteristic':'L^3*(L^2-6)','DX_positive_eigenvalue':'-1+sqrt(6)'}

def check_cross_reduction():
    R,H,pars=models.cross_unreduced();v={n:R.var(n) for n in R.names}
    x,y,z,a,b=[v[n] for n in ('x','y','z','a','b')]
    al,be,ga,ep,ze,et1,et2,t=[v[n] for n in ('al','be','ga','ep','ze','et1','et2','t')]
    vs=('x','y','z','a','b')
    trace=sum(H[i].diff(vs[i]) for i in range(5))
    require(trace==(al+et1)*x*y+(be+et2)*x*x*z,'Unreduced trace identity')
    ht=[p.subst(R,{'et1':-al,'et2':-be}) for p in H]
    J=[[p.diff(name) for name in vs] for p in ht]
    tr=sum(J[i][i] for i in range(5))
    e2=(tr*tr-sum(J[i][j]*J[j][i] for i in range(5) for j in range(5)))*Fraction(1,2)
    require(e2.diff('a')==-x*(al+t*ep+ze),'First auxiliary relation in second coefficient')
    require(e2.diff('b')==-x*x*(3*t+ga-t*al-be),'Second auxiliary relation in second coefficient')
    RR,HH,_=models.cross();vv={n:RR.var(n) for n in RR.names}
    images={'et1':-vv['al'],'et2':-vv['be'],'ze':-vv['al']-vv['t']*vv['ep'],'ga':vv['be']+vv['t']*vv['al']-3*vv['t']}
    require(all(p.subst(RR,images)==q for p,q in zip(H,HH)),'Unreduced-to-certificate-model bridge')
    return {'unreduced_parameters':len(pars),'trace_identity':True,'second_coefficient_relations':2,'substituted_model_match':True}

def cross_sources():
    R,H,pars=models.cross();cp=charpoly(H,('x','y','z','a','b'),specialise_x=True)-R.var('L')**5
    eq=coefficients(cp,('L','y','z','a','b'),pars)
    pure=[]
    for p in eq:
        if all(not any(m[5:]) for m in p.d):pure.append(p.reindex(Ring(pars[:5])))
    R0=Ring(pars[:5]);al,be,de,ep,t=R0.vars()
    require((-al**2-al*ep*t-be*ep).key() in {p.key() for p in pure},'Missing alpha-square coefficient')
    reduced=[p.subst(Ring(('de','t')),{'al':0,'be':0,'ep':0}) for p in pure]
    D,T=Ring(('de','t')).vars()
    require((8*D-128*T).key() in {p.key() for p in reduced},'Missing delta relation')
    Rf=Ring(pars[5:]+('t',));images={'al':0,'be':0,'ep':0,'de':16*Rf.var('t')}
    final=unique([p.subst(Rf,images) for p in eq]);final=[p for p in final if p]
    return pure,final,R0,Rf,len(eq)

def check_cross():
    bridge=check_cross_reduction()
    pure,final,R0,Rf,n=cross_sources();al,be,de,ep,t=R0.vars()
    c1=check_certificate(load('cross_pure_certificate.json'),pure,[be**2,ep**2])
    # Over a characteristic-zero field, be^2=ep^2=0 implies be=ep=0.
    # The explicitly checked coefficient then gives al^2=0, so al=0.
    # The checked linear relation gives de=16*t. These exhaust all solutions.
    c2=check_certificate(load('cross_final_certificate.json'),final,[Rf.const(1)])
    return {'reduction_bridge':bridge,'nonzero_coefficient_equations':n,'square_certificate':c1,'specialised_unit_certificate':c2,'conclusion':'no parameters over any characteristic-zero field'}

def chain_sources():
    R,H,pars=models.chain();cp=charpoly(H,('x','y','z','a','b'),specialise_x=True)-R.var('L')**5
    eq=coefficients(cp,('L','y','z','a','b'),pars);R0=Ring(pars[:4]);t,r,ss,j=R0.vars()
    pure=[p.reindex(R0) for p in eq if all(not any(m[4:]) for m in p.d)]
    needed=[-2*j-r*t-ss-3*t,r*t+ss,-45*r*t+248*r-8*ss*t+42*ss,32*j+16*r*t+3*ss*t,-t*(16*r+3*ss)]
    keys={p.key() for p in pure}
    require(all(p.key() in keys for p in needed),'Missing chain reduction equations')
    require(needed[3]+needed[4]==32*j,'Chain proof j=0')
    require(needed[0]+needed[1]==-2*j-3*t,'Chain proof t=0')
    # These imply j=t=ss=0; the third then implies r=0.
    Rf=Ring(pars[4:]);final=unique([p.subst(Rf,{'t':0,'r':0,'ss':0,'j':0}) for p in eq]);final=[p for p in final if p]
    return final,Rf,len(eq)

def check_chain():
    source,R,n=chain_sources();proof=check_certificate(load('chain_final_certificate.json'),source,[R.const(1)])
    return {'nonzero_coefficient_equations':n,'unit_certificate':proof,'conclusion':'no parameters over any characteristic-zero field'}

def check_hessian():
    cert=load('hessian_trace_certificate.json');R,Psi=models.hessian_potential();pairs=[(i,j) for i in range(5) for j in range(i,5)]
    require(cert['pairs']==[list(p) for p in pairs],'Hessian index list')
    hs=[Psi.diff(R.names[i]).diff(R.names[j])*(1 if i==j else 2) for i,j in pairs]
    mons=[tuple(m) for m in cert['monomials']]
    require(len(mons)==15 and all(any(m) for m in mons),'Need 15 nonconstant coefficient rows')
    M=[[int(p.d.get(m,0)) for p in hs] for m in mons]
    require(M==cert['matrix'],'Hessian coefficient matrix does not match polynomial')
    det=bareiss(M);require(det!=0 and str(det)==cert['determinant'],'Hessian full-rank certificate')
    return {'rank':15,'minor_determinant':str(det),'conclusion':'constant trace of B*Hess(Psi) forces symmetric part of B to vanish, hence trace zero'}

def main(report=None):
    t=time.time();res={'status':'partial; no five-dimensional counterexample constructed','arithmetic':'exact integers and fractions; standard-library implementation','checks':{}}
    for name,fn in [('six_dimensional_baseline',check_six),('naive_compression_rejected',check_naive),('parallel_auxiliaries',check_cross),('unequal_weight_chain',check_chain),('hessian_trace_obstruction',check_hessian)]:
        tick=time.time();r=fn();r['seconds']=round(time.time()-tick,6);res['checks'][name]=r;print('PASS:',name,flush=True)
    res['seconds']=round(time.time()-t,6)
    if report:Path(report).write_text(json.dumps(res,indent=2)+'\n')
    print('ALL EXACT CHECKS PASSED')
    return res
if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--report');args=parser.parse_args();main(args.report)
