"""Explicit polynomials audited in this research note."""
from poly import Ring

def seed_terms(R):
    x,y,z=[R.var(n) for n in ('x','y','z')];q=y-6*x*z
    f=x*x*(3*q+2*x*z)
    return q,f,24*x*q*(q+x*z),2*q*(2*q+3*x*z),12*q*q*f,2*x*q*q*(7*q+6*x*z),4*q**3*f

def six():
    R=Ring(('x','y','z','a','b','c','L'));x,y,z,a,b,c,L=R.vars()
    q,f,h2,h3,h4,h5,h6=seed_terms(R)
    return R,[f,a-h2,b-h3,-h4,c-h5,-h6]

def naive_five():
    R=Ring(('x','y','z','a','b','L'));x,y,z,a,b,L=R.vars()
    q,f,h2,h3,h4,h5,h6=seed_terms(R)
    return R,[f,3*a-h2,b-h3,-4*q*q*f,q*a-h5]

def cross():
    pars=('al','be','de','ep','t')+tuple('p'+str(i) for i in range(4))+tuple('r'+str(i) for i in range(4))
    R=Ring(('x','y','z','a','b','L')+pars);v={n:R.var(n) for n in R.names}
    x,y,z,a,b=[v[n] for n in ('x','y','z','a','b')];al,be,de,ep,t=[v[n] for n in pars[:5]]
    q,f,h2,h3,_,_,_=seed_terms(R)
    P=sum(v['p'+str(k)]*x**(2+k)*y**(3-k)*z**k for k in range(4))
    Q=sum(v['r'+str(k)]*x**(1+k)*y**(3-k)*z**k for k in range(4))
    A=al*x*y+be*x*x*z
    B=(be+t*al-3*t)*x*x*y+de*x**3*z
    C=ep*y-(al+t*ep)*x*z
    H=[f,a+t*x*b-h2,b-h3,-P+A*a+B*b,-Q+C*a-A*b]
    return R,H,pars

def chain():
    pars=('t','r','ss','j')+tuple('p'+str(i) for i in range(4))+tuple('v'+str(i) for i in range(5))+tuple('w'+str(i) for i in range(3))
    R=Ring(('x','y','z','a','b','L')+pars);v={n:R.var(n) for n in R.names}
    x,y,z,a,b=[v[n] for n in ('x','y','z','a','b')];t,r,ss,j=[v[n] for n in pars[:4]]
    q,f,h2,h3,_,_,_=seed_terms(R)
    P=sum(v['p'+str(k)]*x**(1+k)*y**(3-k)*z**k for k in range(4))
    Q=sum(v['v'+str(k)]*x**(2+k)*y**(4-k)*z**k for k in range(5))
    S=sum(v['w'+str(k)]*x**(2+k)*y**(2-k)*z**k for k in range(3))
    T=r*x*y+ss*x*x*z
    H=[f,t*x*a-h2,a-h3,b+P+T*a,Q+S*a+j*x*x*a*a-T*b]
    return R,H,pars

def hessian_potential():
    R=Ring(('x','y','a','b','c'));x,y,a,b,c=R.vars();u=1+x*y
    A=a*u**3+3*x*b*u**2-x**3*c
    B=a*y*y*u*(4+3*x*y)+b*(y+3*x*y*y*(4+3*x*y))+c*(2*x-3*x*x*y)
    return R,A*A+13*A+2*B

def cross_unreduced():
    """The full 17-parameter parallel family, before trace reductions."""
    pars=('al','be','ga','de','ep','ze','et1','et2','t')+tuple('p'+str(i) for i in range(4))+tuple('r'+str(i) for i in range(4))
    R=Ring(('x','y','z','a','b','L')+pars);v={n:R.var(n) for n in R.names}
    x,y,z,a,b=[v[n] for n in ('x','y','z','a','b')]
    al,be,ga,de,ep,ze,et1,et2,t=[v[n] for n in pars[:9]]
    _,f,g2,g3,_,_,_=seed_terms(R)
    P=sum(v['p'+str(k)]*x**(2+k)*y**(3-k)*z**k for k in range(4))
    Q=sum(v['r'+str(k)]*x**(1+k)*y**(3-k)*z**k for k in range(4))
    A=al*x*y+be*x*x*z;B=ga*x*x*y+de*x**3*z
    C=ep*y+ze*x*z;D=et1*x*y+et2*x*x*z
    return R,[f,a+t*x*b-g2,b-g3,-P+A*a+B*b,-Q+C*a+D*b],pars
