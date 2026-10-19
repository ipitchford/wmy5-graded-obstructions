(* Supplemental exact CAS searches. Not required by verify.py. *)
ClearAll["Global`*"];
RunParallel[a0_,b0_,cross_]:=Module[{tt,eta,theta,hs,mat,cp,eqs,pars,ans},
 tt=If[cross,t,0];eta=-2a0-2c-g-aa;theta=-3b0-d-2h-bb;
 pars=Join[{p0,p1,p2,p3,r0,r1,r2,r3,c,d,e,f,g,h,aa,bb,dd,ee},If[cross,{t},{}]];
 hs={a0 x^2 y+b0 x^3 z,a+tt x b+x(c y^2+d x y z+e x^2 z^2),b+f y^2+g x y z+h x^2 z^2,
 p0 x^2 y^3+p1 x^3 y^2 z+p2 x^4 y z^2+p3 x^5 z^3+(aa x y+bb x^2 z)a+((-theta-tt eta-a0 tt)x^2 y+dd x^3 z)b,
 r0 x y^3+r1 x^2 y^2 z+r2 x^3 y z^2+r3 x^4 z^3+(ee y+(-aa-tt ee)x z)a+(eta x y+theta x^2 z)b};
 mat=Table[D[hs[[i]],{x,y,z,a,b}[[j]]],{i,5},{j,5}]/.x->1;
 cp=Expand[CharacteristicPolynomial[mat,ll]-ll^5];
 eqs=DeleteDuplicates[Last/@CoefficientRules[cp,{ll,y,z,a,b}]];
 ans=GroebnerBasis[eqs,pars,MonomialOrder->DegreeReverseLexicographic];
 Print[{"parallel",a0,b0,cross,ans}];ans];
RunChain[]:=Module[{m,n,hs,mat,cp,eqs,pars,ans},
 pars={p0,p1,p2,p3,v0,v1,v2,v3,v4,w0,w1,w2,j,c,d,e,f,g,h,r,s,t};
 m=-2-2c-g-r;n=-3-d-2h-s;
 hs={x^2 y+x^3 z,t x a+x(c y^2+d x y z+e x^2 z^2),a+f y^2+g x y z+h x^2 z^2,
 b+p0 x y^3+p1 x^2 y^2 z+p2 x^3 y z^2+p3 x^4 z^3+(r x y+s x^2 z)a,
 v0 x^2 y^4+v1 x^3 y^3 z+v2 x^4 y^2 z^2+v3 x^5 y z^3+v4 x^6 z^4+
 (w0 x^2 y^2+w1 x^3 y z+w2 x^4 z^2)a+j x^2 a^2+(m x y+n x^2 z)b};
 mat=Table[D[hs[[i]],{x,y,z,a,b}[[j]]],{i,5},{j,5}]/.x->1;
 cp=Expand[CharacteristicPolynomial[mat,ll]-ll^5];eqs=DeleteDuplicates[Last/@CoefficientRules[cp,{ll,y,z,a,b}]];
 ans=GroebnerBasis[eqs,pars,MonomialOrder->DegreeReverseLexicographic];Print[{"chain",ans}];ans];
RunParallel[1,1,True];
RunParallel[1,0,False];
RunParallel[0,1,False];
RunChain[];
