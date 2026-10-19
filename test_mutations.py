#!/usr/bin/env python3
"""Deliberate corruptions must be rejected; no files are modified."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import copy,json
from verify import cross_sources,load,check_certificate,require

def main():
    pure,final,R0,Rf,_=cross_sources();al,be,de,ep,t=R0.vars()
    tests=[]
    c=copy.deepcopy(load('cross_pure_certificate.json'));c['inputs'][0][0][1]+=1
    tests.append(('input_not_from_jacobian',c,pure,[be**2,ep**2]))
    c=copy.deepcopy(load('cross_pure_certificate.json'));c['steps'][0]['terms'][0][2]+=1
    tests.append(('false_intermediate_identity',c,pure,[be**2,ep**2]))
    c=copy.deepcopy(load('cross_final_certificate.json'));c['targets'][0]['target'][0][1]+=1
    tests.append(('changed_unit_target',c,final,[Rf.const(1)]))
    c=copy.deepcopy(load('cross_final_certificate.json'));c['targets'][0]['terms'][0][2]+=1
    tests.append(('false_final_derivation',c,final,[Rf.const(1)]))
    out=[]
    for name,cert,source,target in tests:
        try:check_certificate(cert,source,target)
        except (AssertionError,ValueError,ArithmeticError) as ex:
            out.append({'test':name,'rejected':True,'reason':str(ex)});print('REJECTED:',name)
        else:raise AssertionError('Mutation accepted: '+name)
    Path('mutation_report.json').write_text(json.dumps({'tests':out,'all_rejected':True},indent=2)+'\n')
    print('ALL CORRUPTIONS REJECTED')
if __name__=='__main__':main()
