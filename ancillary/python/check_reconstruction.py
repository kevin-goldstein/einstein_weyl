#!/usr/bin/env python3
"""Run portable rational checks and inspect the reconstructed core output.

The supplied core output is a historical comparison input, not a certified
enclosure established by this program. Only the explicitly selected regenerated
core output is passed to the final matching calculation. A successful run does
not replace independent review of the generating C++ validation algorithm.
"""
import argparse
from decimal import Decimal, localcontext, ROUND_CEILING
from fractions import Fraction
from pathlib import Path
import os
import re
import subprocess
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parent
NAMES=('A','p','C','d','Y','J')

def parse_centers(path):
    text=path.read_text()
    centers={}
    for name in NAMES:
        match=re.search(rf'^{name} center \[([^,]+),([^\]]+)\]',text,re.M)
        if match is None:
            raise ValueError(f'Missing center interval for {name} in {path}')
        lo,hi=Fraction(match[1]),Fraction(match[2])
        if lo>hi:
            raise ValueError(f'Reversed {name} interval in {path}')
        centers[name]=(lo,hi)
    match=re.search(r'^t \[([^,]+),([^\]]+)\]',text,re.M)
    if match is None or not Fraction(match[1])<=40<=Fraction(match[2]):
        raise ValueError(f'{path} is not an endpoint output at t=40')
    return centers

def number(q):
    if q == 0:
        return '0'
    with localcontext() as context:
        context.prec=25
        return f'{Decimal(q.numerator)/Decimal(q.denominator):.16E}'

def upper_number(q):
    with localcontext() as context:
        context.prec=17
        context.rounding=ROUND_CEILING
        return f'{Decimal(q.numerator)/Decimal(q.denominator):.16E}'

def compare(core,reference):
    new,old=parse_centers(core),parse_centers(reference)
    print('\nCore center comparison (diagnostic only):',flush=True)
    print('  Regenerated input:',core,flush=True)
    print('  Supplied historical reference:',reference,flush=True)
    print('  The reference intervals are not established as rigorous by this comparison.',flush=True)
    for name in NAMES:
        lo,hi=new[name];ol,oh=old[name]
        relation=('inside supplied interval' if ol<=lo<=hi<=oh else
                  'contains supplied interval' if lo<=ol<=oh<=hi else
                  'overlaps supplied interval' if max(lo,ol)<=min(hi,oh) else
                  'DISJOINT from supplied interval')
        delta=(lo+hi-ol-oh)/2
        ratio=(hi-lo)/(oh-ol) if oh>ol else None
        print(f'  {name}: {relation}; midpoint difference={number(delta)}; '
              f'new width={number(hi-lo)}; '
              f'width ratio={number(ratio) if ratio is not None else "undefined"}',flush=True)

def check_initial_second_variations(horizon,weights):
    """Check all initial scaled Hessian norms against the core's Binit.

    The small relative reduction below conservatively encloses the weight's
    binary128 decimal conversion followed by one downward nextafter operation.
    For the checked normal positive weights, this is much larger than the
    binary128 rounding allowance. All calculations here are exact rationals.
    """
    header=weights.read_text()
    first=re.search(r'EW_Q\[EW_QN\]\[6\]\s*=\s*\{\s*\{([^}]+)\}',header)
    if first is None:
        raise ValueError(f'Cannot read first weight vector from {weights}')
    q=[Fraction(x) for x in re.findall(r'"([^"]+)"',first[1])]
    if len(q)!=6 or not all(Fraction(1,2**1000)<x<2**1000 for x in q):
        raise ValueError('Initial weights must be six normal positive binary128-scale values')
    qlo=[x*(1-Fraction(1,2**100)) for x in q]
    text=horizon.read_text()
    bounds={key:[] for key in ('00','01','11')}
    for i,name in enumerate(NAMES):
        block=re.search(rf'^{name}\n(.*?)(?=^[A-Za-z]+\n|\Z)',text,re.M|re.S)
        if block is None:
            raise ValueError(f'Missing Hessian block {name} in {horizon}')
        for key in bounds:
            match=re.search(rf'^ h{key} \[([^,]+),([^\]]+)\]',block[1],re.M)
            if match is None:
                raise ValueError(f'Missing h{key} interval for {name}')
            lo,hi=Fraction(match[1]),Fraction(match[2])
            if lo>hi:
                raise ValueError(f'Reversed h{key} interval for {name}')
            bounds[key].append(max(abs(lo),abs(hi))/qlo[i])
    binit=Fraction('1e-50')
    print('\nInitial weighted second-variation bound (exact rational check):',flush=True)
    for key,components in bounds.items():
        norm=max(components)
        print(f'  h{key}: norm <= {upper_number(norm)}; norm / Binit <= {upper_number(norm/binit)}',flush=True)
        if not norm<binit:
            raise ValueError(f'Initial h{key} weighted norm is not bounded by Binit=1e-50')
    print('  PASS: Binit=1e-50 bounds all three initial weighted Hessian norms.',flush=True)

def run(script,*arguments):
    print('\nRunning:',script.relative_to(HERE),*[str(a) for a in arguments],flush=True)
    environment=os.environ.copy()
    environment['PYTHONDONTWRITEBYTECODE']='1'
    result=subprocess.run([sys.executable,str(script),*[str(a) for a in arguments]],
                          cwd=ROOT,env=environment)
    if result.returncode:
        raise RuntimeError(f'{script.name} failed with exit status {result.returncode}')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--core',type=Path,default=ROOT/'logs/core_reconstructed.txt',
                        help='Regenerated t=40 core output to pass explicitly to matching')
    parser.add_argument('--reference',type=Path,default=ROOT/'reference/core_tm_final2_t40.txt',
                        help='Supplied output used only for diagnostic comparison')
    parser.add_argument('--matching-only',action='store_true',
                        help='Skip the independent horizon, tail, and constraint Python checks')
    parser.add_argument('--horizon-second',type=Path,default=ROOT/'logs/horizon_second_derivative.txt',
                        help='Regenerated horizon second-derivative output')
    parser.add_argument('--weights',type=Path,default=ROOT/'include/fixed_weight_schedule.hpp',
                        help='Replacement weight schedule used by the reconstructed core')
    args=parser.parse_args()
    core=args.core.resolve();reference=args.reference.resolve()
    if not core.is_file():
        raise FileNotFoundError(f'Regenerated core output is missing: {core}. Run the rebuilt core first.')
    if core==reference:
        print('NOTICE: core input equals the supplied reference. This is a script smoke test, '
              'not validation of a reconstructed core.',flush=True)
    compare(core,reference)
    check_initial_second_variations(args.horizon_second.resolve(),args.weights.resolve())
    if not args.matching_only:
        for relative in ('horizon/horizon_majorant_certificate.py',
                         'horizon/horizon_tail_bounds.py',
                         'horizon/horizon_derivation_checks.py',
                         'horizon/horizon_initialization_bounds.py',
                         'tail/centered_tail_uniform_certificate.py',
                         'matching/constraint_v_uniqueness.py'):
            run(HERE/relative)
    run(HERE/'matching/poincare_certificate.py',core)
    print('\nPASS: the selected Python checks and exact-rational matching calculation completed.',flush=True)
    print('The generating C++ validation and its assumptions remain part of the proof obligation.',flush=True)

if __name__=='__main__':
    try:main()
    except (FileNotFoundError,ValueError,RuntimeError) as error:
        print(f'ERROR: {error}',file=sys.stderr)
        raise SystemExit(1)
