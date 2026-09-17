#!/usr/bin/env python3
"""Reproduce the repaired, readable core; never overwrite the supplied source.

These are newly written repairs and readability changes, not recovered author
code. --check-output compares executable tokens without changing any files or
requiring a formatter. --output-source regenerates the production copy and two
reviewable diffs, using the recorded clang-format configuration.
"""
from pathlib import Path
import argparse
import difflib
from prepare_readable_core import readable, format_cpp, code_tokens


PATCH_DIR = Path(__file__).resolve().parent
ROOT_DIR = PATCH_DIR.parent
DEFAULT_SOURCE = ROOT_DIR / "reference/verified_core_taylor_model.cpp"


def replace_once(text, old, new, description):
    count = text.count(old)
    if count != 1:
        raise ValueError(f"{description}: expected one source match, found {count}")
    return text.replace(old, new, 1)


def patched(text):
    # Both the historical reference and the portable source copy are accepted.
    for header in ('mpfr_interval.hpp','fixed_weight_schedule.hpp'):
        absolute=f'#include "/mnt/data/ew_work2/{header}"'
        portable=f'#include "{header}"'
        if absolute in text:
            text=replace_once(text,absolute,portable,f'portable include for {header}')
        elif text.count(portable)!=1:
            raise ValueError(f'Expected exactly one include for {header}')
    changes = [
        ('// TM_INSERT_MARKER\n','','remove obsolete development marker'),
        (
            'QI r;if(b.lo>0)r=QI(qdn(1/b.hi),qup(1/b.lo),true);else r=QI(qdn(1/b.lo),qup(1/b.hi),true);return a*r;',
            'QI r(qdn(1/b.hi),qup(1/b.lo),true);return a*r;',
            'correct reciprocal endpoints on negative intervals',
        ),
        (
            'static __float128 qsumup(__float128 a,__float128 b){return qup(a+b);}',
            '''static __float128 qsumup(__float128 a,__float128 b){return qup(a+b);}
// Enclose distances from the actual stored center, including subtraction rounding.
static __float128 qradius_about(QI a,__float128 center){
 return fmaxq(qup(center-a.lo),qup(a.hi-center));
}
// These numbers are subsequently parsed as certified lower endpoints.
static std::string mdownstr(const Real&a,int digits){
 std::vector<char> out((size_t)digits+128);
 int n=mpfr_snprintf(out.data(),out.size(),"%.*RDg",digits,a.x);
 if(n<0)throw std::runtime_error("MPFR downward formatting failed");
 if((size_t)n>=out.size()){
  out.resize((size_t)n+1);
  mpfr_snprintf(out.data(),out.size(),"%.*RDg",digits,a.x);
 }
 return std::string(out.data());
}''',
            'introduce directed radius and lower-endpoint formatting helpers',
        ),
        (
            'QI v=sp[i][a][K-1];for(int n=K-2;n>=0;n--)v=v*tau+sp[i][a][n];Pr[i][a]=v;',
            'QI v=si[i][a][K-1];for(int n=K-2;n>=0;n--)v=v*tau+si[i][a][n];Pr[i][a]=v;',
            'use actual-Jacobian Taylor enclosure in sensitivity tube',
        ),
        (
            'if(need.lo<St[i][a].lo||need.hi>St[i][a].hi)stable=false;',
            'if(!qstrict_subset(need,St[i][a]))stable=false;',
            'require strict sensitivity-tube inclusion',
        ),
        (
            'static QI qexpand(QI x,__float128 e)',
            'static bool qstrict_subset(QI inner,QI outer){return inner.lo>outer.lo && inner.hi<outer.hi;}\nstatic QI qexpand(QI x,__float128 e)',
            'name the strict binary128 containment test for production and tests',
        ),
        (
            '__float128 l=fminq(St[i][a].lo,need.lo),u=fmaxq(St[i][a].hi,need.hi),rr=qup((u-l)/2*1.02Q+1e-80Q),cc=(u+l)/2;St[i][a]=QI(qdn(cc-rr),qup(cc+rr),true);',
            '__float128 l=fminq(St[i][a].lo,need.lo),u=fmaxq(St[i][a].hi,need.hi),cc=(u+l)/2;cc=fminq(u,fmaxq(l,cc));__float128 rr=qsumup(qmulup(qradius_about(QI(l,u,true),cc),qup(1.02Q)),qup(1e-80Q));St[i][a]=QI(qdn(cc-rr),qup(cc+rr),true);',
            'make sensitivity tube inflation contain its hull about rounded center',
        ),
        (
            'qdivup(fmaxq(fabsq(vi.lo-pc),fabsq(vi.hi-pc)),qlo[i])',
            'qdivup(qradius_about(vi,pc),qlo[i])',
            'round each endpoint distance outward in local sensitivity error',
        ),
        (
            'qdivup(qmulup(0.5Q,qup(aa.hi-aa.lo)),ql[i])',
            'qdivup(qradius_about(aa,S[i][0]),ql[i])',
            'initialize first sensitivity radius about stored midpoint',
        ),
        (
            'qdivup(qmulup(0.5Q,qup(bb.hi-bb.lo)),ql[i])',
            'qdivup(qradius_about(bb,S[i][1]),ql[i])',
            'initialize second sensitivity radius about stored midpoint',
        ),
        (
            'str(minA,30)<<" "<<str(minD,30)<<" "<<str(minU,30)',
            'mdownstr(minA,30)<<" "<<mdownstr(minD,30)<<" "<<mdownstr(minU,30)',
            'publish MPFR positivity minima with rounding toward minus infinity',
        ),
    ]
    for old, new, description in changes:
        text = replace_once(text, old, new, description)
    return text


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source', nargs='?', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--output-source', type=Path)
    parser.add_argument('--check-output', type=Path,
                        help='Check generated C++ tokens without requiring clang-format')
    parser.add_argument('--clang-format', default=None,
                        help='Formatter executable for --output-source maintenance runs')
    args = parser.parse_args()
    source = args.source.resolve()
    original = source.read_text()
    revised = patched(original)
    annotated = readable(revised)
    if args.check_output:
        dest = args.check_output.resolve()
        if code_tokens(dest.read_text()) != code_tokens(annotated):
            raise ValueError('Production source differs from the reproducible repair/readability transformations; run make regenerate-core')
        print('PASS: production core tokens match the repaired reference and named readability transformations')
        return
    diff = ''.join(difflib.unified_diff(
        original.splitlines(keepends=True), revised.splitlines(keepends=True),
        fromfile=str(source.relative_to(ROOT_DIR)) if source.is_relative_to(ROOT_DIR) else source.name,
        tofile='src/verified_core_taylor_model.repaired.cpp'))
    diff_path = PATCH_DIR / 'core_validation.patch'
    diff_path.write_text(diff)
    print(f'Wrote reviewable diff: {diff_path}')
    if args.output_source:
        dest = args.output_source.resolve()
        if dest == source:
            raise ValueError('Refusing to overwrite input source')
        dest.parent.mkdir(parents=True, exist_ok=True)
        formatted = format_cpp(annotated, ROOT_DIR, args.clang_format)
        if code_tokens(formatted) != code_tokens(annotated):
            raise ValueError('Formatter changed C++ tokens; refusing to write output')
        dest.write_text(formatted)
        readability_diff = ''.join(difflib.unified_diff(
            revised.splitlines(keepends=True), formatted.splitlines(keepends=True),
            fromfile='repaired-before-readability.cpp',
            tofile='src/verified_core_taylor_model.repaired.cpp'))
        (PATCH_DIR/'core_readability.patch').write_text(readability_diff)
        print(f'Wrote separate repaired copy: {dest}')


if __name__ == '__main__':
    main()
