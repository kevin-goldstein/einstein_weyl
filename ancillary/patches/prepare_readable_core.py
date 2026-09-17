#!/usr/bin/env python3
"""Deterministic readability transformations for the distributed core source.

No arithmetic expression is rewritten. The transformations remove four unused
development helpers, rename identifiers, add proof references, make main's
implicit successful return explicit, and format C++.
The token checker used by make verify does not require clang-format.
"""
from pathlib import Path
import os
import re
import shutil
import subprocess

# These names identify computational stages rather than changing their order.
NAMES = {
    "QI": "QuadInterval", "QAD2": "StateHessianDual",
    "MTS": "MpfrStateTaylor", "QTS": "QuadStateTaylor",
    "QJS": "QuadJacobianTaylor", "QSS": "QuadSensitivityTaylor",
    "MD": "HorizonDual", "MEV": "HorizonDualJets", "EV": "HorizonJets",
    "mhull": "mpfr_hull", "minflate": "mpfr_inflate",
    "msubset": "strict_mpfr_subset", "mrangeprod": "mpfr_forward_range_product",
    "rhs": "core_rhs", "mconv": "mpfr_convolution",
    "mcoef": "mpfr_state_taylor_coefficients", "mstep": "validated_center_step",
    "qconv": "quad_convolution", "qjacseries": "quad_jacobian_taylor_coefficients",
    "qsenscoef": "quad_sensitivity_taylor_coefficients",
    "qstatecoef": "quad_state_taylor_coefficients",
    "qevalstate": "evaluate_center_at_offset",
    "qstep_centeredQ_TM": "validated_parameter_substep",
    "qstrict_subset": "strict_quad_subset", "qrhs_hessian": "interval_rhs_hessian",
    "qweights": "load_outward_weights", "qindex": "weight_index",
    "meval": "evaluate_horizon_dual_series",
    "mhorizon": "initialize_horizon_first_variations",
    "evalp": "evaluate_horizon_series",
    "powu": "mpfr_positive_integer_power",
    "tailb": "horizon_derivative_tail",
    "center_horizon": "initialize_horizon_center",
}
# Retain the C++ lexical tokens used by this source; omit comments/whitespace.
TOKEN = re.compile(
    r'//[^\n]*|/\*[\s\S]*?\*/|\s+|'
    r'"(?:\\[\s\S]|[^"\\])*"|\'(?:\\[\s\S]|[^\'\\])*\'|'
    r'[A-Za-z_][A-Za-z_0-9]*|'
    r'(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eEpP][+-]?[0-9]+)?[A-Za-z_0-9]*|'
    r'>>=|<<=|->\*|::|->|\+\+|--|&&|\|\||<=|>=|==|!=|<<|>>|'
    r'\+=|-=|\*=|/=|%=|&=|\|=|\^=|\.{3}|[^\s]'
)
def code_tokens(source):
    return [m.group() for m in TOKEN.finditer(source)
            if not m.group().isspace() and not m.group().startswith(("//","/*"))]

def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"Readability anchor must occur once: {old!r}")
    return text.replace(old,new,1)

def remove_function(text, signature):
    if text.count(signature) != 1:
        raise ValueError(f"Unused helper must occur once: {signature}")
    start=text.index(signature); body=text.index("{",start); depth=0
    for m in TOKEN.finditer(text,body):
        tok=m.group()
        if tok=="{": depth+=1
        elif tok=="}":
            depth-=1
            if depth==0: return text[:start]+text[m.end():]
    raise ValueError(f"Unbalanced helper body: {signature}")

def readable(text):
    for signature in (
        "template<int K>static QTS<K>qtruncate(",
        "static unsigned long long qbinom(",
        "template<int K>static QTS<K> qshiftstate(",
        "static std::string qs(QI a,int n=38)",
    ):
        text=remove_function(text,signature)
    text=replace_once(text,
        'std::cout<<" core_quadratic_remainder "<<qs(rr,20)<<"\\n";}\n}',
        'std::cout<<" core_quadratic_remainder "<<qs(rr,20)<<"\\n";}\n return 0;\n}')
    references = (
        ("// MPFR interval helpers for the rigorously validated center.",
         "// Center-state interval arithmetic.\n"
         "// Appendix C labels: eq:core-picard-inclusion, eq:core-order56-step."),
        ("template<class T>std::array<T,6>rhs",
         "// Six-state vector field, in the order (A,p,C,d,Y,J).\n"
         "// Main equations eq:F and eq:first-order-rest; Appendix C uses this\n"
         "// same expression for state, Jacobian and Hessian calculations.\n"
         "template<class T>std::array<T,6>rhs"),
        ("template<int K>static MTS<K>mcoef",
         "// Taylor coefficients obtained directly from the ODE recurrence.\n"
         "// Appendix C: eq:core-order56-step.\n"
         "template<int K>static MTS<K>mcoef"),
        ("template<int K>static bool mstep",
         "// Strict Picard tube, followed by the order-56 Lagrange-remainder step.\n"
         "// Appendix C: eq:core-picard-inclusion, eq:core-order56-step.\n"
         "template<int K>static bool mstep"),
        ("struct QAD2 {",
         "// Interval automatic differentiation in the six state coordinates.\n"
         "// Appendix C: eq:core-hessian-source.\n"
         "struct QAD2 {"),
        ("static bool qstrict_subset",
         "// Equality at either endpoint is rejected, including for zero-width\n"
         "// candidates. This is the strict tube test used in the production step.\n"
         "static bool qstrict_subset"),
        ("template<int K>static QJS<K>qjacseries",
         "// Enclose the time coefficients of the true center-orbit Jacobian.\n"
         "// Appendix C: eq:core-first-variation, eq:core-sensitivity-recurrence.\n"
         "template<int K>static QJS<K>qjacseries"),
        ("template<int K>static QSS<K>qsenscoef",
         "// Linear variational recurrence, including interval coefficient widths.\n"
         "// Appendix C: eq:core-sensitivity-recurrence.\n"
         "template<int K>static QSS<K>qsenscoef"),
        ("template<int K>static bool qstep_centeredQ_TM",
         "// One sensitivity substep and the simultaneous parameter bootstrap.\n"
         "// Appendix C: eq:core-sensitivity-tube through eq:core-error-update;\n"
         "// eq:core-hessian-source through eq:core-bootstrap-inequalities.\n"
         "template<int K>static bool qstep_centeredQ_TM"),
        (" QJS<K>Mp;",
         " // A midpoint Jacobian is only a predictor for the next stored center.\n"
         " // The inclusion polynomial below always uses the true interval M.\n"
         " QJS<K>Mp;"),
        (" bool tube_ok=false;",
         " // Prove strict inclusion before using this tube for the remainder.\n"
         " // Appendix C: eq:core-sensitivity-tube.\n"
         " bool tube_ok=false;"),
        (" auto H=qrhs_hessian(zbox,qT);",
         " // Hessian source and logarithmic norm on the enlarged state box.\n"
         " // Appendix C: eq:core-hessian-source, eq:core-second-bound.\n"
         " auto H=qrhs_hessian(zbox,qT);"),
        (" __float128 tr=0;",
         " // Transfer every weighted bound to the next positive weight vector.\n"
         " // Appendix C: eq:core-error-update.\n"
         " __float128 tr=0;"),
        ("// MPFR first-order horizon dual, only to initialize the two scaled columns.",
         "// Horizon initialization with first derivatives in the two scaled parameters.\n"
         "// Appendix C: eq:core-initial-map, eq:initial-dual-padding."),
        ("// Center horizon series with explicit majorant tail (from the previous certificate).",
         "// Center horizon series with the explicit analytic majorant tail.\n"
         "// Appendix A: eq:general-horizon-tail; Appendix C: eq:core-initial-map."),
        ("int main(int argc,char**argv){",
         "// Execute the declared coarse/substep schedule and serialize the proof data.\n"
         "// Appendix C: eq:core-coarse-schedule, eq:core-parameter-remainder.\n"
         "int main(int argc,char**argv){"),
    )
    for old,new in references: text=replace_once(text,old,new)
    # Rename only identifier tokens, never literals or comments.
    text=TOKEN.sub(lambda m: NAMES.get(m.group(),m.group()),text)
    return ("// Distributed production core. Generated by patches/patch_core_validation.py.\n"
            "// Four unused development helpers are omitted; numerical stage order is unchanged.\n"
            "// See patches/READABILITY.md for the source map and regeneration/check commands.\n"
            +text)

def format_cpp(text, root, requested=None):
    executable=requested or os.environ.get("CLANG_FORMAT") or shutil.which("clang-format")
    if not executable:
        candidate=Path("/Library/Developer/CommandLineTools/usr/bin/clang-format")
        executable=str(candidate) if candidate.is_file() else None
    if not executable:
        raise RuntimeError("clang-format is needed only to regenerate readable C++; set CLANG_FORMAT")
    result=subprocess.run([executable,"--style=file", "--assume-filename="+str(root/"src/core.cpp")],
                          input=text,text=True,capture_output=True,check=True,cwd=root)
    print(subprocess.check_output([executable,"--version"],text=True).strip())
    return result.stdout
