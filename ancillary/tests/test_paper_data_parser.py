#!/usr/bin/env python3
"""Regression checks for exact decimal transfer from formatted C++ source."""
from fractions import Fraction
import importlib.util
from pathlib import Path

root = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location(
    'generate_paper_data', root/'python/matching/generate_paper_data.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
read = module.source_decimal

for source in ('I b0("0.125");', 'I b0(\n    "0.125"\n);',
               'I b0 \t( "0.125" );', 'I other("9"), b0("0.125");'):
    assert read(source, 'b0') == Fraction(1, 8)
assert read('I R0(\n"-1.25e-3"\n);', 'R0') == Fraction(-1, 800)
for source in ('I other("0.125");', 'I b0("0.125"); I b0("0.25");'):
    try:
        read(source, 'b0')
    except ValueError:
        pass
    else:
        raise AssertionError('Missing or ambiguous initializers must be rejected')
actual = (root/'src/verified_core_taylor_model.repaired.cpp').read_text()
assert Fraction(36, 100) < read(actual, 'b0') < Fraction(37, 100)
assert Fraction(69, 100) < read(actual, 'R0') < Fraction(70, 100)
print('PASS: exact initializer parsing across C++ line wrapping; missing and ambiguous inputs rejected')
