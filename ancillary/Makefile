# Revised certificate chain. Run `make verify` with GCC, MPFR, GMP,
# libquadmath, Python and SymPy. Optional construction tools use mpmath.
# GNU Make's built-in CXX=g++ can be Apple Clang, so prefer installed MacPorts GCC.
ifeq ($(origin CXX),default)
CXX := $(shell command -v g++-mp-15 2>/dev/null || command -v g++-mp-14 2>/dev/null || command -v g++)
endif
PYTHON ?= python3
MPFR_PREFIX ?= $(shell if test -f /opt/local/include/mpfr.h; then echo /opt/local; elif test -f /opt/homebrew/include/mpfr.h; then echo /opt/homebrew; else echo /usr; fi)
CPPFLAGS += -Iinclude -I$(MPFR_PREFIX)/include
CXXFLAGS += -std=gnu++17 -O2 -fno-fast-math -ffp-contract=off
LDFLAGS += -L$(MPFR_PREFIX)/lib
LDLIBS += -lmpfr -lgmp

CORE = build/core_reconstructed
HORIZON = build/horizon_positivity build/horizon_second_derivative
TESTS = build/test_mpfr_interval build/test_float128_compat build/check_core_repairs

.PHONY: all test certificates verify symbolic paper-data weights check-core-source regenerate-core clean
all: $(CORE) $(HORIZON) $(TESTS)

build logs:
	mkdir -p $@

check-core-source:
	$(PYTHON) patches/patch_core_validation.py --check-output src/verified_core_taylor_model.repaired.cpp

# Maintenance only: ordinary verification checks the distributed source without
# changing it, and therefore does not require a formatter.
regenerate-core:
	$(PYTHON) patches/patch_core_validation.py --output-source src/verified_core_taylor_model.repaired.cpp

$(CORE): src/verified_core_taylor_model.repaired.cpp src/mpfr_float128_compat.cpp include/mpfr_interval.hpp include/fixed_weight_schedule.hpp | build check-core-source
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) src/verified_core_taylor_model.repaired.cpp src/mpfr_float128_compat.cpp $(LDFLAGS) $(LDLIBS) -lquadmath -o $@

build/horizon_%: src/horizon_%.cpp include/mpfr_interval.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

build/test_mpfr_interval: tests/test_mpfr_interval.cpp include/mpfr_interval.hpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $< $(LDFLAGS) $(LDLIBS) -o $@

build/test_float128_compat: src/test_float128_compat.cpp src/mpfr_float128_compat.cpp | build
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $^ $(LDFLAGS) $(LDLIBS) -lquadmath -o $@

build/check_core_repairs: patches/check_core_repairs.cpp src/verified_core_taylor_model.repaired.cpp src/mpfr_float128_compat.cpp include/mpfr_interval.hpp include/fixed_weight_schedule.hpp | build check-core-source
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) patches/check_core_repairs.cpp src/mpfr_float128_compat.cpp $(LDFLAGS) $(LDLIBS) -lquadmath -o $@

test: $(TESTS) | logs
	build/test_mpfr_interval > logs/test_mpfr_interval.txt
	build/test_float128_compat > logs/test_float128_compat.txt
	build/check_core_repairs > logs/check_core_repairs.txt
	$(PYTHON) tests/test_paper_data_parser.py > logs/test_paper_data_parser.txt

certificates: $(CORE) $(HORIZON) | logs
	build/horizon_positivity > logs/horizon_positivity.txt
	build/horizon_second_derivative > logs/horizon_second_derivative.txt
	build/core_reconstructed 40 > logs/core_reconstructed.txt 2> logs/core_reconstructed.log

verify: test certificates
	$(PYTHON) python/check_reconstruction.py > logs/verification.txt
	$(PYTHON) python/symbolic/verify_constraints.py > logs/symbolic_constraints.txt
	$(PYTHON) python/matching/generate_paper_data.py > logs/paper_data_generation.txt
	@echo "PASS: arithmetic, rebuilt horizon/core, horizon derivations and initialization, tail, exact matching, symbolic identities and outward paper data. See logs/."

symbolic: | logs
	$(PYTHON) python/symbolic/verify_constraints.py > logs/symbolic_constraints.txt

paper-data: | logs
	$(PYTHON) python/matching/generate_paper_data.py > logs/paper_data_generation.txt

weights:
	$(PYTHON) weights/generate_fixed_weights.py

clean:
	rm -rf build
