#!/usr/bin/env bash
# Console-level regression tests. Run in Git Bash after building program.exe.

set -uo pipefail

project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
program="$project_dir/program.exe"
if [[ ! -f "$program" ]]; then
    echo 'Compile program.exe before running test_cli.sh' >&2
    exit 1
fi

passed=0
failed=0

pass() { printf '[PASS] %s\n' "$1"; passed=$((passed + 1)); }
fail() { printf '[FAIL] %s\n' "$1" >&2; failed=$((failed + 1)); }

expect_success() {
    local name="$1" required="$2" output
    shift 2
    if ! output=$("$program" "$@" 2>&1); then
        fail "$name (unexpected error: $output)"
    elif [[ "$output" == *"$required"* ]]; then
        pass "$name"
    else
        fail "$name (missing text: $required)"
    fi
}

expect_failure() {
    local name="$1" output
    shift
    if output=$("$program" "$@" 2>&1); then
        fail "$name (unexpected success: $output)"
    else
        pass "$name"
    fi
}

expect_success 'quadratic: textbook formula (2.23)' \
    'Formula (2.23):    0.33333333333333331' x2 0 1 100 1e-9
expect_success 'quadratic: trapezoid approximation' \
    'Trapezoidal (2.21): 0.333349' x2 0 1 100 1e-9
expect_success 'tiny quadratic integral: adaptive scale' \
    'Adaptive:          3.33333' x2 0 1e-12 100 1e-6
expect_success 'tiny negative integral: adaptive scale' \
    'Adaptive:          -3.33333' x2 1e-12 0 100 1e-6
expect_success 'sine: adaptive approximate value' \
    'Adaptive:          1.99' sin 0 3.141592653589793 100 1e-8
expect_success 'constant: exact integral' \
    'Formula (2.23):    1' one 0 1 5 1e-9
expect_success 'reversed bounds: negative integral' \
    'Formula (2.23):    -0.333' x2 1 0 5 1e-8
expect_success 'equal bounds: zero integral' \
    'Adaptive:          0' x2 7 7 5 1e-8

expect_failure 'missing arguments' x2 0 1
expect_failure 'unknown function' unknown 0 1 100 1e-8
expect_failure 'n too small for (2.23)' x2 0 1 4 1e-8
expect_failure 'negative n' x2 0 1 -5 1e-8
expect_failure 'text instead of n' x2 0 1 2x 1e-8
expect_failure 'NaN bound' x2 nan 1 100 1e-8
expect_failure 'infinite bound' x2 0 inf 100 1e-8
expect_failure 'zero tolerance' x2 0 1 100 0
expect_failure 'NaN tolerance' x2 0 1 100 nan
expect_failure 'unrepresentably small tolerance' x2 0 1 100 1e-400
expect_failure 'non-finite integrand' exp 0 1000 100 1e-8
expect_failure 'overflowing integral' one -1e308 1e308 100 1e-8

printf 'CLI checks: %d; failures: %d\n' "$((passed + failed))" "$failed"
((failed == 0))
