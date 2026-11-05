#!/usr/bin/env bash
# End-to-end CLI checks for Task 2.3. Run with Git Bash after compiling program.exe.

set -uo pipefail

project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
program="$project_dir/program.exe"
if [[ ! -f "$program" ]]; then
    echo 'Compile program.exe before running test_cli.sh' >&2
    exit 1
fi

test_dir=$(mktemp -d) || exit 1
cleanup() {
    rm -f -- "$test_dir/three.txt" "$test_dir/empty.txt" \
       "$test_dir/one.txt" "$test_dir/constant.txt" \
       "$test_dir/bad.txt" "$test_dir/nonfinite.txt" \
       "$test_dir/numbers with spaces.txt" "$test_dir/extreme.txt"
    rmdir -- "$test_dir"
}
trap cleanup EXIT

printf '1\n2\n3\n' > "$test_dir/three.txt"
printf '' > "$test_dir/empty.txt"
printf '42\n' > "$test_dir/one.txt"
printf '5\n5\n5\n' > "$test_dir/constant.txt"
printf '1\nwrong\n3\n' > "$test_dir/bad.txt"
printf '1\nnan\n3\n' > "$test_dir/nonfinite.txt"
printf '1\n2\n3\n' > "$test_dir/numbers with spaces.txt"
printf '1e200\n-1e200\n' > "$test_dir/extreme.txt"

passed=0
failed=0

pass() { printf '[PASS] %s\n' "$1"; passed=$((passed + 1)); }
fail() { printf '[FAIL] %s\n' "$1" >&2; failed=$((failed + 1)); }

expect_success() {
    local name="$1" required="$2" output
    shift 2
    if ! output=$("$program" "$@" 2>&1); then
        fail "$name (unexpected failure: $output)"
    elif [[ "$output" == *"$required"* ]]; then
        pass "$name"
    else
        fail "$name (missing expected text: $required)"
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

expect_success 'three values: variance with N-1' 'Variance (N-1): 1' "$test_dir/three.txt"
expect_success 'three values: excess kurtosis' 'Excess kurtosis: -2.333' "$test_dir/three.txt"
expect_success 'filename with spaces' 'Count: 3' "$test_dir/numbers with spaces.txt"
expect_success 'single value: undefined moments' 'Variance (N-1): undefined' "$test_dir/one.txt"
expect_success 'constant data: undefined skew' 'Skewness: undefined' "$test_dir/constant.txt"
expect_failure 'empty file rejected' "$test_dir/empty.txt"
expect_failure 'invalid token rejected' "$test_dir/bad.txt"
expect_failure 'NaN rejected' "$test_dir/nonfinite.txt"
expect_failure 'variance outside double range rejected' "$test_dir/extreme.txt"
expect_failure 'missing file rejected' "$test_dir/not-found.txt"
expect_failure 'missing filename rejected'
expect_failure 'too many command-line arguments rejected' "$test_dir/three.txt" extra

if output=$(printf '1\n2\n3\n' | "$program" - 2>&1) \
   && [[ "$output" == *'Count: 3'* ]] \
   && [[ "$output" == *'Variance (N-1): 1'* ]]; then
    pass 'stdin pipe: processing without rewind'
else
    fail 'stdin pipe: processing without rewind'
fi

printf 'CLI checks: %d; failures: %d\n' "$((passed + failed))" "$failed"
((failed == 0))
