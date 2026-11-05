#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

# The library is compiled separately, as requested in section 2.5.
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic -c linal.c -o linal.o
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic main.c linal.o -lm -o program.exe
gcc -std=c11 -O2 -Wall -Wextra -Wpedantic tests.c linal.o -lm -o tests.exe

./tests.exe

# Strip possible CR characters from Windows console output (MinGW / Git Bash).
actual="$(./program.exe | tr -d '\r')"
expected='A (2x3):
 1  2  3
 4  5  6
B (3x2):
 7  8
 9  10
 11  12
C = A*B (2x2):
 58  64
 139  154
|u| = 5
u dot v = 11'

if [[ "$actual" != "$expected" ]]; then
    printf '%s\n' 'FAIL: demonstration output differs from expected:' >&2
    printf '%s\n' "$actual" >&2
    exit 1
fi

printf '%s\n' 'CLI demonstration: PASS'
