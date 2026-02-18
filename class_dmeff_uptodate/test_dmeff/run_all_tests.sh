#!/bin/bash
#
# run_all_tests.sh — Run all dmeff validation tests and report results.
#
# Usage (from the class_dmeff_uptodate/ directory):
#     bash test_dmeff/run_all_tests.sh
#
# This script:
#   1. Verifies the CLASS executable exists
#   2. Runs each test case (deleting stale outputs first)
#   3. Compares outputs against reference baselines
#   4. Prints a summary table
#
# Exit code: 0 if all tests pass, 1 if any test fails.

set -e

# Determine paths
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CLASS_DIR="$(dirname "$SCRIPT_DIR")"
CLASS_EXE="${CLASS_DIR}/class"
COMPARE="${SCRIPT_DIR}/compare_outputs.py"
OUTPUT_DIR="${CLASS_DIR}/output"

# Colors for terminal output (disable if not a terminal)
if [ -t 1 ]; then
    GREEN='\033[0;32m'
    RED='\033[0;31m'
    YELLOW='\033[0;33m'
    BOLD='\033[1m'
    NC='\033[0m'
else
    GREEN='' RED='' YELLOW='' BOLD='' NC=''
fi

echo "============================================"
echo "  dmeff Validation Test Suite"
echo "============================================"
echo ""

# Check prerequisites
if [ ! -x "$CLASS_EXE" ]; then
    echo -e "${RED}ERROR${NC}: CLASS executable not found at ${CLASS_EXE}"
    echo "Build CLASS first:  make clean && make"
    exit 1
fi

if [ ! -f "$COMPARE" ]; then
    echo -e "${RED}ERROR${NC}: Comparison script not found at ${COMPARE}"
    exit 1
fi

# Ensure output directory exists
mkdir -p "$OUTPUT_DIR"

# Define test cases
TESTS="test_coulomb test_constant test_electron test_mixed test_multi test_vanilla"

PASSED=0
FAILED=0
TOTAL=0

for TEST in $TESTS; do
    TOTAL=$((TOTAL + 1))
    INI="${SCRIPT_DIR}/${TEST}.ini"

    echo ""
    echo -e "${BOLD}[$TOTAL] Running ${TEST}...${NC}"

    if [ ! -f "$INI" ]; then
        echo -e "  ${RED}SKIP${NC}: INI file not found: ${INI}"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Delete stale output files
    rm -f "${OUTPUT_DIR}/${TEST}_"*

    # Run CLASS
    if ! "$CLASS_EXE" "$INI" > /dev/null 2>&1; then
        echo -e "  ${RED}FAIL${NC}: CLASS returned an error"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Verify output files were created
    OUTPUT_COUNT=$(ls "${OUTPUT_DIR}/${TEST}_"* 2>/dev/null | wc -l)
    if [ "$OUTPUT_COUNT" -eq 0 ]; then
        echo -e "  ${RED}FAIL${NC}: No output files generated"
        FAILED=$((FAILED + 1))
        continue
    fi

    # Run comparison
    if python "$COMPARE" "$TEST" --check all --tolerance 0.1 2>&1; then
        PASSED=$((PASSED + 1))
    else
        FAILED=$((FAILED + 1))
    fi
done

# Print summary
echo ""
echo ""
echo "============================================"
echo "  Summary"
echo "============================================"
echo ""

for TEST in $TESTS; do
    INI="${SCRIPT_DIR}/${TEST}.ini"
    if [ ! -f "$INI" ]; then
        echo -e "  ${TEST}: ${YELLOW}SKIP${NC} (no INI file)"
        continue
    fi
    OUTPUT_COUNT=$(ls "${OUTPUT_DIR}/${TEST}_"* 2>/dev/null | wc -l)
    if [ "$OUTPUT_COUNT" -gt 0 ]; then
        echo -e "  ${TEST}: ${GREEN}RAN${NC} ($OUTPUT_COUNT output files)"
    else
        echo -e "  ${TEST}: ${RED}NO OUTPUT${NC}"
    fi
done

echo ""
echo "============================================"
if [ "$FAILED" -eq 0 ]; then
    echo -e "  ${GREEN}${BOLD}OVERALL: ${PASSED}/${TOTAL} tests PASSED${NC}"
else
    echo -e "  ${RED}${BOLD}OVERALL: ${FAILED}/${TOTAL} tests FAILED${NC}"
fi
echo "============================================"
echo ""

# Remind about known baselines
echo -e "${YELLOW}Note:${NC} ~0.20% differences at l=2000 and ~0.30% in Tb(z=10)"
echo "are expected CLASS v2.9.4 -> v3.3.4 baselines (HyRec update),"
echo "not dmeff porting errors. These appear identically in all tests"
echo "including vanilla."

if [ "$FAILED" -gt 0 ]; then
    exit 1
fi
exit 0
