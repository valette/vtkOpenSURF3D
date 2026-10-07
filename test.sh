#!/bin/bash
#
# vtkOpenSURF3D test suite.
#
# Runs surf3d (and match3d when available) through several scenarios and
# validates that the expected output files are produced. All artifacts are
# written into the ./test_results directory to keep the repo tree clean.
#
# Usage:
#   ./test.sh            run the default test battery
#   ./test.sh all        run the test battery against every image in niivue-images/
#
# Exit code is 0 only if every test passed.

set -u

TESTDIR="test_results"
IMG_DIR="niivue-images"
DEFAULT_IMG="$IMG_DIR/CT_Abdo.nii.gz"

PASS=0
FAIL=0
FAILED_TESTS=()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# run_test <name> <expected_suffix...> -- <command...>
# Executes <command...> (writing outputs into $TESTDIR/<name>), then checks
# that every file "$TESTDIR/<name><suffix>" exists and is non-empty.
run_test() {
    local name="$1"; shift
    local expected=()
    while [ "$1" != "--" ]; do
        expected+=("$1"); shift
    done
    shift  # drop the "--"

    local outfile="$TESTDIR/$name"
    local logfile="$TESTDIR/$name.log"
    local ok=1

    # Redirect the run's own output to a per-test log for inspection.
    "$@" >"$logfile" 2>&1
    local ret=$?

    if [ "$ret" -ne 0 ]; then
        echo "[FAIL] $name : command exited with code $ret"
        ok=0
    fi

    for suf in "${expected[@]}"; do
        if [ ! -s "$outfile$suf" ]; then
            echo "[FAIL] $name : missing or empty expected output '$outfile$suf'"
            ok=0
        fi
    done

    if [ "$ok" -eq 1 ]; then
        PASS=$((PASS+1))
        echo "[PASS] $name"
    else
        FAIL=$((FAIL+1))
        FAILED_TESTS+=("$name")
    fi
}

# expect_fail <name> -- <command...>
# Runs <command...> and requires it to return a non-zero exit code.
expect_fail() {
    local name="$1"; shift
    shift  # drop the "--"

    local logfile="$TESTDIR/$name.log"
    "$@" >"$logfile" 2>&1
    local ret=$?

    if [ "$ret" -ne 0 ]; then
        PASS=$((PASS+1))
        echo "[PASS] $name (failed as expected, code $ret)"
    else
        FAIL=$((FAIL+1))
        FAILED_TESTS+=("$name")
        echo "[FAIL] $name : expected failure but command succeeded"
    fi
}

summary() {
    echo
    echo "=================== SUMMARY ==================="
    echo "Passed : $PASS"
    echo "Failed : $FAIL"
    if [ "$FAIL" -ne 0 ]; then
        printf 'Failed tests:\n'
        for t in "${FAILED_TESTS[@]}"; do
            printf '  - %s\n' "$t"
        done
    fi
    echo "Artifacts are in ./$TESTDIR/"
    [ "$FAIL" -eq 0 ]
}

# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

if [ ! -f "$DEFAULT_IMG" ]; then
    echo "ERROR: default image '$DEFAULT_IMG' not found (niivue-images/ missing?)" >&2
    exit 2
fi

rm -rf "$TESTDIR"
mkdir -p "$TESTDIR"

# ---------------------------------------------------------------------------
# Image selection
# ---------------------------------------------------------------------------

images=()
if [ "${1:-}" = "all" ]; then
    while IFS= read -r f; do
        images+=("$f")
    done < <(find "$IMG_DIR" -name '*.nii.gz' | sort)
    echo "Testing ${#images[@]} images from $IMG_DIR"
else
    images=("$DEFAULT_IMG")
fi

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

# --- 1. Per-image baseline: default detection on every selected image ------
for img in "${images[@]}"; do
    base="$(basename "$img" .nii.gz)"
    run_test "baseline_$base" ".csv.gz" ".json" -- ./surf3d "$img" -t 0.001 -o "$TESTDIR/baseline_$base"
done

# --- 2. Output formats ------------------------------------------------------
run_test "fmt_json" ".json" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -json 1 -csv 0 -csvgz 0 -bin 0 -o "$TESTDIR/fmt_json"
run_test "fmt_csv" ".csv" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -json 0 -csv 1 -csvgz 0 -bin 0 -o "$TESTDIR/fmt_csv"
run_test "fmt_bin" ".bin" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -json 0 -csv 0 -csvgz 0 -bin 1 -o "$TESTDIR/fmt_bin"
run_test "fmt_all" ".json" ".csv" ".bin" ".csv.gz" -- \
    ./surf3d "$DEFAULT_IMG" -t 0.001 -json 1 -csv 1 -csvgz 1 -bin 1 -o "$TESTDIR/fmt_all"

# --- 3. Threshold sweep -----------------------------------------------------
run_test "thresh_low" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.0001 -o "$TESTDIR/thresh_low"
run_test "thresh_med" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001  -o "$TESTDIR/thresh_med"
run_test "thresh_high" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.01   -o "$TESTDIR/thresh_high"

# --- 4. Descriptor types ----------------------------------------------------
run_test "desc_surf3d" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -type 0 -o "$TESTDIR/desc_surf3d"
run_test "desc_haar" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -type 1 -r 3 -o "$TESTDIR/desc_haar"
run_test "desc_raw" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -type 2 -r 3 -o "$TESTDIR/desc_raw"

# --- 5. Threading -----------------------------------------------------------
run_test "threads_1" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -nt 1 -o "$TESTDIR/threads_1"
run_test "threads_4" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -nt 4 -o "$TESTDIR/threads_4"

# --- 6. Resampling ----------------------------------------------------------
run_test "resample_spacing" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -s 2.0 -o "$TESTDIR/resample_spacing"
run_test "resample_maxsize" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -d 100 -o "$TESTDIR/resample_maxsize"

# --- 7. Normalize on/off ----------------------------------------------------
run_test "normalize_on" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -normalize 1 -o "$TESTDIR/normalize_on"
run_test "normalize_off" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -normalize 0 -o "$TESTDIR/normalize_off"

# --- 8. Padding -------------------------------------------------------------
run_test "padding_5" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -pad 5 -o "$TESTDIR/padding_5"

# --- 9. csv.gz precision and compression -----------------------------------
run_test "precision_3" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -precision 3 -o "$TESTDIR/precision_3"
run_test "gz_opts" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -gz 6 -o "$TESTDIR/gz_opts"

# --- 10. Pointfile reuse (descriptor-only pass, expects a CSV) --------------
if [ -s "$TESTDIR/fmt_csv.csv" ]; then
    run_test "pointfile" ".csv.gz" -- \
        ./surf3d "$DEFAULT_IMG" -p "$TESTDIR/fmt_csv.csv" \
        -csv 0 -csvgz 1 -bin 0 -o "$TESTDIR/pointfile"
fi

# --- 11. Limit number of points ---------------------------------------------
run_test "maxpoints" ".csv.gz" -- ./surf3d "$DEFAULT_IMG" -t 0.001 -n 10 -o "$TESTDIR/maxpoints"

# --- 12. Error handling -----------------------------------------------------
expect_fail "err_missing_image" -- ./surf3d "$TESTDIR/does_not_exist.nii.gz"
expect_fail "err_missing_value" -- ./surf3d "$DEFAULT_IMG" -t

# --- 13. match3d (only if built and loadable) -------------------------------
# match3d writes its result to a hard-coded "transform.json" in the CWD, so it
# is executed from within $TESTDIR to keep the repo tree clean.
if ./match3d 2>&1 | grep -q "Usage"; then
    run_test "transform" ".json" -- \
        bash -c 'cd "$1" && ../match3d ../fmt_json.json ../fmt_json.json -i' _ "$TESTDIR"
else
    echo "[SKIP] match3d present but not loadable in this environment (missing shared libs?)"
fi

# ---------------------------------------------------------------------------
summary
