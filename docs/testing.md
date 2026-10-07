# Testing

> **Last updated:** 2026-10-07

vtkOpenSURF3D is validated by a shell-based test suite, `test.sh`. There is no
CTest integration; the suite is invoked directly from the shell.

---

## Prerequisites
- A built `surf3d` binary at the repository root.
- The `niivue-images` test-data directory containing at least
  `niivue-images/CT_Abdo.nii.gz`. This is provided as a **git submodule** and
  is needed only for testing:

```bash
git submodule update --init --recursive
```

If the directory is missing, `test.sh` prints an error and exits (code 2).

---

## Running

```bash
# Default battery (uses the single default image)
./test.sh

# Test against every image in niivue-images/
./test.sh -a

# Verbose mode: stream each program's console output to the terminal
./test.sh -v

# Options are combinable and order-independent
./test.sh -a -v
```

Exit code is `0` only if every test passed.

### Options

| Option | Meaning |
|--------|---------|
| `-a` | run against every image in `niivue-images/` |
| `-v` | verbose mode: stream program console output to the terminal (also useful in CI) |

Unknown options print the usage line and exit with a non-zero code.

### Modes

| Mode | Coverage |
|------|----------|
| `./test.sh` | 24 tests against `CT_Abdo.nii.gz` |
| `./test.sh -a` | baseline detection on every image in `niivue-images/` + the 24 standard tests |

---

## Test categories

The suite (`test.sh`) exercises:

1. **Baseline** — default detection on each selected image
   (`-t 0.001`, expects `.csv.gz` + `.json`).
2. **Output formats** — JSON, CSV, BIN, and all-at-once.
3. **Threshold sweep** — `-t 0.0001 / 0.001 / 0.01`.
4. **Descriptor types** — `-type 0 / 1 / 2`.
5. **Threading** — `-nt 1` vs `-nt 4`.
6. **Resampling** — `-s 2.0` and `-d 100`.
7. **Normalization** — `-normalize 1 / 0`.
8. **Padding** — `-pad 5`.
9. **Precision / compression** — `-precision 3`, `-gz 6`.
10. **Pointfile reuse** — `-p` descriptor-only pass.
11. **Point limiting** — `-n 10`.
12. **Error handling** — missing image, missing option value (must fail).
13. **match3d** — runs `match3d` if it is present and loadable.

---

## How a test is validated

Two helpers implement the checks:

- `run_test <name> <suffixes...> -- <command...>` — runs the command, then
  asserts every file `test_results/<name><suffix>` exists and is non-empty.
- `expect_fail <name> -- <command...>` — asserts the command returns
  non-zero (used for error cases).

All command output is captured to `test_results/<name>.log` for inspection,
and produced artifacts stay inside `test_results/` (which is git-ignored),
keeping the repository tree clean.

In **verbose mode** (`-v`), the command output is additionally streamed to
the terminal via `tee`, so it appears live on the console while still being
written to the per-test log. This is convenient for debugging and for CI
logs.

### match3d special handling

`match3d` writes its output to a **hard-coded** `transform.json` in the
current working directory. The test therefore runs it inside a subshell
`cd`'d into `test_results/`. If the `match3d` binary exists but cannot load
(missing shared libraries), the test is skipped with a notice rather than
failed.

---

## Artifacts

Everything is written under `test_results/`:

```
test_results/
  baseline_CT_Abdo.csv.gz, .json      (per-image baseline)
  fmt_json.json
  fmt_csv.csv
  fmt_bin.bin
  ...
  transform.json                       (from match3d, if run)
  *.log                                (captured command output)
```

`test_results/` is listed in `.gitignore`.
