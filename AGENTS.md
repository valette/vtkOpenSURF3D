# AGENTS.md — Guidance for coding agents

> **Last updated:** 2026-10-07

This file helps coding agents (LLMs / AI assistants) work effectively in this
repository. Read it before making changes.

---

## Project at a glance

vtkOpenSURF3D is a C++17, VTK-based implementation of **3D SURF** (Speeded Up
Robust Features) for volumetric images. It detects blob keypoints, computes
local descriptors, and can match two keypoint sets to estimate a similarity
transform.

Two executables:
- **`surf3d`** — detect + describe keypoints, write output (built by default).
- **`match3d`** — match two JSON keypoint sets, RANSAC → `transform.json`
  (optional build, `BUILD_SURFMATCH=ON`).

There is no test framework beyond `test.sh`. There is **no linter, formatter,
or typecheck command** configured; the C++ is checked by the compiler and by
running the test suite.

---

## Essential commands

```bash
# Configure + build surf3d
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel

# Build with debug symbols and warnings (recommended for development)
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_FLAGS="-Wall -Wextra"
cmake --build build --parallel

# Build match3d too
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_SURFMATCH=ON
cmake --build build --parallel

# Fetch test data (git submodule)
git submodule update --init --recursive

# Run the test suite (must pass before/after changes)
./test.sh          # 24 tests, default image
./test.sh -a       # 40 tests, every image in niivue-images/
./test.sh -v       # verbose mode: stream program output to the terminal
```

> There is **no `make test` / `ctest`** — the check command is `./test.sh`.

---

## Repository layout

| Path | Purpose |
|------|---------|
| `surf3d.cxx` | `surf3d` CLI: args, preprocessing, output |
| `vtk3DSURF.h/.cxx` | pipeline orchestration + all writers |
| `fasthessian.h/.cxx` | scale-space detection (`FastHessian`) |
| `responselayer.h` | response layers |
| `integral.h/.cxx` | integral image + `BoxIntegral` |
| `surf.h/.cxx` | Haar descriptors (`Surf`) |
| `ipoint.h` | `Ipoint` struct + distance |
| `vtkRobustImageReader.h` | image reader (flips, NIFTI rescale) |
| `MainMatch.cxx`, `MatchPoint.h/.cxx` | `match3d`: matching + RANSAC |
| `picojson.h` | vendored single-header JSON |
| `CMakeLists.txt` | build |
| `test.sh` | test suite |
| `docs/` | full technical documentation |
| `TODO.md` | audit of known bugs / improvements |

Read `docs/README.md` first for a map, then `docs/architecture.md`,
`docs/modules.md`, and `docs/algorithm.md`.

---

## Code conventions

- **Language / standard:** C++17. Style is mixed/legacy; new code should use
  modern C++ (`std::vector`, smart pointers, `auto`, range-for) where it does
  not clash with surrounding code.
- **VTK objects** use `vtkNew`/`vtkSmartPointer` and `vtkTypeMacro`.
- **Do not add comments** unless they genuinely clarify non-obvious logic;
  match the surrounding style.
- **Naming:** camelCase members/functions; file-scope constants are
  `UPPER_CASE` in headers (e.g. `OCTAVES`, `INTERVALS`, `THRES`,
  `INIT_SAMPLE` in `fasthessian.h`).
- **Do not reformat unrelated code** — keep diffs minimal.

---

## Critical gotchas

1. **Keypoint coordinates/scale are in VOXELS** of the internal `Cast` image
   until export. Exporters convert to physical space. If you add a writer or
   reader, you must apply the same origin/spacing conversion.
2. **`BoxIntegralOptim` is unchecked** (no bounds clamping). Callers must
   guarantee in-bounds boxes. The clamped `BoxIntegral` is the safe default.
3. **`match3d` hard-codes its output to `transform.json` in the CWD.**
   Tests run it from a subshell to keep the tree clean.
4. **`match3d` rotation is currently `Identity`** — the SVD/Umeyama rotation
   is commented out. Don't rely on it producing a real rotation.
5. **`-nt` only controls VTK filters**, not OpenMP. OpenMP thread counts are
   set by `OMP_NUM_THREADS`.
6. **The root `match3d` binary is stale** and may not load against current
   VTK. Rebuild with `BUILD_SURFMATCH=ON` before relying on it.
7. **`test_results/` is git-ignored** — never commit test artifacts.
8. **`niivue-images/` is a git submodule** — `git clone` alone won't populate
   it; use `git submodule update --init`.

---

## Known bugs & open items

A comprehensive, severity-ordered audit lives in **`TODO.md`**. Highlights an
agent should be aware of (details + file:line in TODO.md):

- NIFTI intensity rescale is skipped when no axis flip is needed.
- JSON keypoint `scale` is written as 0 when `-s` is not given.
- `BoxIntegralOptim`/Haar paths can read out of bounds near borders.
- Response-layer buffers may contain uninitialized border voxels.
- The `match3d` rotation is fixed to identity (see gotcha 4).
- `responselayer.h` lacks an include guard.
- The CLI parsers can read past `argc` (segfault on a trailing flag).

Before touching any of these, read the matching section of `TODO.md`.

---

## Workflow for changes

1. Read the relevant module in `docs/modules.md` first.
2. Make the change, keeping conventions above.
3. Build with `-Wall -Wextra` and fix new warnings in the touched code.
4. Run `./test.sh` (and `./test.sh -a` if feasible) and confirm all pass.
5. Do not commit unless asked.

---

## Verification checklist

- [ ] `cmake --build build --parallel` succeeds.
- [ ] No new compiler warnings in modified files.
- [ ] `./test.sh` → all PASS.
- [ ] If you changed matching: `./test.sh -a` also passes.
- [ ] If you added/removed a CLI option: update `docs/usage.md` and `test.sh`.
