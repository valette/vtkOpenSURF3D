# Codebase Analysis, Bugs & Improvements (TODO)

This document provides a comprehensive technical audit of the **vtkOpenSURF3D** codebase. It outlines identified bugs, memory leaks, undefined behaviors, mathematical inconsistencies, architecture issues, and recommended improvements, prioritized by severity.

> **Status legend:** `[STILL OPEN]` = issue remains in the code · `[PARTIAL]` =
> only partially addressed. Resolved issues have been removed from this audit.

---

## Table of Contents
1. [Memory Management & Resource Leaks (High Priority)](#1-memory-management--resource-leaks-high-priority)
2. [Algorithmic & Mathematical Inconsistencies (Medium Priority)](#2-algorithmic--mathematical-inconsistencies-medium-priority)
3. [Safety, Robustness & Error Handling (Medium Priority)](#3-safety-robustness--error-handling-medium-priority)
4. [Code Quality, Modern C++ & Refactoring (Medium Priority)](#4-code-quality-modern-c--refactoring-medium-priority)
5. [Build System, CI & Repository Hygiene (Low Priority)](#5-build-system-ci--repository-hygiene-low-priority)
6. [Documentation & CLI Consistency (Low Priority)](#6-documentation--cli-consistency-low-priority)
7. [Remediation Plan (Proposed Work Order)](#7-remediation-plan-proposed-work-order)

---

## 1. Memory Management & Resource Leaks (High Priority)

- [ ] **[M1][STILL OPEN (minor)] Per-candidate heap allocations in `interpolateStep` (`fasthessian.cxx:631-657`)**
  - **Issue**: `deriv4D()` and `hessian4D()` dynamically allocate `new cv::Matx41d()` and `new cv::Matx44d()`.
  - **Impact**: Unnecessary heap allocation and deallocation overhead for thousands of candidate points.
  - **Fix**: Return `cv::Matx41d` and `cv::Matx44d` by value on the stack.

---

## 2. Algorithmic & Mathematical Inconsistencies (Medium Priority)

- [ ] **[A1][STILL OPEN] Rigid registration in `MatchPoint::getRT` hardcodes `R = Identity` (`MatchPoint.cxx:323-344`)**
  - **Issue**: The Umeyama SVD rotation calculation is commented out. `Transform.get_rotation() = R;` assigns identity.
  - **Impact**: `match3d` cannot estimate 3D rotations, only isotropic scale and translation.
  - **Fix**: Uncomment and validate the SVD-based rotation estimator, ensuring proper reflection handling when $\det(U V^T) < 0$.

- [ ] **[A2][STILL OPEN] Scale divisor check in `MatchPoint::getRT` (`MatchPoint.cxx:321`)**
  - **Issue**: `double s = sqrt((EigCa * EigCb) / (EigCa * EigCa));`
  - **Impact**: If points in set A are collinear or identical, `EigCa * EigCa` can be 0, causing division by zero.
  - **Fix**: Guard against near-zero denominator before computing scale.

- [ ] **[A3][STILL OPEN] Inconsistent coordinate conventions across function signatures**
  - **Issue**: Some methods use `(column, row, layer)` = $(x, y, z)$, while others pass `(r, c, d)` = $(y, x, z)$ or `(d, r, c)` = $(z, y, x)$ (e.g., `interpolateExtremum(d, r, c, ...)` vs `isExtremum(r, c, d, ...)`).
  - **Impact**: Highly error-prone for future maintainers.
  - **Fix**: Adopt a consistent parameter naming and ordering convention (e.g. `x, y, z` or `row, col, slice`).

---

## 3. Safety, Robustness & Error Handling (Medium Priority)

- [ ] **[S1][STILL OPEN] Operator `-` semantics in `Ipoint` (`ipoint.h:35-45`)**
  - **Issue**: `float operator-(const Ipoint &rhs)` computes the Euclidean distance between descriptor vectors.
  - **Impact**: Overloading `operator-` to return a scalar distance is counter-intuitive in C++ (where subtraction usually returns a difference vector or offset). Furthermore, it does not verify that `this->descriptor.size() == rhs.descriptor.size()`.
  - **Fix**: Replace with a named method `float distanceTo(const Ipoint& rhs) const` and assert matching descriptor lengths.

---

## 4. Code Quality, Modern C++ & Refactoring (Medium Priority)

- [ ] **[C1][STILL OPEN] Replace global macro definitions with inline functions**:
  - `integ(x,y,z)` and `sourc(x,y,z)` in `integral.h:6-7` pollute the global preprocessor namespace and rely on implicitly named local pointers (`pin_tegral`, `int_egral_incs`).
  - `PRINT`, `CHECKMAT`, `PRINTIF` in `fasthessian.cxx`.
  - **Fix**: Replace macros with typed inline helper functions or lambda accessors.

- [ ] **[C2][STILL OPEN] Remove dead / commented-out code**:
  - `Surf::ThreadedDesc` declared in `surf.h:41` but never implemented.
  - `FastHessian::EigenValue` defined in `fasthessian.cxx:487-520` but never used.
  - `FastHessian::FittingQuadric` (~200 lines) dead code in `fasthessian.cxx`.
  - Unused variables generating compiler warnings (`imin2`, `temp2`, `test`, `bounds`, `img_pointer`, `layer_min`, `layer_max`).

- [ ] **[C3][STILL OPEN] Eliminate compiler warnings under `-Wall -Wextra`**:
  - Constructor member initialization order warnings (`-Wreorder`) in `Ipoint`, `FastHessian`.
  - Signed vs. unsigned comparisons (`-Wsign-compare`) across loops in `vtk3DSURF.cxx`, `surf.cxx`, `fasthessian.cxx`.
  - Unused parameters (`-Wunused-parameter`) in `BoxIntegralOptim` and `buildResponseLayer`.

---

## 5. Build System, CI & Repository Hygiene (Low Priority)

- [ ] **[B1][PARTIAL] Modernize `CMakeLists.txt`**:
  - Move `cmake_minimum_required(VERSION 3.20)` to line 1 outside any `if()` blocks (currently requires 3.30 inside an `if` statement).
  - Explicitly set `set(CMAKE_CXX_STANDARD 17)` and `set(CMAKE_CXX_STANDARD_REQUIRED ON)`.
  - Modernize target linking: use `target_link_libraries(surf3d PRIVATE VTK::... OpenMP::OpenMP_CXX ZLIB::ZLIB ${OpenCV_LIBS})` instead of global `include_directories` and modifying `CMAKE_CXX_FLAGS`.
  - Add `find_package(LAPACK)` / `find_package(BLAS)` and TooN include path resolution when `BUILD_SURFMATCH=ON`.
  - Add `install(TARGETS surf3d RUNTIME DESTINATION bin)`.
  - **Partial**: `CMAKE_CXX_STANDARD 17` / `REQUIRED ON` set and `target_link_libraries(... PRIVATE ...)` used for `surf3d`. Still open: `cmake_minimum_required` remains inside the `if()` block (line 2), global `include_directories` still used, no `find_package(LAPACK/BLAS)`/TooN resolution, no `install()`.

- [ ] **[B2][PARTIAL] Update GitHub Actions CI (`.github/workflows/ci.yml`)**:
  - Upgrade `actions/checkout@v2` to `@v4` to prevent deprecation warnings.
  - Add a testing step with sample image data to verify feature detection outputs (`ctest`).
  - **Partial**: A `./test.sh -v` testing step was added (ci.yml:63-67). Still open: `actions/checkout@v2` not yet upgraded to `@v4`.

---

## 6. Documentation & CLI Consistency (Low Priority)

- [ ] **[D1][STILL OPEN] Document missing command-line options in `Readme.md` and usage printout**:
  - Supported CLI options implemented in `surf3d.cxx` but omitted from `Readme.md` or `surf3d` usage text:
    - `-d <size>`: set maximum dimension size
    - `-m <mask>`: path to binary mask image
    - `-nt <threads>`: number of OpenMP/VTK execution threads
    - `-p <pointfile>`: pre-computed keypoints file for descriptor computation only
    - `-pad <val>`: apply mirror padding before feature extraction
    - `-gz <opts>`: gzip compression options
    - `-precision <n>`: coefficient precision for floating-point CSV.GZ export
  - **Fix**: Align `surf3d.cxx` usage output and `Readme.md` option tables.

---

## 7. Remediation Plan (Proposed Work Order)

Ordered by priority. Each item lists the current `file:line` and the proposed fix.

### 7.1 Matching correctness
- [ ] **[A1] Re-enable rotation in `getRT`** (`MatchPoint.cxx:323-344`) — uncomment and
  validate the SVD/Umeyama rotation with proper reflection handling
  (`det(U V^T) < 0`).

### 7.2 Safety / robustness
- [ ] **[A2] Guard `getRT` scale divisor** (`MatchPoint.cxx:321`) — skip/bail when
  `EigCa*EigCa` is near zero.
- [ ] **[S1] Rename `Ipoint::operator-`** (`ipoint.h:35-45`) — replace with a named
  `distanceTo(const Ipoint&) const` and assert matching descriptor lengths.

### 7.3 Cleanup / minor
- [ ] **[C2] Remove dead code** (`surf.h:41` `ThreadedDesc`; `fasthessian.cxx:479`
  `EigenValue`; `fasthessian.cxx:716` `FittingQuadric`; unused variables).
- [ ] **[M1] Return matrices by value** (`fasthessian.cxx:623-624,665,hessian4D`) —
  have `deriv4D`/`hessian4D` return `cv::Matx41d`/`cv::Matx44d` on the stack.
- [ ] **[C1] Replace macros** (`integral.h:6-7` `integ/sourc`; `fasthessian.cxx`
  `PRINT`/`CHECKMAT`/`PRINTIF`) with inline functions.
- [ ] **[C3] Triage `-Wall -Wextra` warnings** — `-Wreorder`, `-Wsign-compare`,
  `-Wunused-parameter` across the sources.

### 7.4 Build / CI / docs (low priority)
- [ ] **[B1] Modernize `CMakeLists.txt`** — move `cmake_minimum_required` to line 1,
  replace global `include_directories`, add `find_package(LAPACK/BLAS)` + TooN
  path for `BUILD_SURFMATCH`, add `install(TARGETS surf3d ...)`.
- [ ] **[B2] Upgrade CI** (`ci.yml:28`) — `actions/checkout@v2` → `@v4`.
- [ ] **[D1] Document missing CLI options** — add `-d`, `-m`, `-nt`, `-p`, `-pad`,
  `-gz`, `-precision` to `Readme.md` option table.

