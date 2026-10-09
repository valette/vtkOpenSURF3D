# Codebase Analysis, Bugs & Improvements (TODO)

This document provides a comprehensive technical audit of the **vtkOpenSURF3D** codebase. It outlines identified bugs, memory leaks, undefined behaviors, mathematical inconsistencies, architecture issues, and recommended improvements, prioritized by severity.

> **Status legend:** `[STILL OPEN]` = issue remains in the code · `[FIXED]` = issue
> resolved (see note) · `[PARTIAL]` = only partially addressed.
> Checked boxes `[x]` are resolved items kept for the audit trail.

---

## Table of Contents
1. [Critical Bugs & Undefined Behavior (High Priority)](#1-critical-bugs--undefined-behavior-high-priority)
2. [Memory Management & Resource Leaks (High Priority)](#2-memory-management--resource-leaks-high-priority)
3. [Algorithmic & Mathematical Inconsistencies (Medium Priority)](#3-algorithmic--mathematical-inconsistencies-medium-priority)
4. [Safety, Robustness & Error Handling (Medium Priority)](#4-safety-robustness--error-handling-medium-priority)
5. [Code Quality, Modern C++ & Refactoring (Medium Priority)](#5-code-quality-modern-c--refactoring-medium-priority)
6. [Build System, CI & Repository Hygiene (Low Priority)](#6-build-system-ci--repository-hygiene-low-priority)
7. [Documentation & CLI Consistency (Low Priority)](#7-documentation--cli-consistency-low-priority)
8. [Remediation Plan (Proposed Work Order)](#8-remediation-plan-proposed-work-order)

---

## 1. Critical Bugs & Undefined Behavior (High Priority)

- [ ] **[STILL OPEN] NIFTI intensity scaling bypassed when no flip is needed (`vtkRobustImageReader.h:118-120`)**
  - **Issue**: If an image requires no axis flips (`!flip[0] && !flip[1] && !flip[2]`), the reader executes `Reader->Delete(); return;` early.
  - **Impact**: Lines 138–156 (which check and apply `GetRescaleSlope()` and `GetRescaleIntercept()`) are completely bypassed. Most standard NIFTI medical images have slopes/intercepts but positive axes, causing incorrect voxel intensities to be processed.
  - **Fix**: Move the rescale slope/intercept check before or outside the early-return condition, or execute flipping and rescaling independently.

- [x] **[FIXED] Scale is always zero in JSON export (`vtk3DSURF.cxx:368`)**
  - **Issue**: In `WritePoints()`, keypoint scale is written as `point.scale * Spacing`.
  - **Impact**: `Spacing` is the member variable storing the `-s` command line option (default: `0.0`). In contrast, `WritePointsCSV`, `WritePointsCSVGZ`, and `WritePointsBinary` all calculate `Sspacing = pow(spacing[0]*spacing[1]*spacing[2], 1.0/3.0)` and output `point.scale * Sspacing`. Consequently, when `-s` is not supplied, all scales in exported JSON files evaluate to `0.0`.
  - **Fix**: Compute geometric mean `Sspacing` from image spacing in `WritePoints()`, matching the CSV and binary exporters.
  - **Resolved**: `WritePoints()` now computes `Sspacing = pow(spacing[0]*spacing[1]*spacing[2], 1.0/3)` and writes `point.scale * Sspacing` (`vtk3DSURF.cxx:359,368`).

- [x] **[FIXED] Double-sized descriptor vector with leading zeros in `MatchPoint::Parse` (`MatchPoint.cxx:53-54`)**
  - **Issue**: `point.allocate(t.size());` calls `std::vector::resize(size)`, initializing `t.size()` elements to `0.0f`. Then the loop calls `point.descriptor.push_back(m.get<double>())`.
  - **Impact**: Descriptors end up with $2 \times$ the intended dimension (e.g., 96 floats instead of 48), where the first 48 entries are `0.0f`. Euclidean distance computations (`operator-`) compare 48 zeros plus shifted values, distorting feature matching.
  - **Fix**: Either use `point.descriptor.reserve(t.size())` with `push_back()`, or keep `allocate()` and assign by index: `point.descriptor[k] = m.get<double>()`.
  - **Resolved**: `Parse()` now uses `descriptor.reserve(t.size())` + `push_back()` (`MatchPoint.cxx:54-55`).

- [ ] **[STILL OPEN] Scale-space limit integer division bug (`fasthessian.cxx:206-209`)**
  - **Issue**: The code tests lower scale boundaries with `r * (int)(t->width / b->width) < LimDownScale`.
  - **Impact**: `t` is the top (coarser) layer and `b` is the bottom (finer) layer, so `t->width < b->width`. Because integer division `(int)(smaller / larger)` truncates to `0`, `0 < LimDownScale` evaluates to `true` whenever octaves change. `param &= FIRST_SCALE` is therefore erroneously executed, disabling downward scale extrema tests.
  - **Fix**: Invert to `(int)(b->width / t->width)` or compute scale factors using floating-point values / layer steps.

- [x] **[FIXED] Out-of-bounds array access in `FastHessian::getIpoints()` (`fasthessian.cxx:163`)**
  - **Issue**: `int layer_max = filter_map[octaves][2];` is executed.
  - **Impact**: `filter_map` has dimensions `[OCTAVES][INTERVALS]` (where `OCTAVES = 5`). If `octaves == 5`, accessing `filter_map[5][2]` is an out-of-bounds read and undefined behavior.
  - **Fix**: Remove `layer_max` (and unused `layer_min`) or bound index to `octaves - 1`.
  - **Resolved**: `layer_max`/`layer_min` removed; `saveParameters()` clamps `octaves`/`intervals`/`init_sample` (`fasthessian.cxx:99-110`), so `filter_map[o][...]` stays in bounds.

- [x] **[FIXED] Command-line argument parser crashes / segfaults (`surf3d.cxx:77-82`, `MainMatch.cxx:27-39`)**
  - **Issue**: Both CLI tools parse arguments with `char *value = argv[argumentsIndex + 1]` without verifying `argumentsIndex + 1 < argc`.
  - **Impact**: Passing a flag without a value or having a trailing option (e.g. `./surf3d file.mhd -bin` or `./match3d f1 f2 -i`) reads `argv[argc]` (`NULL`), leading to immediate segmentation fault. In `MainMatch.cxx`, flag `-i` decrements `argumentsIndex` by 1 after accessing `argv[argumentsIndex + 1]`.
  - **Fix**: Validate `argumentsIndex + 1 < argc` before reading values; handle boolean flags without consuming a following argument.
  - **Resolved**: Both parsers now check `argumentsIndex + 1 >= argc` before reading a value, and `-i` in `MainMatch.cxx` is handled as a boolean flag (`surf3d.cxx:77-82`, `MainMatch.cxx:30-39`).

- [ ] **[STILL OPEN] Missing include guard in `responselayer.h` (`responselayer.h:1-22`)**
  - **Issue**: `responselayer.h` has no `#ifndef RESPONSELAYER_H` header guard or `#pragma once`.
  - **Impact**: Including this header in multiple files or transitively causes redefinition errors.
  - **Fix**: Add standard `#ifndef RESPONSELAYER_H` / `#define RESPONSELAYER_H` guards.

- [x] **[FIXED] Uninitialized class members & heap memory**
  - **`ipoint.h:33`**: `Ipoint` default constructor initializes `response(0), laplacian(0), scale(0)`, but leaves `x, y, z` uninitialized.
  - **`responselayer.h:44-52`**: `cornerResponses`, `laplacian`, and `isblob` buffers are allocated with `new[]` without zero-initialization (`memset` is commented out). Border voxels read during extrema checks contain uninitialized values, yielding non-deterministic keypoint detection.
  - **`MatchPoint.h:29-35`**: `computeBoundingBoxes`, `maxinlier`, `nbPointInA`, and `nbPointInB` are uninitialized primitives.
  - **Resolved**: All member variables are now initialized in constructors (`ipoint.h:32`; `MatchPoint.h:27-37`); `ResponseLayer` buffers use value-initialized allocation `new float[total]()` etc. (`responselayer.h:48-52`).

- [x] **[FIXED] Mask filtering bug & potential NULL pointer dereference (`fasthessian.cxx:245-265`)**
  - **Issue**: Coordinates are calculated as `it->x + i * 0 * it->scale` (`* 0 *` zeroes the offset).
  - **Impact**: The $3 \times 3 \times 3$ loop redundantly checks the exact same voxel 27 times. If coordinates are out of bounds, `Mask->GetScalarPointer(...)` returns `nullptr`, which is dereferenced directly (`*static_cast<unsigned char*>(...)`), crashing with SIGSEGV.
  - **Fix**: Remove the `* 0 *` typo, check image bounds before calling `GetScalarPointer`, and handle scalar types dynamically.
  - **Resolved**: `* 0 *` removed (`fasthessian.cxx:252-254`) and image bounds are checked before dereferencing `GetScalarPointer` (`fasthessian.cxx:256-259`).

- [ ] **[STILL OPEN] Unsigned integer underflow in 3D Box Integrals (`integral.h:63, 109-116`)**
  - **Issue**: `b222 - b221 - b212 - b122 + b112 + b121 + b211 - b111` is evaluated using `unsigned long long`.
  - **Impact**: Intermediate or final negative quantities wrap around modulo $2^{64}$ to $\approx 1.8 \times 10^{19}$. In `BoxIntegral`, `std::max((unsigned long long)0, ...)` fails to clamp because wrapped values are already $> 0$.
  - **Fix**: Cast operands to signed `long long` before subtraction, or accumulate positive and negative terms separately before subtracting.

---

## 2. Memory Management & Resource Leaks (High Priority)

- [x] **[FIXED] Leaks in `vtkRobustImageReader.h`**:
  - `vtkMatrix4x4 *formMatrix = vtkMatrix4x4::New()` (line 95) is never deleted.
  - `vtkImageFlip *flip = vtkImageFlip::New()` (line 128) is allocated inside the axis flip loop and never deleted.
  - `vtkImageShiftScale *shiftScale = vtkImageShiftScale::New()` (line 146) is never deleted.
  - `FileName` uses raw `vtkSetMacro(FileName, char*)` instead of `vtkSetStringMacro`, leaking strings and failing to free memory in the destructor.
  - **Fix**: Use `vtkNew<T>` or `vtkSmartPointer<T>` for all VTK pipeline objects and `vtkSetStringMacro`/`vtkGetStringMacro` for string attributes.
  - **Resolved**: VTK objects now use `vtkNew`/`vtkSmartPointer`; `FileName` uses `vtkSetStringMacro`/`vtkGetStringMacro` (`vtkRobustImageReader.h`).

- [x] **[FIXED] Catastrophic memory leak in `EXTRACT_CORNER` (`fasthessian.cxx:468-469`)**
  - **Issue**: `CvMat* evec = cvCreateMat(3,3,CV_32FC1);` and `CvMat* eval = cvCreateMat(3,1,CV_32FC1);` are allocated on every voxel iteration inside a 3-level nested loop running over millions of voxels.
  - **Impact**: Leaks gigabytes of RAM in seconds, causing OOM termination.
  - **Fix**: Allocate temporary matrices once outside the voxel loops or use OpenCV stack structures (`cv::Matx33f`, `cv::Vec3f`, `cv::eigen`).
  - **Resolved**: Uses stack `cv::Matx33f` + `cv::Vec3f`/`cv::eigen` (`fasthessian.cxx:461-470`).

- [x] **[FIXED] Rule of Three / Five violation in `ResponseLayer` (`responselayer.h:21-60`)**
  - **Issue**: Manages raw heap arrays (`responses`, `cornerResponses`, `laplacian`, `isblob`) with `new[]` and `delete[]`, but has no copy constructor or copy assignment operator.
  - **Impact**: Copying an instance performs shallow pointer copying, leading to double-free errors upon destruction.
  - **Fix**: Delete copy operations (`ResponseLayer(const ResponseLayer&) = delete;`) or replace raw pointers with `std::vector<float>`, `std::vector<uint8_t>`, and `std::vector<bool>`.
  - **Resolved**: Copy ctor and copy assignment are deleted (`responselayer.h:35-36`).

- [ ] **[STILL OPEN (minor)] Per-candidate heap allocations in `interpolateStep` (`fasthessian.cxx:631-657`)**
  - **Issue**: `deriv4D()` and `hessian4D()` dynamically allocate `new cv::Matx41d()` and `new cv::Matx44d()`.
  - **Impact**: Unnecessary heap allocation and deallocation overhead for thousands of candidate points.
  - **Fix**: Return `cv::Matx41d` and `cv::Matx44d` by value on the stack.

---

## 3. Algorithmic & Mathematical Inconsistencies (Medium Priority)

- [ ] **[STILL OPEN] Rigid registration in `MatchPoint::getRT` hardcodes `R = Identity` (`MatchPoint.cxx:323-344`)**
  - **Issue**: The Umeyama SVD rotation calculation is commented out. `Transform.get_rotation() = R;` assigns identity.
  - **Impact**: `match3d` cannot estimate 3D rotations, only isotropic scale and translation.
  - **Fix**: Uncomment and validate the SVD-based rotation estimator, ensuring proper reflection handling when $\det(U V^T) < 0$.

- [ ] **[STILL OPEN] Scale divisor check in `MatchPoint::getRT` (`MatchPoint.cxx:321`)**
  - **Issue**: `double s = sqrt((EigCa * EigCb) / (EigCa * EigCa));`
  - **Impact**: If points in set A are collinear or identical, `EigCa * EigCa` can be 0, causing division by zero.
  - **Fix**: Guard against near-zero denominator before computing scale.

- [x] **[FIXED] `EXTRACT_CORNER` formula typos (`fasthessian.cxx:434, 437, 478`)**
  - **Issue**: In gradient box dimensions, `l*l-1` is passed instead of `2*l-1` for $I_y$ and $I_z$. At line 478, `cvGet2D(&mat, 0, ...)` reads input matrix elements instead of computed eigenvalues.
  - **Fix**: Correct filter dimensions to `2*l - 1` and read eigenvalues from `eval`.
  - **Resolved**: Filter sizes now `2*l-1` for `Iy`/`Iz` (`fasthessian.cxx:439-443`); eigenvalues read from `cv::eigen` output `eigs` (`fasthessian.cxx:466-470`).

- [ ] **[STILL OPEN] Inconsistent coordinate conventions across function signatures**
  - **Issue**: Some methods use `(column, row, layer)` = $(x, y, z)$, while others pass `(r, c, d)` = $(y, x, z)$ or `(d, r, c)` = $(z, y, x)$ (e.g., `interpolateExtremum(d, r, c, ...)` vs `isExtremum(r, c, d, ...)`).
  - **Impact**: Highly error-prone for future maintainers.
  - **Fix**: Adopt a consistent parameter naming and ordering convention (e.g. `x, y, z` or `row, col, slice`).

---

## 4. Safety, Robustness & Error Handling (Medium Priority)

- [ ] **[PARTIAL] Unchecked file operations (`fopen`, `gzopen`, `ifstream`)**:
  - `WritePointsBinary()` (`vtk3DSURF.cxx:504`): Calls `fopen(fileName, "wb")` without checking for `NULL` before `fwrite()`.
  - `WritePointsCSVGZ()` (`vtk3DSURF.cxx:450`): Calls `gzopen()` without checking for `NULL` before `gzprintf()`.
  - `ReadIPoints()` (`vtk3DSURF.cxx:39`): Does not verify `file.is_open()`. `std::stof` throws unhandled exceptions on malformed CSV rows.
  - **Fix**: Add null checks, log errors, and wrap parsing with `try / catch (const std::exception&)`.
  - **Partial**: Null checks added for `fopen` (`vtk3DSURF.cxx:531-534`), `gzopen` (`vtk3DSURF.cxx:472-475`), and `ReadIPoints` now verifies `is_open()` (`vtk3DSURF.cxx:42-45`). Still open: `ReadIPoints` `std::stof` is not wrapped in `try/catch`.

- [ ] **[STILL OPEN] Premature overwriting of output JSON files (`surf3d.cxx:295-318`)**
  - **Issue**: `surf3d.cxx` unconditionally writes image bounding box metadata to `outfilename + ".json"`. If `-json 1` is requested, `SURF->WritePoints()` immediately overwrites that same file with the keypoint list.
  - **Fix**: Merge bounds metadata into the keypoints JSON object, or write bounds to a dedicated file (e.g., `<basename>_bounds.json`).

- [ ] **[STILL OPEN] Global mutable counters in `fasthessian.cxx:24-25, 574`**
  - **Issue**: `int nb_pts = 0; int nb_corner_pts = 0; int fsdi = 0;` are non-const global variables.
  - **Impact**: Not thread-safe; counters accumulate across multiple calls to `FastHessian` within the same process.
  - **Fix**: Convert them to class member variables or local return values.

- [ ] **[STILL OPEN] Operator `-` semantics in `Ipoint` (`ipoint.h:35-45`)**
  - **Issue**: `float operator-(const Ipoint &rhs)` computes the Euclidean distance between descriptor vectors.
  - **Impact**: Overloading `operator-` to return a scalar distance is counter-intuitive in C++ (where subtraction usually returns a difference vector or offset). Furthermore, it does not verify that `this->descriptor.size() == rhs.descriptor.size()`.
  - **Fix**: Replace with a named method `float distanceTo(const Ipoint& rhs) const` and assert matching descriptor lengths.

---

## 5. Code Quality, Modern C++ & Refactoring (Medium Priority)

- [x] **[FIXED] Namespace pollution in header files**:
  - `using namespace picojson;` and `using namespace TooN;` in `MatchPoint.h:13-14`.
  - `using namespace std;` in `vtk3DSURF.h:9`.
  - `using std::cout; using std::endl;` in `vtkRobustImageReader.h:22-23`.
  - **Fix**: Remove `using namespace` directives from all `.h` header files and qualify symbols explicitly.
  - **Resolved**: Removed from `vtk3DSURF.h`, `vtkRobustImageReader.h`, and `MatchPoint.h` (still present in some `.cxx` files, which is acceptable).

- [ ] **[STILL OPEN] Replace global macro definitions with inline functions**:
  - `integ(x,y,z)` and `sourc(x,y,z)` in `integral.h:6-7` pollute the global preprocessor namespace and rely on implicitly named local pointers (`pin_tegral`, `int_egral_incs`).
  - `PRINT`, `CHECKMAT`, `PRINTIF` in `fasthessian.cxx`.
  - **Fix**: Replace macros with typed inline helper functions or lambda accessors.

- [ ] **[STILL OPEN] Remove dead / commented-out code**:
  - `Surf::ThreadedDesc` declared in `surf.h:41` but never implemented.
  - `FastHessian::EigenValue` defined in `fasthessian.cxx:487-520` but never used.
  - `FastHessian::FittingQuadric` (~200 lines) dead code in `fasthessian.cxx`.
  - Unused variables generating compiler warnings (`imin2`, `temp2`, `test`, `bounds`, `img_pointer`, `layer_min`, `layer_max`).

- [ ] **[STILL OPEN] Eliminate compiler warnings under `-Wall -Wextra`**:
  - Constructor member initialization order warnings (`-Wreorder`) in `Ipoint`, `FastHessian`.
  - Signed vs. unsigned comparisons (`-Wsign-compare`) across loops in `vtk3DSURF.cxx`, `surf.cxx`, `fasthessian.cxx`.
  - Unused parameters (`-Wunused-parameter`) in `BoxIntegralOptim` and `buildResponseLayer`.

---

## 6. Build System, CI & Repository Hygiene (Low Priority)

- [x] **[FIXED] Add `.gitignore`**
  - **Issue**: The repository has no `.gitignore`. In-tree build outputs (`CMakeCache.txt`, `CMakeFiles/`, `Makefile`, `surf3d`, `match3d`, `*.csv`, `*.gz`, `*.json`) are currently tracked as untracked files by Git.
  - **Fix**: Create a standard `.gitignore` covering CMake artifacts, build directories, IDE files, and data outputs (`points*`, `transform.json`).
  - **Resolved**: `.gitignore` added (CMake artifacts, binaries, data outputs, IDE files, `test_results/`).

- [ ] **[PARTIAL] Modernize `CMakeLists.txt`**:
  - Move `cmake_minimum_required(VERSION 3.20)` to line 1 outside any `if()` blocks (currently requires 3.30 inside an `if` statement).
  - Explicitly set `set(CMAKE_CXX_STANDARD 17)` and `set(CMAKE_CXX_STANDARD_REQUIRED ON)`.
  - Modernize target linking: use `target_link_libraries(surf3d PRIVATE VTK::... OpenMP::OpenMP_CXX ZLIB::ZLIB ${OpenCV_LIBS})` instead of global `include_directories` and modifying `CMAKE_CXX_FLAGS`.
  - Add `find_package(LAPACK)` / `find_package(BLAS)` and TooN include path resolution when `BUILD_SURFMATCH=ON`.
  - Add `install(TARGETS surf3d RUNTIME DESTINATION bin)`.
  - **Partial**: `CMAKE_CXX_STANDARD 17` / `REQUIRED ON` set and `target_link_libraries(... PRIVATE ...)` used for `surf3d`. Still open: `cmake_minimum_required` remains inside the `if()` block (line 2), global `include_directories` still used, no `find_package(LAPACK/BLAS)`/TooN resolution, no `install()`.

- [ ] **[PARTIAL] Update GitHub Actions CI (`.github/workflows/ci.yml`)**:
  - Upgrade `actions/checkout@v2` to `@v4` to prevent deprecation warnings.
  - Add a testing step with sample image data to verify feature detection outputs (`ctest`).
  - **Partial**: A `./test.sh -v` testing step was added (ci.yml:63-67). Still open: `actions/checkout@v2` not yet upgraded to `@v4`.

---

## 7. Documentation & CLI Consistency (Low Priority)

- [ ] **[STILL OPEN] Document missing command-line options in `Readme.md` and usage printout**:
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

## 8. Remediation Plan (Proposed Work Order)

Ordered by priority. Each item lists the current `file:line` and the proposed fix.

### 8.1 Critical / correctness (do first)
- [ ] **Fix NIFTI rescale bypass** (`vtkRobustImageReader.h:118-120`) — move the
  `GetRescaleSlope`/`GetRescaleIntercept` block (lines 138-156) before or outside
  the early `return`, so scaling applies even when no axis flip is needed.
- [ ] **Fix scale-space limit integer division** (`fasthessian.cxx:206-209`) —
  invert `(int)(t->width/b->width)` to `(int)(b->width/t->width)`, or use
  floating-point scale factors, so `FIRST_SCALE` is applied correctly.
- [ ] **Add include guard to `responselayer.h`** (`responselayer.h:1-22`) —
  wrap in `#ifndef RESPONSELAYER_H` / `#define RESPONSELAYER_H` / `#endif`.
- [ ] **Fix unsigned underflow in box integrals** (`integral.h:63,109-116`) —
  cast operands to signed `long long` before subtraction in both
  `BoxIntegral` and `BoxIntegralOptim`.
- [ ] **Guard `getRT` scale divisor** (`MatchPoint.cxx:321`) — skip/bail when
  `EigCa*EigCa` is near zero.

### 8.2 Matching correctness
- [ ] **Re-enable rotation in `getRT`** (`MatchPoint.cxx:323-344`) — uncomment and
  validate the SVD/Umeyama rotation with proper reflection handling
  (`det(U V^T) < 0`).

### 8.3 Safety / robustness
- [ ] **Wrap `ReadIPoints` parsing** (`vtk3DSURF.cxx:54-77`) — wrap `std::stof` in
  `try/catch` and skip malformed rows.
- [ ] **Fix JSON overwrite** (`surf3d.cxx:295-318`) — merge bounds into the
  keypoints JSON, or write bounds to `<basename>_bounds.json`.
- [ ] **Convert global counters to members** (`fasthessian.cxx:24-25,574`) —
  make `nb_pts`/`nb_corner_pts`/`fsdi` class members or local values.
- [ ] **Rename `Ipoint::operator-`** (`ipoint.h:35-45`) — replace with a named
  `distanceTo(const Ipoint&) const` and assert matching descriptor lengths.

### 8.4 Cleanup / minor
- [ ] **Remove dead code** (`surf.h:41` `ThreadedDesc`; `fasthessian.cxx:479`
  `EigenValue`; `fasthessian.cxx:716` `FittingQuadric`; unused variables).
- [ ] **Return matrices by value** (`fasthessian.cxx:623-624,665,hessian4D`) —
  have `deriv4D`/`hessian4D` return `cv::Matx41d`/`cv::Matx44d` on the stack.
- [ ] **Replace macros** (`integral.h:6-7` `integ/sourc`; `fasthessian.cxx`
  `PRINT`/`CHECKMAT`/`PRINTIF`) with inline functions.
- [ ] **Triage `-Wall -Wextra` warnings** — `-Wreorder`, `-Wsign-compare`,
  `-Wunused-parameter` across the sources.

### 8.5 Build / CI / docs (low priority)
- [ ] **Modernize `CMakeLists.txt`** — move `cmake_minimum_required` to line 1,
  replace global `include_directories`, add `find_package(LAPACK/BLAS)` + TooN
  path for `BUILD_SURFMATCH`, add `install(TARGETS surf3d ...)`.
- [ ] **Upgrade CI** (`ci.yml:28`) — `actions/checkout@v2` → `@v4`.
- [ ] **Document missing CLI options** — add `-d`, `-m`, `-nt`, `-p`, `-pad`,
  `-gz`, `-precision` to `Readme.md` option table.

