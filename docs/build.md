# Build

> **Last updated:** 2026-10-07

vtkOpenSURF3D is built with CMake. It produces two executables:

- `surf3d` — built by default.
- `match3d` — optional (requires TooN + LAPACK/BLAS).

---

## Dependencies

| Dependency | Required for | Notes |
|------------|--------------|-------|
| CMake ≥ 3.30 | both | `cmake_minimum_required(VERSION 3.30)` |
| C++17 compiler | both | GCC/Clang |
| VTK | both | components `CommonSystem`, `ImagingColor`, `IOExport` |
| OpenCV (core) | `surf3d` | `find_package(OpenCV REQUIRED COMPONENTS core)` |
| ZLIB | `surf3d` | gzipped CSV output |
| OpenMP | both | optional, auto-detected (`find_package(OpenMP)`) |
| TooN | `match3d` | only when `BUILD_SURFMATCH=ON` |
| LAPACK + BLAS | `match3d` | only when `BUILD_SURFMATCH=ON` |

> **Note:** `match3d` also links VTK, so VTK is needed for both targets.

---

## Basic build (surf3d only)

```bash
git clone https://github.com/valette/vtkOpenSURF3D.git
cd vtkOpenSURF3D
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

The binary is `build/surf3d`.

## Build with match3d

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_SURFMATCH=ON
cmake --build build --parallel
```

> **Known issue:** the committed `match3d` binary at the repository root is
> stale (built 2024-12-05) and may not link against the currently installed
> VTK/TooN. Prefer a fresh build with `BUILD_SURFMATCH=ON`.

---

## CMake options

| Option | Default | Meaning |
|--------|---------|---------|
| `BUILD_SURFMATCH` | `OFF` | Build the `match3d` executable |
| `EXTRACT_CORNER` | `OFF` | Define `EXTRACT_CORNER` to enable corner detection code paths |
| `BUILD_TESTING` | - | Used by CI; no CTest targets are defined in-tree |

The CMake file supports being included as a **subdirectory** (it guards the
project-level calls with
`if( CMAKE_PARENT_LIST_FILE STREQUAL CMAKE_CURRENT_LIST_FILE )`,
`CMakeLists.txt:1`). When included from a parent project it defines the
`surf3d`/`match3d` targets but not a new project.

---

## Recommended build for development

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_FLAGS="-Wall -Wextra"
cmake --build build --parallel -j
```

Adding `-Wall -Wextra` is strongly recommended: the codebase compiles with a
number of warnings (unused variables, implicit conversions) that signal latent
issues. See [TODO](../TODO.md) for details.

---

## Building on the CI image

The repository's CI (`ci.yml`) builds inside the
`kitware/vtk-for-ci:<version>` containers, installs `libopencv-dev`, and runs:

```bash
mkdir build
cmake -S . -B build -DCMAKE_PREFIX_PATH=/opt/vtk/install/ \
      -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel 2
```
