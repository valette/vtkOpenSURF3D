# Architecture

> **Last updated:** 2026-10-07

vtkOpenSURF3D is a C++17 project that detects 3D scale-space features
(keypoints) in volumetric images and computes a local descriptor for each of
them. It is split into **two executables** with overlapping, but distinct,
responsibilities:

- **`surf3d`** — loads an image, runs the detector (`FastHessian`) and the
  descriptor (`Surf` / sub-volume extraction), and writes the keypoints in
  several output formats. Built by default.
- **`match3d`** — reads two JSON keypoint sets, matches their descriptors,
  runs a RANSAC loop to estimate a similarity transform (SIM3), and writes a
  transform JSON. Optional (`BUILD_SURFMATCH=ON`, needs TooN + LAPACK/BLAS).

Both executables share the core data structure `Ipoint` (`ipoint.h`) and the
JSON writer `picojson` (`picojson.h`, vendored single header).

---

## `surf3d` pipeline

The processing flow is orchestrated in `surf3d.cxx::main` and
`vtk3DSURF.cxx::vtk3DSURF::Update`. The stages, in order:

```
 surf3d.cxx                    vtk3DSURF::Update()
 -------------                 ---------------------
 input image file
      │
      ▼
 vtkRobustImageReader::Update()        ── reads file, applies axis flips and
      │                                   NIFTI rescale slope/intercept
      ▼
 optional preprocessing in surf3d.cxx  ── mirror pad (-pad), clamp (-cmin/-cmax),
      │                                   mask loading (-m, dimension check)
      ▼
 vtk3DSURF::Update()
      │
      ├─ 1. color→luminance conversion (if ≥3 components)
      ├─ 2. optional isotropic resampling (-s spacing / -d maxsize)
      ├─ 3. vtkImageCast → int, then vtkImageShiftScale by -Range[0]  (Cast)
      ├─ 4. ComputeIntegral(Cast)            → Integral (uint64 integral image)
      ├─ 5a. if -p pointfile: ReadIPoints()  → points from CSV
      │       else: FastHessian(Integral,...).getIpoints()  → detect keypoints
      ├─ 6. optional top-N selection (-n)
      ├─ 7. descriptors by DescriptorType:
      │        0 → Surf::getDescriptors      (48-dim SURF3D)
      │        1 → Surf::getRawDescriptors   (3·8·r³ HAAR coefficients)
      │        2 → ThreadedSubVolumes        (8·r³ raw voxels)
      └─ 8. return to surf3d.cxx → write JSON/CSV/BIN/CSV.GZ + bounds JSON
```

Key invariants:

- **Integral image is the single source of truth** for both detection and
  description. It stores `unsigned long long` sums so that box integrals do
  not overflow.
- **Keypoint coordinates and scale are stored in *voxel* units** relative to
  the processed (`Cast`) image. All exporters convert back to physical space
  using `Cast` origin/spacing (see [Output formats](output-formats.md)).
- The detector operates on **`Cast`** — the integer, zero-shifted image — not
  the raw input.

---

## Threading model

Three distinct parallelism mechanisms are used:

1. **OpenMP** (`#pragma omp parallel for`)
   - `buildResponseLayer` (`fasthessian.cxx:384`) — parallelizes over the Z
     axis of each response layer. This is the hot loop of detection.
   - `Surf::getDescriptors` / `getRawDescriptors` (`surf.cxx:42,53`) —
     parallelizes over keypoints.
   - `ThreadedIntegral` (`integral.cxx:90,101,111`) — parallelizes the three
     prefix-sum passes.
   - Enabled via `find_package(OpenMP)` in CMake. The `-nt` CLI option only
     controls **VTK** filters (`vtkImageResample`, `vtkImageCast`,
     `vtkImageShiftScale`, and the `vtkMultiThreader` paths); it does *not*
     change OpenMP thread counts.

2. **VTK `vtkMultiThreader`**
   - `vtk3DSURF::ThreadedSubVolumes` (`vtk3DSURF.cxx:274`) — distributes
     keypoints across threads for descriptor type 2 (raw sub-volume).

3. **Serial**
   - Matching (`match3d`) and the RANSAC loop are serial.

> **Note:** The OpenMP loops write to disjoint `responses[index]` /
> `ipts[id]` slots, so there is no shared-memory race in normal operation.
> The global counters `nb_pts` / `nb_corner_pts` in `fasthessian.cxx` are,
> however, incremented inside `interpolateExtremum` (serial section), so they
> are only informational.

---

## Directory / file map

| File | Role |
|------|------|
| `surf3d.cxx` | `surf3d` CLI: arg parsing, preprocessing, output writing |
| `vtk3DSURF.h/.cxx` | `vtk3DSURF` VTK object: pipeline orchestration + writers |
| `fasthessian.h/.cxx` | `FastHessian`: scale-space detection (DoH) |
| `responselayer.h` | `ResponseLayer`: response/laplacian/blob storage |
| `integral.h/.cxx` | integral image computation + `BoxIntegral` helpers |
| `surf.h/.cxx` | `Surf`: Haar-wavelet descriptors |
| `ipoint.h` | `Ipoint` keypoint struct + descriptor distance |
| `vtkRobustImageReader.h` | robust image reader (flips, NIFTI rescale) |
| `MainMatch.cxx` | `match3d` CLI |
| `MatchPoint.h/.cxx` | matching, RANSAC, SIM3 estimation |
| `picojson.h` | vendored JSON parser/writer (single header) |
| `CMakeLists.txt` | build definition |
| `test.sh` | automated test suite |

See [Modules](modules.md) for a detailed per-class walkthrough.
