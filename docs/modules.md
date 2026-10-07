# Modules — Code Reference

> **Last updated:** 2026-10-07

This page documents each source file and its key classes/functions, with
line references (relative to the current source tree). It is intended as a
navigational aid for developers and coding agents.

---

## `ipoint.h` — keypoint data structure

The central data structure shared by all modules.

- `class Ipoint` (line 27):
  - Fields: `x, y, z` (voxel coordinates, `float`), `scale` (`float`),
    `response` (`float`), `laplacian` (`int`), `descriptor`
    (`std::vector<float>`).
  - `operator-` (line 35): Euclidean distance between two descriptors
    (sum over the *common* prefix length, `std::min`).
  - `allocate(int size)` (line 47): `descriptor.resize(size)` (zero-fills).
- `typedef std::vector<Ipoint> IpVec` (line 22).
- `typedef std::vector<std::pair<int,int>> IpPairVec` (line 23).

> **Units:** `x/y/z` and `scale` are in **voxels** of the processed `Cast`
> image. Physical-space conversion happens only at export time.

---

## `vtkRobustImageReader.h` — image loading

A `vtkObject` that reads a NIFTI or MHD image and normalizes its geometry.

- `Update()` (line 32):
  1. Uses `vtkImageReader2Factory` with registered `vtkMetaImageReader` and
     `vtkNIFTIImageReader`.
  2. For MHD: parses `TransformMatrix`/`Orientation`/`Rotation` from the
     header to detect axis flips (negative diagonal element → flip).
  3. For NIFTI: reads `GetQFormMatrix()` (falling back to `GetSFormMatrix()`),
     sets the origin, and handles `GetQFac()`.
  4. If any flip is needed, applies `vtkImageFlip` per axis and fixes the
     origin (line 125–136).
  5. For NIFTI with a rescale slope/intercept (`GetRescaleSlope()` /
     `GetRescaleIntercept()`), applies `vtkImageShiftScale` (line 138–156).

> **Known issue:** the early `return` at line 118–120 (when no flip is
> needed) bypasses the NIFTI rescale step. See
> [TODO](../TODO.md) (H1).

---

## `integral.h` / `integral.cxx` — integral image

- `ComputeIntegral(vtkImageData*)` (`integral.cxx:11`): builds a
  `VTK_UNSIGNED_LONG_LONG` integral image by three sequential prefix-sum
  passes (X, Y, Z), each OpenMP-parallelized (`ThreadedIntegral`,
  `integral.cxx:71`).
- `BoxIntegral(img, dim0, dim1, dim2, size0, size1, size2)`
  (`integral.h:31`): **bounds-clamped** box sum using 8 corner reads of the
  integral image; returns a non-negative `unsigned long long`. Safe to call
  anywhere (used by the non-optimized Haar functions).
- `BoxIntegralOptim(...)` (`integral.h:68`): **unchecked** fast variant using
  raw increments and pointer arithmetic; used in the hot loops
  (`buildResponseLayer`, `Surf::*Optim`). Assumes the box lies inside the
  image (guaranteed by the caller's limit computations).
- Macros `integ(x,y,z)` / `sourc(x,y,z)` (line 6–7) index into the integral /
  source buffers using VTK increments.

> **Caution:** `BoxIntegralOptim` performs **no bounds clamping**. If a box
> extends outside the image it reads out of bounds. The caller is responsible
> for staying in-bounds. See [TODO](../TODO.md) (M1).

---

## `responselayer.h` — response layers

`class ResponseLayer` (line 22) stores one scale layer of the scale-space:

- Fields: `width, height, depth, step, filter`, plus parallel buffers
  `responses` (`float*`), `laplacian` (`unsigned char*`),
  `cornerResponses` (`float*`), `isblob` (`bool*`).
- Constructor (line 38) allocates and **zero-initializes** all buffers with
  `new float[total]()` etc.
- Indexing is `column + row*width + layer*width*height`
  (i.e. `x + y*width + z*width*height`).
- Accessors have two overloads: direct (`getResponse(r,c,l)`) and
  scale-adapted (`getResponse(r,c,l,src)`) where `scale = this->width/src->width`
  maps a coarser layer index to a finer layer.

> **Known issue:** no include guard. See [TODO](../TODO.md) (L1).

---

## `fasthessian.h/.cxx` — scale-space detection

`class FastHessian` performs blob detection via the determinant of Hessian
(DoH) over a multi-octave scale space. See [Algorithm](algorithm.md) for
details. Key members:

- Constants (`fasthessian.h:28-31`): `OCTAVES=5`, `INTERVALS=4`,
  `THRES=0.0004f`, `INIT_SAMPLE=2`.
- Constructor with image (line 53) adjusts the number of octaves based on the
  smallest image dimension (line 65–84) so that small images do not require
  non-existent layers.
- `getIpoints()` (`fasthessian.cxx:141`):
  1. `buildResponseMap()` builds the response layers.
  2. Iterates `(octave, interval)` triplets `(b, m, t)` and searches for
     scale-space maxima (`isExtremum`).
  3. Refines each extremum by 4D sub-pixel interpolation
     (`interpolateExtremum` → `interpolateStep`, using `deriv4D` and
     `hessian4D` with SVD).
  4. Applies the mask (if any) to remove points outside the mask.
- `buildResponseMap()` (`fasthessian.cxx:287`): allocates `ResponseLayer`s for
  octaves 1..5 with filter sizes 9/15/21/27, 39/51, 75/99, 147/195, 291/387.
- `buildResponseLayer(id, decile)` (`fasthessian.cxx:348`): computes the DoH
  approximation (`Dxx, Dyy, Dzz, Dxy, Dyz, Dxz`) from box integrals, plus
  `Sdet2p`, `Trace`, `Det`, the blob/laplacian flags. OpenMP-parallelized
  over the Z axis. Optional corner response under `EXTRACT_CORNER`.
- `isExtremum` / `isCornerExtremum` (line 515/543): 3×3×3×3 non-maximum
  suppression across `(b, m, t)` layers.
- `interpolateExtremum` (line 577): accepts a point if the 4D offset is < 1
  voxel/filter-step; pushes an `Ipoint` with `x,y,z = (c+xX, r+xY, d+xZ)*step`
  and `scale = 0.1333*(m->filter + xS*filterStep)`.
- `interpolateStep` (line 618): Newton step `X = -H⁻¹ ∇` computed via SVD
  with eigenvalue truncation.
- `deriv4D` (line 665), `hessian4D` (line 949): finite-difference 1st and
  2nd derivatives in `(x,y,z,s)`.
- `FittingQuadric` (line 716): a full quadric fit alternative, **not used**
  (`interpolateExtremum` calls `interpolateStep`; the `FittingQuadric` call
  is commented out at line 588).
- `WriteResponseMap()` (line 1040): debug helper that writes each response
  layer as a Meta image under `FH/`.

> **Globals:** `nb_pts` / `nb_corner_pts` (`fasthessian.cxx:24-25`) and
> `fsdi` (line 574) are file-scope globals used only for diagnostics.
> See [TODO](../TODO.md) (L3).

---

## `surf.h/.cxx` — Haar descriptors

`class Surf` computes the keypoint descriptors from the integral image.

- `getDescriptors(radius=5, normalize=true)` (line 40): 48-dimensional SURF3D
  descriptor. OpenMP-parallelized over keypoints.
- `getRawDescriptors(radius=5)` (line 51): `3·8·radius³` raw Haar-coefficient
  descriptor. OpenMP-parallelized.
- `getDescriptor` (line 63): divides the neighbourhood into 8 sub-blocks
  (2 per axis); for each sub-block accumulates the Haar X/Y/Z responses and
  their absolute values (6 values per sub-block → 48 total). Optionally
  normalizes the vector to unit length.
- `getRawDescriptor` (line 161): writes raw `haarX/Y/ZOptim` values at every
  sample (3 per voxel).
- `gaussian` (line 227): denormalized 3D Gaussian weighting, used to weight
  Haar responses (spread `2.5*scale`).
- `haarX/Y/Z` (lines 235–257): box-integral Haar wavelets using the
  **clamped** `BoxIntegral`.
- `haarX/Y/ZOptim` (lines 262–284): same wavelets using the **unchecked**
  `BoxIntegralOptim` for speed (used in the descriptor loops).

> **Caution:** the optimized Haar functions read unclamped boxes whose
> coordinates are computed from `ipt->x + u*scale` etc. For keypoints near
> the image border these can go out of bounds. See [TODO](../TODO.md) (H4).

---

## `vtk3DSURF.h/.cxx` — pipeline orchestration

`class vtk3DSURF : public vtkObject` glues everything together and writes
outputs.

- Constructor defaults (line 85): `Normalize=true`, `Threshold=0.0004`,
  `NbThread=-1`, `MaxSize=100`, `DescriptorType=0`, `NumberOfPoints=-1`,
  `Spacing=0`, `SubVolumeRadius=5`.
- `Update()` (line 81): runs the full pipeline (see
  [Architecture](architecture.md)). Notably:
  - Color→luminance conversion (line 92) if ≥3 scalar components.
  - Resampling (line 105–158) when `Spacing` or `MaxSize` is set (also
    resamples the mask with nearest-neighbour interpolation).
  - Cast to `int` (`vtkImageCast`) and zero-shift (`vtkImageShiftScale` by
    `-Range[0]`) producing `Cast` (line 160–178).
  - `ComputeIntegral` (line 181).
  - Detection (`FastHessian`, line 191) or `ReadIPoints` from `-p` (line 185).
  - Top-N selection (line 207).
  - Descriptor dispatch (line 225–244): type 0/1 via `Surf`, type 2 via
    `ThreadedSubVolumes`.
  - Prints descriptor range (line 249–258).
- `ReadIPoints()` (line 35): reads a **CSV** point file (the `-p` option),
  converting physical coords back to voxels using `Cast` origin/spacing, and
  scale to voxel units via `Sspacing`.
- `ThreadedSubVolumes(void*)` (line 274): static VTK thread function for
  descriptor type 2. For each keypoint it crops a `2*radius` cube from
  `Resized` (via `vtkImageResize`) around the point, resampled to
  `2*radius` voxels, and stores the `8*radius³` raw voxels as the descriptor.
- `WritePoints(fileName)` (line 349): JSON keypoints + bounds + origin +
  spacing + dimensions.
- `WritePointsCSV` (line 417), `WritePointsCSVGZ` (line 456, uses zlib),
  `WritePointsBinary` (line 518).
- `compareResponses` (line 33): sort predicate (descending response).

---

## `surf3d.cxx` — `surf3d` CLI

`main()` (line 19):

- Usage/help text (line 22–49).
- Argument parsing loop (line 77–167). See [Usage](usage.md) for the option
  list. Note the parser reads `argv[argumentsIndex+1]` without a bounds check
  on the last argument (see [TODO](../TODO.md)).
- Loads image via `vtkRobustImageReader` (line 174–181), returns exit code 5
  on failure.
- Optional mask loading (line 184–218) with dimension/component checks
  (exit 6 on mismatch).
- Mirror padding (line 220–239) via `vtkImageMirrorPad` + `vtkImageTranslateExtent`.
- Intensity clamping (line 241–262) via `vtkImageThreshold`.
- Configures `vtk3DSURF` and calls `Update()` (line 267–293).
- Writes the **bounds JSON** (`<out>.json`, line 296–309) — this is the file
  used by `match3d` for its `bounds` field.
- Writes JSON/BIN/CSV/CSV.GZ per flags (line 311–345).

---

## `MainMatch.cxx` — `match3d` CLI

`main()` (line 6):

- Requires ≥ 2 keypoint JSON files.
- Parses options: `-i` (write inliers), `-Rd` (RANSAC dist, default 50),
  `-Ri` (RANSAC min inliers, default 2), `-Md` (matching dist, default 0.95),
  `-Md2` (second-best ratio, default 0.98), `-Ms` (matching scale, default
  1.5), `-b` (compute bounding boxes), `-bb` (base64-encoded bbox JSON).
- Parses both JSON files, optionally a bounding box, computes matches,
  computes the transform, and writes `transform.json` (hard-coded filename
  in the current working directory).

---

## `MatchPoint.h/.cxx` — matching & transform

`class MatchPoint` (see [SurfMatch](surfmatch.md) for algorithm details):

- `Parse(fileName, id)` (`MatchPoint.cxx:16`): reads a JSON keypoint file into
  `points[id]`.
- `computeMatches()` (line 62): descriptor matching with laplacian gating,
  scale-ratio gating, best/second-best ratio test.
- `computeTransform()` (line 157): RANSAC over `matches`; each hypothesis
  picks `N=3` matches, calls `getRT`, counts inliers, and keeps the best
  transform.
- `getRT(pairs, Transform)` (line 271): estimates a SIM3 (translation,
  rotation, scale) via an Umeyama-style least-squares fit. **The rotation is
  currently hard-coded to `Identity`** (line 344) — the SVD-based rotation
  estimation is commented out (line 323–343).
- `WriteTransform(fileName, writeInliers)` (line 353): writes
  `transform.json` with translation/rotation/scale/inliers/fail/bboxes.
- `BboxParse` (line 426): parses a base64-encoded JSON bounding box for
  constraining matching.
- `Box3` (line 499): axis-aligned bounding box helper.
- `base64_decode` (line 535): helper for `-bb`.

---

## `picojson.h` — JSON

Vendored single-header JSON parser/serializer (public domain style).
Used for keypoint JSON export (`vtk3DSURF::WritePoints`), JSON import
(`MatchPoint::Parse`), and the bounds file. No external dependency.
