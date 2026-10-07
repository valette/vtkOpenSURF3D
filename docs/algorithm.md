# Algorithm — 3D SURF Detection and Description

> **Last updated:** 2026-10-07

This page explains the computer-vision algorithm implemented by vtkOpenSURF3D.
It is a 3D generalization of the 2D SURF method (Bay et al., 2008), following
the SURF3D paper by Agier et al. (2016).

---

## 1. Overview

SURF detects **blob-like** interest points that are stable across scale, then
characterizes each with a **rotation-invariant-free** (upright) local
descriptor. The pipeline has four conceptual stages:

1. **Integral image** — enables O(1) box-sum queries.
2. **Scale-space construction** — a pyramid of response layers, each holding
   an approximation of the determinant of the Hessian (DoH).
3. **Extrema detection** — find local maxima of the DoH in the 3D+scale space.
4. **Sub-pixel refinement & description** — refine coordinates and compute
   Haar-wavelet descriptors.

vtkOpenSURF3D does **not** estimate orientation; it produces *upright* SURF
features (the optional orientation computation is removed). The output keypoint
therefore carries position, scale, sign of the Laplacian, and descriptor, but
no orientation.

---

## 2. Integral image

Given the (integer, zero-shifted) image `I`, the integral image `S` is:

```
S(x,y,z) = Σ_{i≤x, j≤y, k≤z} I(i,j,k)
```

computed by three sequential prefix-sum passes over X, Y, Z
(`ComputeIntegral`, `integral.cxx:11`). Any axis-aligned box sum can then be
computed from 8 corner samples:

```
BoxSum = S(x2,y2,z2) - S(x2,y2,z1) - S(x2,y1,z2) - S(x1,y2,z2)
       + S(x1,y1,z2) + S(x1,y2,z1) + S(x2,y1,z1) - S(x1,y1,z1)
```

This is implemented as `BoxIntegral` (clamped) and `BoxIntegralOptim`
(unchecked, fast). Storing 64-bit sums prevents overflow.

---

## 3. Scale-space construction

SURF builds a scale-space by applying **box filters of increasing size** to
the integral image, rather than resizing the image and using a fixed filter.
This is faster and equivalent to Gaussian scale-space for blob detection.

### Response layers

A `ResponseLayer` is defined by a `filter` size (in voxels) and a `step`
(subsampling). The filter sizes grow roughly geometrically across octaves:

| Octave | Filters            |
|--------|--------------------|
| 1      | 9, 15, 21, 27      |
| 2      | 39, 51             |
| 3      | 75, 99             |
| 4      | 147, 195           |
| 5      | 291, 387           |

Each successive octave also halves the resolution (`step *= 2`,
dimensions / 2). The number of octaves actually used is limited by the image
size (`fasthessian.cxx:65-84`).

### Approximated determinant of Hessian

For each voxel the Hessian is approximated from box integrals:

- `Dxx, Dyy, Dzz`: second derivatives along each axis.
- `Dxy, Dyz, Dxz`: mixed partial derivatives.

The blob response is the (normalized) determinant:

```
Sdet2p = Dyy*Dzz + Dxx*Dyy + Dxx*Dzz − 0.8330*(Dxy² + Dxz² + Dyz²)
Det    = Dxx*Dyy*Dzz + 2*0.7603*Dxy*Dyz*Dxz
         − 0.8330*(Dxx*Dyz² + Dyy*Dxz² + Dzz*Dxy²)
response = |Det| / filterVolume³
```

A voxel is flagged as a **blob** if `Sdet2p > 0` and `Trace*Det > 0`
(`fasthessian.cxx:452`). The **Laplacian sign** is `sign(Dxx+Dyy+Dzz)`
(line 458) and is used later for fast matching.

---

## 4. Extrema detection (non-maximum suppression)

`getIpoints()` walks triples of adjacent response layers `(b, m, t)` where
`b` = bottom (finer/wider), `m` = middle, `t` = top (coarser/narrower). For
each voxel `(r, c, d)` in the middle layer it performs a **3×3×3×3** search:

- Compare the candidate response against the 3×3×3 neighbours in `b`, `m`,
  and `t`.
- The candidate is kept only if it is strictly greater than all of them
  (subject to the `FIRST_SCALE` / `LAST_SCALE` boundary conditions that relax
  the comparison at the extremes of the octave).

This is `isExtremum` (`fasthessian.cxx:515`).

---

## 5. Sub-pixel / sub-scale refinement

An accepted extremum is refined by a Newton step in 4D space
`(x, y, z, s)` where `s` is scale:

```
X = −H⁻¹ ∇
```

where `∇` is the 4-vector of first derivatives (`deriv4D`) and `H` is the
4×4 Hessian (`hessian4D`), both estimated by finite differences. The inverse
is computed via SVD with small-eigenvalue truncation for numerical stability
(`interpolateStep`, `fasthessian.cxx:618`).

The point is accepted if all four offsets are within `[−1, 1]`, i.e. the
extremum is within one voxel and one filter-step of the detected location.
The final `Ipoint` is then:

```
x = (c + xX) * step
y = (r + xY) * step
z = (d + xZ) * step
scale = 0.1333 * (m->filter + xS * filterStep)
```

(`interpolateExtremum`, `fasthessian.cxx:577`).

---

## 6. Description

Two families of descriptors are available.

### 6.1 SURF3D descriptor (type 0, 48 dimensions)

For each keypoint the neighbourhood is divided into **8 sub-blocks**
(2 along each axis, block size = `radius`). Within each sub-block the Haar
wavelet responses in X, Y, Z and their absolute values are accumulated,
Gaussian-weighted by distance from the keypoint:

```
dx, dy, dz, |dx|, |dy|, |dz|
```

Each sub-block contributes 6 values; 8 sub-blocks → 48 values
(`Surf::getDescriptor`, `surf.cxx:63`). The vector is optionally normalized
to unit length.

### 6.2 Raw Haar descriptor (type 1, 3·8·radius³)

Every sample voxel contributes its raw Haar X/Y/Z responses
(`Surf::getRawDescriptor`, `surf.cxx:161`), producing `3·8·radius³` values.

### 6.3 Raw sub-volume descriptor (type 2, 8·radius³)

A `2·radius` cube around the keypoint is cropped from the resized image and
resampled to `2·radius` voxels per axis; the raw voxel values (8·radius³)
become the descriptor (`vtk3DSURF::ThreadedSubVolumes`,
`vtk3DSURF.cxx:274`).

---

## 7. Matching (match3d)

Matching is described separately in [SurfMatch](surfmatch.md). In short:
descriptors are matched by Euclidean distance with Laplacian-sign and
scale-ratio gating plus a best/second-best ratio test, then a RANSAC loop
estimates a similarity transform (SIM3).
