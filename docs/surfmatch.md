# SurfMatch (`match3d`)

> **Last updated:** 2026-10-07

`match3d` aligns two sets of keypoints by finding correspondences and
estimating a **similarity transform** (SIM3: rotation, translation, uniform
scale) between them. It is implemented by `MatchPoint` (`MatchPoint.h/.cxx`)
with the CLI in `MainMatch.cxx`.

---

## 1. Input

Two JSON keypoint files (from `surf3d -json 1`), each parsed into an `IpVec`
(`MatchPoint::Parse`, `MatchPoint.cxx:16`). Fields read per point:
`laplacian`, `response`, `scale`, `x`, `y`, `z`, and the `descriptor` array.

---

## 2. Matching (`computeMatches`, `MatchPoint.cxx:62`)

For each point in set 1, find its best match in set 2, subject to:

- **Laplacian gate** (`line 111`): the Laplacian sign must match. This is the
  first cheap rejection test.
- **Scale gate** (`line 115`): the scale ratio
  `max(s1/s2, s2/s1) <= MatchingScale` (default 1.5).
- **Descriptor distance** (`line 120`): Euclidean distance between
  descriptors.
- **Best/second-best ratio** (`line 134`): a match is kept only if
  `d1/d2 < MatchingDist2Second` (default 0.98) **and** `d1 < MatchingDist`
  (default 0.95). This is the standard "ratio test" of Lowe that rejects
  ambiguous matches.

Optionally, matching is restricted to a user-supplied bounding box
(`useBBoxin`, via `-bb`).

The result is `IpPairVec matches` of `(i, bestJ)` index pairs.

---

## 3. Transform estimation (`computeTransform`, `MatchPoint.cxx:157`)

A **RANSAC** loop robustly estimates the SIM3:

- Repeatedly (up to 8000 iterations, max 10000 failures):
  1. Randomly sample **3** matches (`N = 3`).
  2. Estimate a SIM3 from them (`getRT`, Umeyama-style least-squares).
  3. Count **inliers**: matches whose transformed point lies within
     `RansacDist` (default 50) of its partner, using the inverse transform.
  4. Track the hypothesis with the most inliers.
- Requires `matches.size() > 3` and final `maxinlier >= RansacMinInliers`
  (default 2), else the run is marked as `echec = true` (failure).

### `getRT` (`MatchPoint.cxx:271`)

Estimates a SIM3 from a set of point pairs:

- Computes centroids of both sets.
- Computes covariance matrices `Ca = AᵀA`, `Cb = BᵀB`.
- Estimates scale as `s = sqrt((EigCa·EigCb)/(EigCa·EigCa))` from the SVD
  eigenvalues.
- **Rotation:** currently hard-coded to `Identity`
  (`MatchPoint.cxx:344`). The SVD-based rotation estimation (Umeyama) is
  present but commented out (`MatchPoint.cxx:323-343`).
- Translation = `centroidA − R·centroidB`.

> **Known issue:** because rotation is fixed to Identity, the estimated
> transform currently ignores rotational alignment. See
> [TODO](../TODO.md) (H6).

---

## 4. Output (`WriteTransform`, `MatchPoint.cxx:353`)

Writes `transform.json` (hard-coded name in the CWD):

```json
{
  "translation": [tx, ty, tz],
  "rotation": [[r00,r01,r02],[r10,r11,r12],[r20,r21,r22]],
  "scale": 1.0,
  "inliers": 42,
  "fail": false,
  "bboxA": {"min": [...], "max": [...]},
  "bboxB": {"min": [...], "max": [...]},
  "nbPointInA": 100,
  "nbPointInB": 100,
  "allInliers": [[i0,j0],[i1,j1], ...]   // only with -i
}
```

- `bboxA`/`bboxB` are the bounding boxes of the inliers in each set (empty
  unless bounding boxes were requested).
- `allInliers` lists the matched index pairs (only when `-i` is given).

---

## Tuning parameters (defaults)

| Parameter | Default | Meaning |
|-----------|---------|---------|
| `-Rd` RansacDist | 50 | inlier distance threshold |
| `-Ri` RansacMinInliers | 2 | minimum inliers to accept |
| `-Md` MatchingDist | 0.95 | max descriptor distance |
| `-Md2` MatchingDist2Second | 0.98 | best/second-best ratio |
| `-Ms` MatchingScale | 1.5 | scale-ratio gate |
