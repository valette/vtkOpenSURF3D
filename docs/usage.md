# Usage

> **Last updated:** 2026-10-07

This page documents the command-line interfaces of `surf3d` and `match3d`.

---

## `surf3d`

```
surf3d file [options]
```

`file` is a 3D image. Supported formats: **NIFTI**, **MHD** (via VTK readers).

### Options

| Option | Argument | Description | Default |
|--------|----------|-------------|---------|
| `-bin` | `0/1` | write points as a `.bin` file | `0` |
| `-cmin` | value | clamp values lower than `value` | – |
| `-cmax` | value | clamp values larger than `value` | – |
| `-csv` | `0/1` | write points as a `.csv` file | `0` |
| `-csvgz` | `0/1` | write points as a `.csv.gz` file | `1` |
| `-json` | `0/1` | write points as a `.json` file | `0` |
| `-n` | number | maximum number of points | `-1` (all) |
| `-normalize` | `0/1` | normalize descriptors | `1` |
| `-o` | basename | output file base name | `points` |
| `-r` | radius | descriptor volume radius | `5` |
| `-s` | spacing | resample to isotropic spacing | `0` (off) |
| `-t` | threshold | detector threshold | `0` (CLI) / `0.0004` (class) |
| `-type` | `0/1/2` | descriptor type | `0` |
| `-d` | maxsize | max size for isotropic resampling | `0` |
| `-gz` | opts | gzip options (e.g. compression level) | – |
| `-m` | maskfile | mask image file | – |
| `-nt` | threads | number of threads for VTK filters | `-1` (auto) |
| `-p` | pointfile | reuse a pre-existing points file (CSV) | – |
| `-pad` | value | mirror padding in voxels | `0` |
| `-precision` | n | precision of csv.gz coefficients | `-1` (full) |

### Descriptor types (`-type`)

| Value | Description | Descriptor size |
|-------|-------------|-----------------|
| `0` | SURF3D descriptor (default) | 48 |
| `1` | sub-volume HAAR coefficients | `24 · radius³` |
| `2` | sub-volume raw voxels | `8 · radius³` |

### Examples

```bash
# Detect points, write default csv.gz + bounds json
./surf3d img.nii.gz

# Write JSON keypoints with a higher threshold
./surf3d img.nii.gz -t 0.01 -json 1 -csvgz 0

# Reuse existing points (descriptor-only pass), limit to 10 points
./surf3d img.nii.gz -p points.csv -n 10

# Resample to 2mm isotropic before detection
./surf3d img.nii.gz -s 2.0

# Mirror-pad by 5 voxels and clamp intensities
./surf3d img.nii.gz -pad 5 -cmin 0 -cmax 1000
```

### Exit codes

| Code | Meaning |
|------|---------|
| `0` | success |
| `1` | usage error / missing option value |
| `5` | failed to load input image (or mask file) |
| `6` | mask/image dimension or component mismatch |

---

## `match3d`

```
match3d file1.json file2.json [options]
```

Both files are **JSON keypoint files** produced by `surf3d -json 1` (the
`points` array).

### Options

| Option | Argument | Description | Default |
|--------|----------|-------------|---------|
| `-i` | (flag) | write the inlier index pairs (`allInliers`) | off |
| `-Rd` | value | RANSAC distance threshold | `50` |
| `-Ri` | value | RANSAC minimum inliers | `2` |
| `-Md` | value | descriptor matching distance | `0.95` |
| `-Md2` | value | second-best ratio | `0.98` |
| `-Ms` | value | scale-ratio gate | `1.5` |
| `-b` | `0/1` | compute bounding boxes | `0` |
| `-bb` | b64 | base64-encoded bounding box JSON | – |

### Output

Writes **`transform.json`** in the current working directory (hard-coded
filename) containing `translation`, `rotation`, `scale`, `inliers`, `fail`,
and optionally `allInliers` (with `-i`).

### Example

```bash
./match3d a.json b.json -i
```
