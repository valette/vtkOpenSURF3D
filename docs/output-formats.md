# Output Formats

> **Last updated:** 2026-10-07

`surf3d` writes detected keypoints in several formats. The enabled formats
are controlled by `-json`, `-csv`, `-csvgz`, `-bin`. The **base file name** is
set by `-o` (default `points`), and each enabled format appends its own
suffix.

In addition to the keypoint files, `surf3d` always writes a **bounds JSON**
file `<base>.json` describing the physical-space bounds, origin, spacing and
dimensions of the processed image. This is what `match3d` uses (via its
`-bb` option, base64-encoded) to constrain matching.

---

## Bounds file (`<base>.json`)

Written unconditionally by `surf3d.cxx:296`. Structure:

```json
{
  "bounds": [
    [xmin, xmax],
    [ymin, ymax],
    [zmin, zmax]
  ],
  "origin": [ox, oy, oz],
  "spacing": [sx, sy, sz],
  "dimensions": [dx, dy, dz]
}
```

All values are in **physical space** (voxel coordinates multiplied by the
image spacing, relative to the `Cast` image origin).

---

## JSON keypoints (`-json 1` → `<base>.json`)

Written by `vtk3DSURF::WritePoints` (`vtk3DSURF.cxx:349`). Structure:

```json
{
  "points": [
    {
      "x": 12.3, "y": 45.6, "z": 78.9,
      "laplacian": -1,
      "response": 0.0123,
      "scale": 4.56,
      "descriptor": [ ... ]
    }
  ],
  "origin": [ox, oy, oz],
  "spacing": [sx, sy, sz],
  "dimensions": [dx, dy, dz]
}
```

- `x/y/z` are in **physical** coordinates (voxel × spacing).
- `scale` is `point.scale * Sspacing` where
  `Sspacing = (sx·sy·sz)^(1/3)` (geometric mean of spacing).
  > **Known issue:** for the JSON writer `Sspacing` is derived from the `-s`
  > option member (`vtk3DSURF::Spacing`), which is `0.0` when `-s` is not
  > given, producing `scale = 0`. See [TODO](../TODO.md) (H2).
- `descriptor` is a flat array of floats (48 for `-type 0`).

---

## CSV (`-csv 1` → `<base>.csv`)

Written by `vtk3DSURF::WritePointsCSV` (`vtk3DSURF.cxx:417`).

Header: `x,y,z,scale,response,laplacian,descriptor`

One row per keypoint. The descriptor is appended as space-separated values
after the laplacian.

---

## CSV.GZ (`-csvgz 1` → `<base>.csv.gz`)

Written by `vtk3DSURF::WritePointsCSVGZ` (`vtk3DSURF.cxx:456`) using ZLIB.
Identical content to the CSV, gzip-compressed (default level; `-gz` option
and `-precision` control compression/formatting).

---

## Binary (`-bin 1` → `<base>.bin`)

Written by `vtk3DSURF::WritePointsBinary` (`vtk3DSURF.cxx:518`).

A raw binary dump:

```
int    count                      (number of keypoints)
for each keypoint:
  float x, y, z, scale, response
  int   laplacian
  float descriptor[48]
```

Note: the binary format currently always writes a **48-element** descriptor
regardless of `-type`.

---

## Point file input (`-p`)

The `-p` option tells `surf3d` to **reuse existing keypoints** instead of
running the detector. The input must be a **CSV** file matching the CSV
output format above. Each row's `x/y/z` (physical) and `scale` are converted
back to voxel units using the `Cast` image origin/spacing
(`vtk3DSURF::ReadIPoints`, `vtk3DSURF.cxx:35`), and descriptors are
(re)computed according to `-type`.
