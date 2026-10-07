# vtkOpenSURF3D — Developer Documentation

> **Last updated:** 2026-10-07
> **Source commit:** `24b550f` (docs describe this commit; update this
> reference whenever the documentation is refreshed to track code changes)

This directory contains the technical documentation for **vtkOpenSURF3D**, a
C++ implementation of 3D Speeded Up Robust Features (SURF) built on VTK.

## Contents

| Document | Description |
|----------|-------------|
| [Architecture](architecture.md) | High-level pipeline, data flow and threading model |
| [Algorithm](algorithm.md) | The 3D SURF detection/description algorithm in detail |
| [Modules](modules.md) | Per-file / per-class reference of the codebase |
| [Build](build.md) | Dependencies and building the two executables |
| [Usage](usage.md) | `surf3d` and `match3d` command-line reference |
| [Output formats](output-formats.md) | JSON / CSV / CSV.GZ / binary keypoint formats |
| [SurfMatch](surfmatch.md) | Feature matching, RANSAC and transform estimation |
| [Testing](testing.md) | The `test.sh` suite and how to run it |

There is also a companion [Code Audit / TODO](../TODO.md) at the repository
root describing known bugs, memory issues and recommended improvements,
organized by severity.

## Quick links for coding agents

- Repo root: `AGENTS.md` contains the essential commands, conventions and
  gotchas an agent needs to work on this codebase.
- The two executables are `surf3d` (detection + description) and `match3d`
  (matching + transform estimation, optional build).
- Detection core: `FastHessian` (`fasthessian.h/.cxx`) + `ResponseLayer`
  (`responselayer.h`) + integral image (`integral.h/.cxx`).
- Description: `Surf` (`surf.h/.cxx`) and `vtk3DSURF::ThreadedSubVolumes`.
- Orchestration / CLI / IO: `vtk3DSURF.h/.cxx`, `surf3d.cxx`.
- Image loading: `vtkRobustImageReader.h`.
- Matching: `MatchPoint.h/.cxx`, `MainMatch.cxx`.
