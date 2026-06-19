# par_to_poni.py

Convert ImageD11 `.par` → pyFAI `.poni` geometry files.

## Usage

```python
import par_to_poni as pp

shape = (2167, 2070)   # Eiger4M: (slow, fast) — pyFAI C-order convention

# par → poni (new pyFAI — orientation chosen per flip)
par = pp.read_par("geometry.par")
poni = pp.par_to_poni(par, detector_shape=shape)
pp.write_poni(poni, "geometry.poni")

# par → poni (old pyFAI — force orientation 3, no pixel flips)
poni3 = pp.par_to_poni(par, detector_shape=shape, force_orient3=True)
pp.write_poni(poni3, "geometry_orient3.poni")

# poni → par (works for both; recovers original flip automatically)
poni = pp.read_poni("geometry.poni")
par = pp.poni_to_par(poni, detector_shape=shape)
pp.write_par(par, "geometry.par")
```

`force_orient3` is for old pyFAI versions that predate flip support.
It always outputs orientation 3 with positive distance. When the
ImageD11 flip is non‑native, the rotation angles are compensated so
that 2θ and azimuth remain correct.  The original flip parameters
`(o11,o12,o21,o22)` are stored in the poni metadata so `poni_to_par` can
exactly recover them — so old‑pyFAI files round‑trip cleanly back
to ImageD11.  For native‑flip `(1,0,0,−1)` the two paths are identical.

## Conversion formula

Simple direct algebra (no solver, no iteration):

| pyFAI rot | ImageD11 tilt |
|-----------|---------------|
| `rot1` | `−tilt_z` |
| `rot2` | `tilt_y` |
| `rot3` | `tilt_x` |

The flip matrix is mapped to a pyFAI orientation number:

| Flip `(o11,o12,o21,o22)` | pyFAI orientation | Description |
|---------------------------|-------------------|-------------|
| `(1, 0, 0, -1)` | 3 | native, no flip |
| `(-1, 0, 0, 1)` | 1 | flip both axes |
| `(-1, 0, 0, -1)` | 2 | flip slow axis |
| `(1, 0, 0, 1)` | 4 | flip fast axis |

All 8 flips (the 4 above + the 4 transpose flips) also work with
`force_orient3=True`, mapping to orientation 3 with rotation
compensation (see [Classic mode](#classic-mode-force_orient3true) below).

Then PONI coordinates follow from orientation-dependent formulas
that account for the 0.5-pixel offset and pixel reordering.

## Azimuth mapping (χ ↔ η)

| Orient | χ = | sin(χ) | cos(χ) |
|--------|-----|--------|--------|
| 3 | 90° − η | +cos(η) | +sin(η) |
| 2 | η − 90° | −cos(η) | +sin(η) |
| 4 | η + 90° | +cos(η) | −sin(η) |
| 1 | 270° − η | −cos(η) | −sin(η) |

```python
eta = pp.chi_to_eta(chi_rad, orientation=3)
chi = pp.eta_to_chi(eta_rad, orientation=3)

# Orientation can come from a par or poni dict:
eta = pp.chi_to_eta(chi_rad, par)
chi = pp.eta_to_chi(chi_rad, poni)
```

When using `chi_to_eta` / `eta_to_chi` with a par or poni dict, the
orientation is extracted automatically — modern (flip‑matched) or
classic (always orient 3) as appropriate.

## Classic mode (`force_orient3=True`)

For old pyFAI versions that only accept orientation 3 (no pixel flips).
Always outputs orientation 3 with positive distance.  Non‑native flips
get rotation compensation via the unified formula

```
R_comp = M · R_tilt · Z
   M = diag(−1, 1, 1)  if det(Z) < 0,  else I
   Z = [[o11, o12, 0], [−o21, −o22, 0], [0, 0, 1]]
```

so that 2θ and azimuth stay correct.  The original `(o11,o12,o21,o22)`
are stored in the poni metadata for exact round‑trip back to ImageD11.

### Azimuth

In classic mode the azimuth mapping is **always** `χ = 90° − η`
(orientation 3), regardless of the original ImageD11 flip:

| ImageD11 flip | Modern χ = | Classic χ = |
|---------------|------------|-------------|
| `(1,0,0,−1)` → orient 3 | 90° − η | 90° − η |
| `(−1,0,0,1)` → orient 1 | 270° − η | 90° − η |
| `(−1,0,0,−1)` → orient 2 | η − 90° | 90° − η |
| `(1,0,0,1)` → orient 4 | η + 90° | 90° − η |

For the native flip `(1,0,0,−1)` the two modes are identical.
For all other flips, modern and classic χ differ — both are
internally consistent, and `poni_to_par` reverses each correctly.

### Transpose flips (axis swap)

Flips with off‑diagonal entries `(o12≠0 or o21≠0)` represent 90° detector
rotations.  These require `force_orient3=True` — pyFAI has no native
orientation for axis swaps.  All 4 transpose flips map to orientation 3
with the same `M·R_tilt·Z` compensation used for all non‑native flips:

| Flip `(o11,o12,o21,o22)` | Orientation | Notes |
|---|---|---|
| `(0, 1, -1, 0)` | 3 | 90° rotation, det(Z)=+1 |
| `(0, -1, 1, 0)` | 3 | −90° rotation, det(Z)=+1 |
| `(0, 1, 1, 0)` | 3 | 90° rot + reflection, det(Z)=−1 |
| `(0, -1, -1, 0)` | 3 | −90° rot + reflection, det(Z)=−1 |

All 4 have correct round‑trip and machine‑precision 2θ (same as
non‑transpose).  Usage example:

```python
par["o11"], par["o12"], par["o21"], par["o22"] = 0, 1, 1, 0  # axis swap
poni = pp.par_to_poni(par, detector_shape=shape, force_orient3=True)
```

## Tested with

| Package | Version |
|---------|---------|
| pyFAI | 2026.5.0 |
| ImageD11 | 2.1.3 |

10 test classes, 35 tests, 148 subtests.
All 2θ ≤ 1e-7 rad for all 8 flips. Round‑trip ≤ 5e-13 m.

---

*Author: DeepSeek V4 Pro (opencode), June 2026, for Jon Wright*
