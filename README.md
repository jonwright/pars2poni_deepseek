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
`(o11,o22)` are stored in the poni metadata so `poni_to_par` can
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
```

## Details

See `par_to_poni.py` source and `detailed_analysis/mapping.md`
for the complete mathematical derivation and solver history.

## Tested with

| Package | Version |
|---------|---------|
| pyFAI | 2026.5.0 |
| ImageD11 | 2.1.3 |

(CI: pyFAI 2026.6.0a0, ImageD11 2.1.5)

9 test classes, ~70 subtest variations including 60° tilts.
All 2θ ≤ 1e-7 rad, round-trip ≤ 5e-13 m, azimuth sin/cos ≤ 1e-7.

---

*Author: DeepSeek V4 Pro (opencode), June 2026, for Jon Wright*
