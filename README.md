# par_to_poni.py

Convert ImageD11 `.par` → pyFAI `.poni` geometry files.

## Usage

```python
import par_to_poni as pp

shape = (2167, 2070)   # Eiger4M: (slow, fast) — pyFAI C-order convention

# par → poni
par = pp.read_par("geometry.par")
poni = pp.par_to_poni(par, detector_shape=shape)
pp.write_poni(poni, "geometry.poni")

# poni → par
poni = pp.read_poni("geometry.poni")
par = pp.poni_to_par(poni, detector_shape=shape)
pp.write_par(par, "geometry.par")
```

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

7 test classes, ~60 subtest variations including 60° tilts.
All 2θ ≤ 1e-7 rad, round-trip ≤ 5e-13 m, azimuth sin/cos ≤ 1e-7.

---

*Author: DeepSeek V4 Pro (opencode), June 2026, for Jon Wright*
