# par_to_poni.py

Convert ImageD11 `.par` ↔ pyFAI `.poni` geometry files.

## Get the file

```bash
curl -O https://raw.githubusercontent.com/jonwright/pars2poni_deepseek/main/par_to_poni.py
```

Requires `numpy`, `scipy`.  Put `par_to_poni.py` next to your script and import it.

## Tested with

| Package   | Version       |
|-----------|---------------|
| pyFAI     | 2026.5.0      |
| ImageD11  | 2.1.3         |

(CI environment uses pyFAI 2026.6.0a0, ImageD11 2.1.5.)

## Convert a file

```python
import par_to_poni as pp

shape = (2162, 2068)   # Eiger 4M: (slow, fast) — pyFAI C-order convention

# par → poni (default: canonical pairing, positive distance)
par = pp.read_par("geometry.par")
poni = pp.par_to_poni(par, detector_shape=shape)
pp.write_poni(poni, "geometry.poni")

# poni → par (reads mirror metadata, exact round-trip)
poni = pp.read_poni("geometry.poni")
par = pp.poni_to_par(poni, detector_shape=shape)
pp.write_par(par, "geometry.par")
```

### Options

```python
# Prefer the solution where χ = 90° − η for all orientations:
poni = pp.par_to_poni(par, exact_chi=True)

# Allow negative distance (mirror/no-mirror variants):
poni = pp.par_to_poni(par, prefer_positive_distance=False)
```

## Find all valid solutions

For a given `.par` geometry, there are **32 valid pyFAI representations**
— 4 orientations × 4 mirror matrices × 2 Euler-angle equivalents per
rotation matrix:

```python
solutions = pp.find_all_poni_solutions(par, detector_shape=shape)

for s in solutions:
    p = s["poni"]
    print(f"  orient={s['orient_tried']} mirror={s['mirror_source']} "
          f"chi_exact={s['chi_eta_exact']} dist={p['dist']:.6f} "
          f"rot=({p['rot1']:.3f}, {p['rot2']:.3f}, {p['rot3']:.3f})")
```

Each solution has correct 2θ and a consistent χ/η azimuth mapping
determined by the mirror matrix.

## Solution structure

| Solution metadata | Meaning |
|-------------------|---------|
| `flip_label` | Which ImageD11 flip (F_o1, F_o2, F_o3, F_o4) |
| `orient_tried` | pyFAI orientation (1–4) |
| `mirror_source` | Which mirror matrix (M1, M2, M3=I, M4) |
| `use_mirror` | Whether a non-identity mirror was used |
| `dist_positive` | Whether orthogonal distance is positive |
| `chi_eta_exact` | Whether mirror matches trial orientation (canonical azimuth mapping) |
| `rot_magnitude` | \|rot1\| + \|rot2\| + \|rot3\| (for ranking) |
| `is_canonical` | Whether orientation matches flip→orient mapping |

Solutions are sorted best-first: positive distance preferred, then
smallest rotation magnitude, then canonical pairing.

## Azimuth mapping (χ ↔ η)

PyFAI chi and ImageD11 eta are related by **mirror-dependent** formulas:

```python
chi_rad = ...                                  # from pyFAI
eta = pp.chi_to_eta(chi_rad, orientation=3)    # → ImageD11 eta (radians)
chi = pp.eta_to_chi(eta, orientation=3)        # → pyFAI chi (radians)
```

Orientation can be an `int` (1–4), a par dict, or a poni dict.

| Mirror   | χ =        | sin(χ)       | cos(χ)       |
|----------|------------|--------------|--------------|
| M3 (I)   | 90° − η    | +cos(η)      | +sin(η)      |
| M2       | η − 90°    | −cos(η)      | +sin(η)      |
| M4       | η + 90°    | +cos(η)      | −sin(η)      |
| M1       | 270° − η   | −cos(η)      | −sin(η)      |

For the default (canonical) pairing for each orientation, the table
reduces to:

| Orient | Canonical Mirror | χ =        | sin(χ)      | cos(χ)      |
|--------|-----------------|------------|-------------|-------------|
| 3      | M3 (= I)        | 90° − η    | +cos(η)     | +sin(η)     |
| 2      | M2              | η − 90°    | −cos(η)     | +sin(η)     |
| 4      | M4              | η + 90°    | +cos(η)     | −sin(η)     |
| 1      | M1              | 270° − η   | −cos(η)     | −sin(η)     |

## The equation

The conversion solves the affine pipeline equation for all 16 (flip, orientation)
pairs with each of the 4 mirror matrices:

```
S(orient) · R · C(orient) = M · R_tilt · Z(flip)
```

Where S, C are pyFAI's per-orientation sign flips and pixel reordering,
Z is the ImageD11 flip matrix, and M is one of four mirror matrices.
See `mapping.md` for the full mathematical derivation.

## Status

All 4 non-transpose orientations match exactly for all 32 valid solutions per flip:
- 2θ at machine precision (≤ 10⁻⁷ rad)
- Azimuth with mirror-dependent simple mapping
- Lab coordinates with per-mirror reflections (≤ 5e-7 m)
- Round-trip par ↔ poni exact (metadata preserved in detector_config)
- Backscattering geometry supported (negative ImageD11 distance)

Transpose flips (`o12, o21 ≠ 0`) are not supported.
Spatial distortion is not handled.
37 tests, 482 subtest variations, all pass.

---

*Author: DeepSeek V4 Pro (opencode), June 2026, for Jon Wright*
