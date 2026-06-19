# par ↔ poni — Key Formulas (simplified, Round 9)

## Coordinate systems

- **pyFAI**: axis1=y_up(slow), axis2=x_starboard(fast), axis3=z_downstream
- **ImageD11**: X=downstream, Y=port, Z=up(slow)
- **Frame transform**: `t_ID11 = G·t_pyFAI` with `G = [[0,0,1],[0,-1,0],[1,0,0]]`

## Rotations

- **pyFAI**: `R = Rz(rot3)·Ry(−rot2)·Rx(−rot1)`
- **ImageD11**: `R_ID11 = Rx(tilt_x)·Ry(tilt_y)·Rz(tilt_z)`

Applying the G-matrix frame transform gives the direct mapping:

```
rot1 = −tilt_z
rot2 =  tilt_y
rot3 =  tilt_x
```

…and the inverse:

```
tilt_x =  rot3
tilt_y =  rot2
tilt_z = −rot1
```

## Orientation mapping

### Modern pyFAI

Each non‑transpose ImageD11 flip `Z = diag(o11, −o22)` matches exactly
one pyFAI orientation:

| Flip | Z | pyFAI orient | Flips |
|------|---|-------------|-------|
| `(1,0,0,-1)` | `diag(1,1)` | 3 | none |
| `(-1,0,0,1)` | `diag(−1,−1)` | 1 | both |
| `(-1,0,0,-1)` | `diag(−1,1)` | 2 | slow |
| `(1,0,0,1)` | `diag(1,−1)` | 4 | fast |

## PONI coordinates

```
dist = Δ·cos(rot1)·cos(rot2)

beam_z = max_d1 − zc  if orient ∈ {2,1}  else  zc
beam_y = max_d2 − yc  if orient ∈ {4,1}  else  yc

poni1 = −Δ·sin(rot2) + pv·(beam_z + 0.5)
poni2 =  Δ·cos(rot2)·sin(rot1) + ph·(beam_y + 0.5)
```

…and inverse:

```
Δ = L / (cos(rot1)·cos(rot2))

zc = max_d1 + 0.5 − (poni1 + L·tan(rot2)/cos(rot1)) / pv   for orient ∈ {2,1}
     (poni1 + L·tan(rot2)/cos(rot1)) / pv − 0.5             for orient ∈ {3,4}

yc = max_d2 + 0.5 − (poni2 − L·tan(rot1)) / ph             for orient ∈ {4,1}
     (poni2 − L·tan(rot1)) / ph − 0.5                       for orient ∈ {2,3}
```

`max_d1 = detector_shape[0] − 1`, `max_d2 = detector_shape[1] − 1`.

## Azimuth

| Orient | χ = | sin(χ) | cos(χ) |
|--------|-----|--------|--------|
| 3 | 90°−η | +cos η | +sin η |
| 2 | η−90° | −cos η | +sin η |
| 4 | η+90° | +cos η | −sin η |
| 1 | 270°−η | −cos η | −sin η |

## force_orient3 — classic pyFAI (orientation 3 only)

Old pyFAI versions accept only orientation 3.  For any flip, the
rotation is compensated so that 2θ and azimuth stay correct while the
output orientation is always 3.

### The formula

For every ImageD11 flip, define the full 3×3 extension of the flip
matrix (embedding the 2×2 `[[o11,o12],[−o21,−o22]]` into a 3×3 with
the z‑axis untouched):

```
Z = [[o11,  o12,  0],
     [−o21, −o22, 0],
     [ 0,    0,   1]]
```

Z is orthogonal: its columns are orthonormal and `Z⁻¹ = Zᵀ`.

When `det(Z) = −1`, Z is a reflection, not a rotation.  Pre‑multiplying
by the mirror M2 = `diag(−1, 1, 1)` yields a proper rotation:

```
R_comp = M · R_tilt · Z      where M = M2 if det(Z) < 0 else I
```

`R_tilt` is the standard tilt rotation from
`rot = (−tilt_z, tilt_y, tilt_x)`.  Euler angles `(rot1,rot2,rot3)`
are extracted from `R_comp` and the positive‑distance equivalent
(`cos(rot1)·cos(rot2) > 0`) is found.  Since `R_comp` is a proper
rotation (det = +1) by construction, a positive‑distance equivalent
always exists.

### All 8 flips

| Flip `(o11,o12,o21,o22)` | Z diag | Z off‑diag | det(Z) | M | orient |
|---|---|---|---|---|---|
| `(1, 0, 0, −1)` (native) | `diag(1,1,1)` | — | +1 | `I` | 3 |
| `(−1, 0, 0, 1)` | `diag(−1,−1,1)` | — | +1 | `I` | 3 |
| `(−1, 0, 0, −1)` | `diag(−1,1,1)` | — | −1 | `M2` | 3 |
| `(1, 0, 0, 1)` | `diag(1,−1,1)` | — | −1 | `M2` | 3 |
| `(0, 1, −1, 0)` | `diag(0,0,1)` | 90° rot | +1 | `I` | 3 |
| `(0, −1, 1, 0)` | `diag(0,0,1)` | −90° rot | +1 | `I` | 3 |
| `(0, 1, 1, 0)` | `diag(0,0,1)` | 90° rot + refl | −1 | `M2` | 3 |
| `(0, −1, −1, 0)` | `diag(0,0,1)` | −90° rot + refl | −1 | `M2` | 3 |

All 8 produce `orientation = 3` in the output.  The original
`(o11,o12,o21,o22)` are stored in the poni metadata for exact
round‑trip.

### PONI coordinates (orientation 3)

Since classic mode always outputs orientation 3, the PONI formulas
simplify to the native case:

```
dist  = Δ·cos(rot1)·cos(rot2)                         (with compensated rot1,rot2)
poni1 = −Δ·sin(rot2) + pv·(zc + 0.5)
poni2 =  Δ·cos(rot2)·sin(rot1) + ph·(yc + 0.5)
```

where `Δ = par["distance"]`, `zc = par["z_center"]`,
`yc = par["y_center"]`, and `pv=zs`, ` ph=ys`.

### Reverse path

```
R_tilt = M · R_comp · Zᵀ     (M is self‑inverse, Z⁻¹ = Zᵀ)
```

Euler angles are extracted from `R_tilt` giving back the original
ImageD11 tilts: `tilt_x = rot3`, `tilt_y = rot2`, `tilt_z = −rot1`.

### Modern vs classic

For modern pyFAI, the non‑transpose flips use the direct mapping
`rot = (−tz, ty, tx)` with the orientation from the table above —
no mirror compensation is needed.  Transpose flips require
`force_orient3=True`.  `par_to_poni` without `force_orient3` raises
`ValueError` for transpose.

## Why this is the full solution

For the modern path (non‑transpose flips with native orientation),
the rotation R equals `R_tilt` directly — no solver needed because
the `_FLIP_TO_ORIENTATION` table correctly matches each flip to the
pyFAI orientation that applies the same pixel transform.

For the classic path (`force_orient3`), the compensated rotation
`R_comp = M·R_tilt·Z` maps any flip to orientation 3.  Z embeds
the ImageD11 flip into a 3×3 orthogonal matrix; M is the identity
when Z is a proper rotation (det = +1) or M2 = diag(−1, 1, 1) when
det(Z) = −1, ensuring `R_comp` is always a proper rotation with a
positive‑distance Euler representation.

The full derivation and historical 32‑solution enumeration are in
`detailed_analysis/mapping.md`.
