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
rotation is compensated via the full 3×3 Z matrix:

```
Z = [[o11,  o12,  0],
     [−o21, −o22, 0],
     [ 0,    0,   1]]
```

When `det(Z) = −1`, pre‑multiplying by mirror M2 = `diag(−1, 1, 1)`
gives a proper rotation:

```
R_comp = M · R_tilt · Z      where M = M2 if det(Z) < 0 else I
```

Euler angles are extracted from `R_comp` and the positive‑distance
equivalent is found.  The reverse path uses `R_tilt = M·R_comp·Z⁻¹`
with `Z⁻¹ = Zᵀ` (Z is orthogonal, M is self‑inverse).

This unified formula handles all 8 flips (4 non‑transpose + 4 transpose).
The original `(o11,o12,o21,o22)` are stored in the poni metadata for
exact round‑trip.

For modern pyFAI, the non‑transpose flips use the direct mapping
`rot = (−tz, ty, tx)` with the orientation from the table above —
no mirror compensation is needed.  Transpose flips require
`force_orient3=True`.  `par_to_poni` without `force_orient3` raises
`ValueError` for transpose.

## Why this is the full solution

The affine pipeline equation `S·R·C = M·R_tilt·Z` has at most
one working (flip, orientation) pairing per flip.  For the modern
path (non‑transpose flips with native orientation), M = I and
the rotation R equals `R_tilt` directly — no solver needed.
For the classic path (`force_orient3`), M is chosen from
{M2, I} based on `det(Z)`, and the compensated rotation is
`R_comp = M·R_tilt·Z`.

For the complete derivation and historical 32‑solution enumeration,
see `detailed_analysis/mapping.md`.
