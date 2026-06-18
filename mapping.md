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

Each ImageD11 flip `Z = diag(o11, −o22)` matches exactly one pyFAI
orientation, which encodes the same pixel-axis flips via its C and S
matrices:

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

## Why this is the full solution

The affine pipeline equation `S·R·C = M·R_tilt·Z` has exactly one
correct (flip, orientation) pairing.  When the flip and orientation match
— with the corrected table above — the right-hand-side mirror M is
the identity and the rotation R equals the uncompensated `R_tilt`.
No solver, no mirror matrices, no rotation compensation is needed.

For the complete derivation and historical 32-solution enumeration,
see `detailed_analysis/mapping.md`.
