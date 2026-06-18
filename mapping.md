# Full Solution-Space Mapping: par ↔ poni

## 1. Notation

| Symbol | pyFAI name | ImageD11 name | Unit |
|---|---|---|---|
| L | `dist` | — | m (orthogonal distance) |
| Δ | — | `distance` | m (along-beam distance) |
| p₁, p₂ | `poni1`, `poni2` | — | m (PONI coordinates) |
| y_c, z_c | — | `y_center`, `z_center` | px (beam center) |
| s_h, s_v | `pixel2`, `pixel1` | `y_size`, `z_size` | m/px |
| θ₁, θ₂, θ₃ | `rot1`, `rot2`, `rot3` | — | rad (pyFAI rotations) |
| θx, θy, θz | — | `tilt_x`, `tilt_y`, `tilt_z` | rad (ImageD11 tilts) |
| O | orientation (1–4) | `o11,o12,o21,o22` | — (flip matrix) |

## 2. Coordinate Systems

**pyFAI lab frame** (left-handed): axis 1 = y_up (slow), axis 2 = x_starboard (fast), axis 3 = z_downstream.

**ImageD11 lab frame** (right-handed): axis X = downstream, axis Y = port, axis Z = up (slow).

Transform between frames: `t_ID11 = G · t_pyFAI` where

```
G = [[0, 0, 1],
     [0,-1, 0],
     [1, 0, 0]]
```
`G² = I`, so the transform is self-inverse.

## 3. Rotation Conventions

**pyFAI rotation** (intrinsic ZYX, left-handed R1/R2, right-handed R3):
```
R_pyFAI(θ₁,θ₂,θ₃) = Rz(θ₃) · Ry(-θ₂) · Rx(-θ₁)
```

**ImageD11 rotation** (all right-handed):
```
R_ID11(θx,θy,θz) = Rx(θx) · Ry(θy) · Rz(θz)
```

Applying the frame transform to pyFAI's rotation gives the **standard (uncompensated)** tilt-to-rot mapping:

```
θx = θ₃      (tilt_x = rot3)
θy = θ₂      (tilt_y = rot2)
θz = -θ₁     (tilt_z = -rot1)
```

These relations hold for the tilts **before** rotation compensation.
After compensation (below) the angle values differ.

## 4. The Affine Pipeline Equation

### pyFAI pipeline (per orientation O)

1. **Pixel reordering** (pre-rotation): `C(O) = diag(c1, c2)` where
   - O=3: c1=+1, c2=+1  (native, no flip)
   - O=2: c1=-1, c2=+1  (d1/slow flipped)
   - O=4: c1=+1, c2=-1  (d2/fast flipped)
   - O=1: c1=-1, c2=-1  (both flipped)

2. **Rotation**: `R = Rz(rot3)·Ry(-rot2)·Rx(-rot1)`

3. **Post-rotation sign flips**: `S(O) = diag(s1, s2, 1)` where
   - O=3: s1=+1, s2=+1  (none)
   - O=2: s1=-1, s2=+1  (t1 flipped, slow axis)
   - O=4: s1=+1, s2=-1  (t2 flipped, fast axis)
   - O=1: s1=-1, s2=-1  (both flipped)

Full pipeline: `t_pyFAI = S(O) · R · C(O) · d + translation`

### ImageD11 pipeline (per flip F)

Flip matrix extended to 3D: `Z(F) = diag(o11, -o22, 1)`.

`t_ID11_py_frame = R_tilt · Z(F) · d + translation`  (in pyFAI lab frame)

### Equating the linear parts

For the 2θ values to match, the full affine transforms must be equivalent:

```
S(O) · R_comp · C(O) = M · R_tilt · Z(F)
```

where M is one of the four mirror matrices (see §5).

## 5. Per-Orientation Mirror Matrices

To keep orthogonal distance positive and preserve 2θ/azimuth, a mirror
matrix M is applied to the right-hand side. The mirror relaxes strict xyz
coordinate matching in the ID11 frame but preserves 2θ and azimuth:

| Mirror | Matrix M | Flips in ID11 frame |
|--------|----------|-------------------|
| M3 (I) | diag( 1,  1, 1) | none (identity) |
| M2     | diag(-1,  1, 1) | axis-0 (maps to Z) |
| M4     | diag( 1, -1, 1) | axis-1 (maps to Y) |
| M1     | diag(-1, -1, 1) | both axes |

Each M is self-inverse (`M² = I`).

## 6. Solving for the Compensated Rotation

For a given (flip F, orientation O, mirror M) triple, the equation

```
S(O) · R_comp · C(O) = M · R_tilt · Z(F)
```

is solved column-by-column for the first two columns of R_comp:

```
R_comp[:,0] = S(O) · M · R_tilt[:,0] · (o11 / c1)
R_comp[:,1] = S(O) · M · R_tilt[:,1] · (-o22 / c2)
```

The third column is the cross product: `R_comp[:,2] = R_comp[:,0] × R_comp[:,1]`.

All scaling factors `o11/c1` and `-o22/c2` have magnitude 1, so the
columns remain orthonormal. The determinant is always +1 (since the
cross product of two orthonormal vectors has magnitude 1 and is orthogonal
to both).  No sign correction is needed.

ZYX Euler angles are then extracted from R_comp via scipy's
`Rotation.from_matrix`.  Each proper rotation matrix has exactly two
ZYX Euler representations (the β-solution pair), differing by a
cyclic ±180° shift.

## 7. The Full 16 × 4 Solution Enumeration

For a given ImageD11 `.par` dictionary (with a specific flip F), the
solution finder enumerates:

- **4 trial pyFAI orientations** O ∈ {1, 2, 3, 4}
- **4 mirror matrices** M ∈ {M1, M2, M3=I, M4}
- **2 ZYX Euler representations** per rotation matrix

→ 4 × 4 × 2 = **32 distinct solutions** per flip.

Each solution has:
- A pyFAI `.poni` dict with orientation O and rotation (rot1, rot2, rot3)
- An associated mirror M that defines the azimuth relationship
- The same 2θ values as the original ImageD11 geometry for all pixels

Half of the 32 solutions have positive orthogonal distance; half negative.

Across all 4 ImageD11 flip matrices, the 128 raw solutions deduplicate
to the same 32 unique (rot, dist) tuples — the solution set is
**flip-independent**.

## 8. Azimuth (chi ↔ eta) Mapping

The relationship between pyFAI chi (χ) and ImageD11 eta (η) depends
**only on the mirror matrix M**, not on the orientation or flip:

```
sin(χ) = M[0,0] · cos(η)
cos(χ) = M[1,1] · sin(η)
```

| Mirror   | sin(χ)      | cos(χ)      | χ = |
|----------|-------------|-------------|-----|
| M3 (I)   | +cos(η)     | +sin(η)     | 90° − η |
| M2       | −cos(η)     | +sin(η)     | η − 90° |
| M4       | +cos(η)     | −sin(η)     | η + 90° |
| M1       | −cos(η)     | −sin(η)     | 270° − η |

When the mirror matrix matches the trial orientation (mirror_orient ==
trial_orient), the combined S·M = I and χ = 90° − η exactly
(chi_eta_exact = True).  For other mirror choices, the azimuth
mapping follows the table above.

Derivation: from t_pyFAI = M · t_id (the mirror-transformed ID11
coordinates in the pyFAI lab frame) and the definitions
χ = atan2(t_pyFAI[0], t_pyFAI[1]), η = atan2(t_id[1], t_id[0]).

## 9. Beam Centre / PONI Formulas

With tilts, the PONI coordinates must account for the 0.5 pixel offset
and the orientation-specific pixel reordering:

**Forward (par → poni):**

```
dist = Δ · cos(θ₁) · cos(θ₂)

poni1 = -Δ · sin(θ₂) + pv · (beam_z_px + 0.5)
poni2 =  Δ · cos(θ₂) · sin(θ₁) + ph · (beam_y_px + 0.5)
```

where `beam_z_px = zc` for O ∈ {3, 4}, and `beam_z_px = max_d1 − zc` for O ∈ {2, 1}.
Similarly `beam_y_px = yc` for O ∈ {3, 2}, and `beam_y_px = max_d2 − yc` for O ∈ {4, 1}.

**Reverse (poni → par):**

```
Δ = L / (cos(θ₁) · cos(θ₂))

zc = (poni1 + L · tan(θ₂)/cos(θ₁)) / pv − 0.5    (O ∈ {3, 4})
zc = max_d1 + 0.5 − (poni1 + L · tan(θ₂)/cos(θ₁)) / pv  (O ∈ {2, 1})

yc = (poni2 − L · tan(θ₁)) / ph − 0.5             (O ∈ {3, 2})
yc = max_d2 + 0.5 − (poni2 − L · tan(θ₁)) / ph    (O ∈ {4, 1})
```

The `max_d1` / `max_d2` are `detector_shape[0]−1` / `detector_shape[1]−1`,
matching pyFAI's C-order convention `(slow, fast)`.

## 10. Backscattering Geometry

When the detector is **upstream** of the sample (ImageD11 distance < 0,
e.g. −0.15 m), the standard tilt mapping produces negative pyFAI distance.
The solution finder's `include_backscattering=True` flag explores π-offset
seed rotations (±π on rot1 or rot2) to discover representations where
pyFAI distance is positive and rot1≈π or rot2≈π.  With zero tilts this
gives exact 2θ matching; with non-zero tilts the match is approximate
(backscattering with significant tilts is physically unrealistic since
the central hole would need to accommodate the beam-spot movement).

## 11. Round-Trip Consistency

Each solution carries metadata in the poni dict:

| Key | Meaning |
|-----|---------|
| `orientation` | pyFAI orientation O (1–4) |
| `_mirror_used` | whether a non-identity mirror was used |
| `_mirror_orient` | which orientation's mirror (1–4) |
| `_forward_o11`, `_forward_o22` | original ImageD11 flip values |

This metadata is persisted in the `Detector_config` JSON block of `.poni`
files, enabling exact round-trip conversion via `poni_to_par()`.

## 12. Implementation

All conversion logic lives in `par_to_poni.py` with the following core
functions:

| Function | Purpose |
|----------|---------|
| `find_all_poni_solutions()` | Enumerate all 32 valid poni solutions for a par |
| `par_to_poni()` | Default conversion (prefers positive distance, canonical pairing) |
| `poni_to_par()` | Reverse conversion (reads mirror metadata from poni) |
| `_compute_compensated_rotation()` | Solve S·R·C = M·R_tilt·Z for R |
| `_compute_id11_from_pyfai()` | Reverse: recover R_tilt from R |
| `_azimuth_factors()` | Return (sf, cf) for a given mirror |
| `_find_all_rot_equivs()` | Find all ZYX Euler representations of a rotation matrix |
| `chi_to_eta()`, `eta_to_chi()` | Public azimuth conversion (canonical orientation) |

## 13. Test Coverage

37 tests covering 482 subtest variations (test_conversion.py):

- Round-trip identity (par→poni→par, poni→par→poni) for all 4 orientations
- 2θ matching ≤ 10⁻⁷ rad for all 32 solutions per flip
- Azimuth matching using correct mirror-dependent factors
- Full xyz lab-coordinate matching (non-square 128×200 detector)
- File I/O with correct metdata round-tripping
- Edge cases: zero tilts, zero pixel sizes, edge beam positions, tilts up to ±π/4
- Backscattering: positive distance, π rotations, all 4 orientations
- Default API unchanged from pre-refactor
- pyFAI rotation matrix validation against actual pyFAI output
- Integration test: write .poni → pyFAI.load() → integrate1d()

## 14. Limitations

- Transpose flips (o12, o21 ≠ 0) are not supported. The 4×4 affine
  analysis shows exact solutions exist, but the mapping to pyFAI's
  internal pixel-reordering requires understanding pyFAI's behaviour
  for transpose orientations.
- Spatial distortion is not handled.
- Backscattering with significant tilts does not preserve 2θ exactly
  (physically unrealistic anyway).
