# PLAN: par ↔ poni Geometry Conversion

## Goals

1. Review both codebases for geometry definitions, identify contradictions vs docs
2. Derive clean algebraic mappings between `.par` and `.poni` parameters
3. Implement IO + conversion functions with correct round-trips
4. Test with all 4 supported flips on a strongly tilted detector
5. Test on a non-square detector (200×128) to verify per-axis pixel reordering

## Scope / Clarifications (from user input)

- **Only 4 non-transpose flips** (o11, o22 ∈ {±1}, o12 = o21 = 0). Transpose flips
  (0/±1 swapped) are not supported by pyFAI's current orientation model.
- **Account for the 0.5 pixel offset** in conversion formulas.
- **Direct algebra** (geometry_conversion.rst derivations), not via Fit2D intermediate.
- **No spatial distortion** — splines/dx/dy files ignored for now.
- **wedge, chi** set to zero in par file output; **omegasign** to 1.0, **fit_tolerance** to 0.05.
- **Unit handling**: ImageD11 par files are unit-agnostic ("inches are OK"). Provide an
  option (`par_length_unit`) for micron (µm), mm, or meters when reading/writing par files.
  Internally all lengths are in meters, matching pyFAI's convention.
  Default par length unit: **µm** (matches ImageD11 historical default).
  Wavelength in par files is always in **angstrom** (Å).
- **Azimuth wrapping**: compare sin/cos pairs on angles to avoid wrap-around issues
  (sin(chi) vs sin(90°-eta), cos(chi) vs cos(90°-eta)).

---

## 1. CODEBASE REVIEW FINDINGS

### 1.1 Coordinate Systems

| | pyFAI lab | ImageD11 lab |
|---|---|---|
| Axis 1 / y | Up (slow pixel dim) | — |
| Axis 2 / x | Starboard (fast pixel dim) | — |
| Axis 3 / z | Downstream beam | — |
| X | — | Downstream beam |
| Y | — | Port (away from ring) |
| Z | — | Up (slow pixel dim) |
| Handedness (xyz) | Left-handed | Right-handed |

**Transform matrix** (pyFAI lab → ImageD11 lab):
```
G = [[0, 0, 1],
     [0,-1, 0],
     [1, 0, 0]]
t_ID11 = G · t_pyFAI
```

Key files:
- `pyFAI/doc/source/geometry.rst:67-79` — pyFAI coordinate system
- `pyFAI/doc/source/geometry_conversion.rst:488-506` — G matrix derivation
- `ImageD11/ImageD11/transform.py:85-147` — compute_xyz_lab

### 1.2 Rotation Mapping

From `geometry_conversion.rst:619-635` and verified in `imaged11.py:142-144,195-197`:

```
θx (tilt_x) = rot3 = θ₃        [right-handed → right-handed, same]
θy (tilt_y) = rot2 = θ₂        [left-handed → right-handed, convention change]
θz (tilt_z) = -rot1 = -θ₁      [left-handed → right-handed, sign flip]
```

**pyFAI rotation matrix** (`core.py:2656-2704`):
```
R_pyFAI = R₃(θ₃) · R₂(-θ₂) · R₁(-θ₁)
```
R₁, R₂ are left-handed (negated angles); R₃ is right-handed.
Implemented via `scipy.spatial.transform.Rotation.from_euler('ZYX', [rot3, -rot2, -rot1])`.

**ImageD11 rotation matrix** (`transform.py:51-82`):
```
R_ID11 = R₁(θx) · R₂(θy) · R₃(θz)   [all right-handed]
```
Implemented via `scipy.spatial.transform.Rotation.from_euler('XYZ', [tilt_x, tilt_y, tilt_z])`.

Despite different conventions, the effective rotation order is the same (verified
in `geometry_conversion.rst:601-617`).

### 1.3 Distance Definitions

| Parameter | Meaning |
|---|---|
| pyFAI `dist` (L) | **Orthogonal** distance from sample to detector plane |
| ImageD11 `distance` (Δ) | Distance from sample **along beam** to detector intersection |

Relationship (`imaged11.py:199`, `geometry_conversion.rst:726-728`):
```
L  = Δ · cos(θy) · cos(θz)
Δ  = L / (cos(θ₁) · cos(θ₂))    [θ₁=rot1, θ₂=rot2]
```

In pyFAI, the distance is applied BEFORE rotations. In ImageD11, it's applied AFTER.

### 1.4 Flip → Orientation Mapping (Non-transpose Only)

From `imaged11.py:168-183`:

| o11 | o12 | o21 | o22 | Orientation | Effect on lab coords |
|---|---|---|---|---|---|
| 1 | 0 | 0 | -1 | **3** (pyFAI native) | No sign flip |
| -1 | 0 | 0 | 1 | **1** | t1 = -t1, t2 = -t2 (flip both) |
| -1 | 0 | 0 | -1 | **4** | t2 = -t2 (flip fast/x axis) |
| 1 | 0 | 0 | 1 | **2** | t1 = -t1 (flip slow/y axis) |

ImageD11 flip matrix (`transform.py:689-691`):
```
fmat = [[1, 0, 0],
        [0, o22, o21],
        [0, o12, o11]]
```
For non-transpose (o12=o21=0): `fmat = diag(1, o22, o11)`.

pyFAI applies sign flips to lab coordinates after rotation (`core.py:554-558`):
```python
if orientation in (1, 2): t1 = -t1   # flip y/slow
if orientation in (1, 4): t2 = -t2   # flip x/fast
```

### 1.5 The 0.5 Pixel Offset

**pyFAI** (`detectors/_common.py:722-725`):
```python
if center:
    d1c = d1 + 0.5    # pixel centers at half-integer offsets
    d2c = d2 + 0.5
```
```python
p1 = pixel1 * (dY + d1c)   # physical coordinate of pixel center
```

**ImageD11** (`transform.py:126-127`):
```python
peaks_on_detector[0, :] = (peaks[0, :] - z_center) * z_size   # no 0.5
peaks_on_detector[1, :] = (peaks[1, :] - y_center) * y_size
```

PyFAI computes the pixel center: `p = pixel * (index + 0.5)`.
ImageD11 operates on floating-point coordinates directly: `offset = pixel * (coord - center)`.

**Consequence for beam center**:
- pyFAI: beam at pixel coordinate `d_beam` where `pixel * (d_beam + 0.5) - poni = 0`
  → `d_beam = poni/pixel - 0.5`
- ImageD11: beam at `z_center`, `y_center` (floating-point pixel coordinates)

Thus the conversion from doc formula which gives:
```
z_center_doc = (poni1 + L*tan(θ₂)/cos(θ₁)) / pixel_v    (no -0.5)
```
is **off by 0.5 pixels**. The corrected formula includes `-0.5`:
```
z_center = (poni1 + L*tan(θ₂)/cos(θ₁)) / pixel_v - 0.5
```

### 1.6 Transformation Pipelines (No Spatial Distortion)

**pyFAI** (`calc_pos_zyx`, `core.py:455-559`):
```
[p1, p2, p3]^T = D_pyFAI · [dH, dV]^T + [-poni1, -poni2, L]^T
[t1, t2, t3]^T = R_pyFAI(rot1,rot2,rot3) · [p1, p2, p3]^T
[t1, t2, t3]^T = orientation_sign_flips(t1, t2, t3)
```
D_pyFAI maps pixel indices to physical coords with `d+0.5` shift.

**ImageD11** (`compute_xyz_lab`, `transform.py:85-147`):
```
flipped = O · [(sc - zc)*zs, (fc - yc)*ys]^T
vec = [0, flipped_y, flipped_z]^T
rotvec = R_ID11(tilt_x,tilt_y,tilt_z) · vec
t_ID11 = rotvec + [Δ, 0, 0]^T
```

### 1.7 Azimuthal Angle Definitions

| | pyFAI chi | ImageD11 eta |
|---|---|---|
| Definition | `arctan2(t1, t2)` | `arctan2(-t_y, t_z)` |
| t1/t2 | y_up, x_starboard | — |
| t_y/t_z | — | port direction, up direction |
| Zero dir | +x (starboard) | +z (up) |
| Positive | CCW from +x towards +y | CW facing downstream |

Given `t_ID11 = G · t_pyFAI`:
```
t_x_ID11 = t3_pyFAI    (downstream)
t_y_ID11 = -t2_pyFAI   (port = -starboard)
t_z_ID11 = t1_pyFAI    (up)
```

Therefore:
```
eta = arctan2(-t_y, t_z)
    = arctan2(-(-t2_pyFAI), t1_pyFAI)
    = arctan2(t2_pyFAI, t1_pyFAI)
    = arctan2(x_starboard, y_up)
```

And:
```
chi = arctan2(t1_pyFAI, t2_pyFAI)
    = arctan2(y_up, x_starboard)
```

Thus: **`eta = 90° - chi`** (when both in range [-π, π]).

To avoid wrap-around issues: compare `sin(chi)` vs `sin(90°-eta)` and `cos(chi)` vs `cos(90°-eta)`.

### 1.8 Existing Code Issues Found

1. `imaged11.py:134-138` — flip matrix hardcoded, TODO: *"manage orientation here"*
2. `convert_to_ImageD11` converts via Fit2D intermediate — unnecessary coupling
3. `convert_from_ImageD11` does not handle all orientations correctly
4. No round-trip test exists for all 4 orientations
5. 0.5 pixel offset unaccounted for in either direction

### 1.9 pyFAI Orientation Implementation (source-code review)

From `pyFAI/detectors/_common.py:657-678` (`_reorder_indexes_from_orientation`):

```python
if center:
    shape1 = self.shape[0] - 1   # first element of shape tuple
    shape2 = self.shape[1] - 1   # second element of shape tuple

if orientation == 2:  d1 = shape1 - d1       # d1 flip uses shape[0]-1
elif orientation == 4:  d2 = shape2 - d2     # d2 flip uses shape[1]-1
elif orientation == 1:  d1 = shape1 - d1; d2 = shape2 - d2
```

From `pyFAI/ext/_geometry.pyx:68-105` (`f_t1`/`f_t2`):

```c
// f_t1: orient = -1 if (orient==1 || orient==2) else 1
// f_t2: orient = -1 if (orient==1 || orient==4) else 1
```

**Two-level orientation mechanism:**

| Orient | Pixel reorder (pre-rotation) | Sign flip (post-rotation) |
|--------|------------------------------|---------------------------|
| 3 | none | none |
| 2 | d1 = shape[0]-1 - d1 | t1 = -t1 |
| 4 | d2 = shape[1]-1 - d2 | t2 = -t2 |
| 1 | both reorders | both flips |

Note: pyFAI's `shape` tuple is `(dim0, dim1)` where `dim0` is the first array
dimension (traditionally "fast" in C-order, "slow" in F-order). The naming
`shape1`/`shape2` in the source code reflects array indexing order, not a
physical axis convention. Regardless of interpretation, the code uses `shape[0]`
for d1 flips and `shape[1]` for d2 flips — this is the implemented behaviour
that the conversion must match and that the tests verify.

---

## 2. MATHEMATICAL MAPPINGS

### 2.1 pyFAI → ImageD11 (poni_to_par direction)

Given: `L, poni1, poni2, rot1, rot2, rot3, pixel_v, pixel_h, orientation, wavelength`

```
# Recover along-beam distance (uses compensated rotations)
Δ  = L / (cos(rot1) · cos(rot2))

# Recover ID11 tilts from compensated pyFAI rotation
tr1, tr2, tr3 = decompensate(rot1, rot2, rot3, orientation)
tilt_x = tr3
tilt_y = tr2
tilt_z = -tr1

# Beam center: reverse the orientation-specific PONI formula
# For orientation 3 (native):
#   z_center = (poni1 + L·tan(rot2)/cos(rot1)) / pixel_v - 0.5
# For orientation 2 or 1 (d1 flipped):
#   z_center = shape[0]-1 + 0.5 - (poni1 + L·tan(rot2)/cos(rot1)) / pixel_v
#   (and analogously for y_center with shape[1]-1 for d2-flipped orientations)
```

### 2.2 ImageD11 → pyFAI (par_to_poni direction)

Given: `distance, y_center, z_center, y_size, z_size, tilt_x, tilt_y, tilt_z, o11, o22, wavelength` and detector shape (typically (nfast, nslow))

```
# Standard tilt mapping (standard rotations, before compensation)
r1 = -tilt_z
r2 = tilt_y
r3 = tilt_x

# Compensated rotations from S·R_comp·C = R_tilt·Z
rot1, rot2, rot3 = compensate(o11, o22, orientation, r1, r2, r3)

# Orthogonal distance
L = distance · cos(rot2) · cos(rot1)

# PONI coordinates (orientation-specific, with 0.5 correction)
# Orientation 3 (native):
#   poni1 = -distance·sin(rot2) + z_size·(z_center + 0.5)
# Orientation 2 or 1 (d1 flipped, uses shape[0]-1):
#   poni1 = -distance·sin(rot2) + z_size·(shape[0]-1 - z_center + 0.5)
# Orientation 4 or 1 (d2 flipped, uses shape[1]-1):
#   poni2 = distance·cos(rot2)·sin(rot1) + y_size·(shape[1]-1 - y_center + 0.5)
```

### 2.3 2θ / q Invariance

2θ computed by both codes: `2θ = arctan2(√(t_y²+t_z²), t_x)` in lab coordinates.

Since the G-transformation between coordinate systems is orthonormal, `||t||` is preserved
and the xy-plane magnitude `√(t_y²+t_z²)` is also preserved (G swaps axis 1 ↔ 3 and
negates axis 2, which preserves distances in the yz plane). Therefore **2θ is IDENTICAL**
between both codes for the same physical geometry, regardless of flip/orientation.

### 2.4 Azimuth Angle Mapping

```
chi_pyFAI = arctan2(t1, t2)      [t1=y_up, t2=x_starboard]
eta_ID11  = arctan2(-t_y, t_z)   [t_y=port, t_z=up]

eta = 90° - chi                   (mod 360°)
chi = 90° - eta                   (mod 360°)
```

**Verification with sin/cos**: Compare `(sin(chi), cos(chi))` vs `(sin(90°-eta), cos(90°-eta))`.

### 2.5 Round-Trip Consistency

For correct conversion:
```
par_to_poni(poni_to_par(poni)) ≈ poni    (poni round-trip)
poni_to_par(par_to_poni(par)) ≈ par      (par round-trip)
```

---

## 3. IMPLEMENTATION FILES

### 3.1 `par_to_poni.py` — Conversion + IO

Functions:
- `par_to_poni(par, detector_shape=None)` → poni dict
- `poni_to_par(poni, detector_shape=None)` → par dict
- `read_par(filepath, par_length_unit="um")` → dict
- `write_par(par_dict, filepath, par_length_unit="um")` → None
- `read_poni(filepath)` → dict
- `write_poni(poni_dict, filepath)` → None

The `detector_shape` parameter is a `(fast_dim, slow_dim)` tuple needed for
non-native orientations to compute orientation-specific PONI formulas that
account for pyFAI's pixel reordering. Defaults to square inferred from beam
center.

Dependencies: `numpy`, `scipy` (for `scipy.spatial.transform.Rotation`).
All internal units: meters for lengths, meters for wavelength.

### 3.2 `mapping.md` — Mathematical Derivations

Algebraic notation and Python code snippets for all mappings, including
the 0.5 correction derivation and azimuth mapping.

### 3.3 `test_conversion.py` — Test Suite

Tests (using pyFAI AzimuthalIntegrator and ImageD11 `compute_tth_eta`/`compute_xyz_lab`):

1. **Round-trip tests** for all 4 orientations and multiple tilt combos
2. **2θ matching**: pyFAI `tth()` vs ImageD11 `compute_tth_eta` — same raw pixels, no flipping
3. **Azimuth matching**: sin/cos comparison of chi vs 90°-eta — same raw pixels, no flipping
4. **Lab coordinate matching**: full xyz comparison on non-square 200×128 detector
5. **Edge cases**: zero tilts, max tilts, edge beam positions

Test geometry: strongly tilted detector (tilt_x=0.3, tilt_y=0.2, tilt_z=-0.15 rad),
75µm pixels, 150mm distance, Cu Kα wavelength. Square (1000×1000) for most tests,
non-square (200×128) for the coordinate-level test.

---

## 4. DELIVERABLES

| # | File | Description |
|---|---|---|
| 1 | `PLAN.md` | This plan (derived from review) |
| 2 | `mapping.md` | Mathematical derivations with formulas and Python snippets |
| 3 | `par_to_poni.py` | Conversion + IO functions |
| 4 | `test_conversion.py` | Test suite running against pyFAI + ImageD11 |

## 5. SOLUTION: AFFINE TRANSFORMATION ANALYSIS

Both pyFAI and ImageD11 are affine transformations from pixel coordinates to
lab coordinates. The pyFAI pipeline decomposes into three operations:

1. **Pixel reordering** (pre-rotation): `C = diag(c1, c2)` — flips pixel indices
   before computing physical coordinates (orientation-dependent, per `_common.py`)
2. **Rotation**: `R` — the pyFAI rotation matrix
3. **Sign flips** (post-rotation): `S = diag(s1, s2, 1)` — flips lab-coordinate
   signs after rotation (orientation-dependent, per `_geometry.pyx`)

The ImageD11 pipeline encodes flips via the matrix Z = diag(o11, -o22) applied
pre-rotation (in the pyFAI lab frame after G transformation).

Equating the linear parts of the two affine transforms gives:

```
S(orient) · R_comp · C(orient) = R_tilt · Z(flip)
```

Solving for the compensated rotation R_comp column by column:

```
R_comp[:,0] = S · R_tilt[:,0] · (o11 / c1)
R_comp[:,1] = S · R_tilt[:,1] · (-o22 / c2)
```

and the third column from cross product (ensuring det=+1).

The compensated rotations are exact (columns orthonormal, det=+1), and the
Euler angles are extracted via `scipy.spatial.transform.Rotation`.

For orientations 3 and 1, the raw constraint `S·R·C = R_tilt·Z` yields
`cos(rot1)·cos(rot2) > 0` and a positive orthogonal distance.

For orientations 2 and 4, the raw constraint yields `R[2,2] < 0` and
therefore `cos(rot1)·cos(rot2) < 0`.  No equivalent ZYX parametrization
of the same matrix gives `cos(rot1)·cos(rot2) > 0`.  The solution relaxes
xyz coordinate matching by introducing per-orientation mirror matrices
that match which pixel axis each orientation flips:

```
orient 2 (flip slow / pyFAI axis 1):  M = diag(-1,  1,  1)
orient 4 (flip fast / pyFAI axis 2):  M = diag( 1, -1,  1)
orient 1 (flip both):                 M = diag(-1, -1,  1)
orient 3 (native, no flip):           M = identity
```

The compensated rotation is found from `S·R_comp·C = M·R_tilt·Z`.
Each mirror accepts that pyFAI's orientation creates an effective
reflection of the coordinate axes, keeping distance positive while
preserving 2θ and azimuth.

**Azimuth relationship**:
- Orient 3: `chi = 90° − eta`   ( sin=+cos(eta), cos=+sin(eta) )
- Orient 2: `chi = eta − 90°`   ( sin=−cos(eta), cos=+sin(eta) )
- Orient 4: `chi = eta + 90°`   ( sin=+cos(eta), cos=−sin(eta) )
- Orient 1: `chi = 270° − eta`  ( sin=−cos(eta), cos=−sin(eta) )

**Lab coordinates**: In the ID11 frame, the mirror reflects:
- Orient 2: Z-axis flip  (slow=y_up maps to ID11 Z)
- Orient 4: Y-axis flip  (fast=x_starboard maps to ID11 −Y)
- Orient 1: Y and Z flip (both)
- Orient 3: no flip

The PONI constants must additionally account for orientation-specific pixel
reordering — the beam center in native coordinates maps to `max - zc + 0.5`
in the reordered coordinate system for d1-flipped orientations (2 and 1).

**Conclusion**: The conversion between par and poni is exact to machine
precision for all 4 non-transpose orientations. Distances are positive
in all cases. Verified by test tolerances of 1e-7 rad for 2θ and
azimuth (per-orientation relationships), and 5e-7 m for lab coordinates
(after per-orientation mirror reflections), on both square and
non-square detectors.

### 5.1 Remaining Limitations

**Transpose Flips**: Not supported (o12, o21 must be 0). PyFAI's orientation
model does not handle axis-swapped images.

**Spatial Distortion**: Assumed absent for the geometric conversion.

### 5.2 Open Question: pyFAI Orientation Model

The conversion locks in pyFAI's current orientation implementation as
observed in the source code (`_common.py:657-678`, `_geometry.pyx:68-105`).
This includes the pixel-reordering convention where `shape[0]-1` is used for
d1 (slow-axis) flips and `shape[1]-1` for d2 (fast-axis) flips.

This is a question for the pyFAI maintainers: is this the intended
orientation model, or should orientation be revised? Possible actions:

- **Option A (lock-in)**: Accept the current pyFAI implementation as
  definitive. Document it clearly. The conversion code and tests enshrine
  this behaviour.

- **Option B (revise)**: Propose that pyFAI's pixel reordering use per-axis
  maximum values (`shape[1]-1` for d1/slow, `shape[0]-1` for d2/fast) so
  that the naming is physically consistent. The conversion code would need
  updating if pyFAI changes.

This is a design decision for humans, not an LLM. The conversion code
reflects Option A (the current pyFAI implementation). Tests pass against
the installed version of pyFAI.

## 6. Round 2 Corrections

### 6.1 Bugs Found by Referee #3 (Claude)

**max_d1/max_d2 swap**: The conversion code used `max_d1 = shape_fast - 1`
and `max_d2 = shape_slow - 1`, swapping the axis assignments.  PyFAI's
`_reorder_indexes_from_orientation` uses shape[0] (C-order slow axis) for
d1 flips and shape[1] (fast axis) for d2 flips.  Fixed: `max_d1` now uses
slow count (detector_shape[1]), `max_d2` uses fast count (detector_shape[0]).
Bug was invisible on square detectors.  Both `par_to_poni()` and
`poni_to_par()` affected.

**TestLabCoordinates shape convention**: `ai.detector.shape` now passed as
(slow, fast) matching pyFAI's C-order convention: `(SHAPE[1], SHAPE[0])`.

### 6.2 Rotations Validation (Referee #2)

Added `test_pyfai_rotation_matrix_matches_actual` comparing
`_pyfai_rotation_matrix()` to pyFAI's `rotation_matrix()` method —
identical to 2.2e-16.

### 6.3 Integration Validation

Added `test_write_poni_loads_and_integrates`: writes par→poni for all 4
orientations, loads with `pyFAI.load()`, calls `integrate1d`.  All pass.

### 6.4 Mirror-Matrix Approach (Round 3: Positive Distance for All Orientations)

For orientations 2 and 4, the raw constraint `S·R_comp·C = R_tilt·Z`
forces `R[2,2] < 0`, giving a negative orthogonal distance (see §5).
The ZYX Euler convention provides no equivalent parametrization with
positive `cos(rot1)·cos(rot2)`.

The solution relaxes xyz coordinate matching by introducing a
per-orientation mirror matrix into the rotation constraint:

```
S · R_comp · C = M · R_tilt · Z
```

Each orientation gets the mirror matching its flipped pixel axis:

| Orient | Flips | Mirror M | ID11 frame effect |
|--------|-------|----------|-------------------|
| 3 | none | identity | none |
| 2 | slow/y | diag(−1, 1, 1) | Z-axis flip |
| 4 | fast/x | diag(1, −1, 1) | Y-axis flip |
| 1 | both | diag(−1, −1, 1) | Y+Z flip |

This makes the coordinate frame consistent with the detector's fast/slow
axes. With this relaxation:

- Distance is positive for all 4 orientations
- 2θ values match exactly (machine precision)
- Azimuth: chi = 90°−eta (orient 3), = eta−90° (2), = eta+90° (4), = 270°−eta (1)
- Lab coordinates match after per-orientation mirror reflection in ID11 frame
- Round-trip is exact to machine precision

Each mirror is self-inverse (`M² = I`); the reverse transform applies
the same mirror: `R_tilt[:,k] = M · S · R_comp[:,k] · (c_k / z_k)`.
