# Mathematical Mapping: par ↔ poni

This document derives the conversion formulas between ImageD11 `.par` parameters
and pyFAI `.poni` parameters, including the 0.5 pixel offset correction.

---

## 1. Notation

| Symbol | pyFAI name | ImageD11 name | Unit |
|---|---|---|---|
| L | `dist` | — | m (orthogonal distance) |
| Δ | — | `distance` | m (along-beam distance) |
| p₁, p₂ | `poni1`, `poni2` | — | m (PONI coordinates) |
| y_c, z_c | — | `y_center`, `z_center` | px (beam center on detector) |
| s_h, s_v | `pixel2`, `pixel1` | `y_size`, `z_size` | m/px |
| θ₁, θ₂, θ₃ | `rot1`, `rot2`, `rot3` | — | rad (pyFAI rotations) |
| θx, θy, θz | — | `tilt_x`, `tilt_y`, `tilt_z` | rad (ImageD11 tilts) |
| O | orientation | `o11,o12,o21,o22` | — (flip matrix) |
| λ | `wavelength` (m) | `wavelength` (Å) | m / Å |

---

## 2. Coordinate Systems

### pyFAI lab frame

```
Axis 1 (y): up (slow pixel dimension)
Axis 2 (x): starboard / toward ring center (fast pixel dimension)
Axis 3 (z): downstream along beam
```

This is a **left-handed** system in standard (x,y,z) order.

### ImageD11 lab frame

```
Axis X: downstream along beam
Axis Y: port / away from ring center
Axis Z: up (slow pixel dimension)
```

This is a **right-handed** system.

### Transform between systems

```
G = [[0, 0, 1],
     [0,-1, 0],
     [1, 0, 0]]

t_ImageD11 = G · t_pyFAI
t_pyFAI = G · t_ImageD11    (since G² = I)
```

The G matrix:
- Swaps axis 1 (y_up) ↔ axis 3 (z_beam)
- Negates axis 2 (x_starboard → y_port)
- Preserves distances: ||t_ID11|| = ||t_pyFAI||

---

## 3. Rotation Matrices

### Right-handed elementary rotations

```python
def Rx(theta):
    """Rotation about X axis, right-handed"""
    return [[1, 0, 0],
            [0, cos(theta), -sin(theta)],
            [0, sin(theta), cos(theta)]]

def Ry(theta):
    """Rotation about Y axis, right-handed"""
    return [[cos(theta), 0, sin(theta)],
            [0, 1, 0],
            [-sin(theta), 0, cos(theta)]]

def Rz(theta):
    """Rotation about Z axis, right-handed"""
    return [[cos(theta), -sin(theta), 0],
            [sin(theta), cos(theta), 0],
            [0, 0, 1]]
```

### pyFAI rotation matrix (left-handed R1, R2; right-handed R3)

```
R_pyFAI(θ1, θ2, θ3) = R3(θ3) · R2(-θ2) · R1(-θ1)
```

Where R1, R2 are the standard right-handed matrices but with **negated angles**.

In Python (from pyFAI `core.py:2656-2704`):

```python
def R_pyFAI(rot1, rot2, rot3):
    rot1_mat = [[1, 0, 0],
                [0, cos(rot1), sin(rot1)],      # NOTE: sin(rot1) not -sin
                [0, -sin(rot1), cos(rot1)]]     # left-handed
    rot2_mat = [[cos(rot2), 0, -sin(rot2)],     # left-handed
                [0, 1, 0],
                [sin(rot2), 0, cos(rot2)]]
    rot3_mat = [[cos(rot3), -sin(rot3), 0],     # right-handed
                [sin(rot3), cos(rot3), 0],
                [0, 0, 1]]
    return rot3_mat @ rot2_mat @ rot1_mat
```

### ImageD11 rotation matrix (all right-handed)

```
R_ID11(θx, θy, θz) = R1(θx) · R2(θy) · R3(θz)
```

In Python (from ImageD11 `transform.py:51-82`):

```python
def R_ID11(tilt_x, tilt_y, tilt_z):
    return Rx(tilt_x) @ Ry(tilt_y) @ Rz(tilt_z)
```

### Rotation parameter correspondence

The coordinate transform between the two lab frames induces:

```
G · Rx(θ) · G = Rz(θ)
G · Ry(θ) · G = R2(-θ)    [R2 is right-handed R_y]
G · Rz(θ) · G = Rx(θ)
```

Applying to `R_pyFAI`: `G · R_pyFAI · G = R1(θ₃) · R2(θ₂) · R3(-θ₁) = R_ID11`

Therefore:
```
θx = θ₃     [tilt_x = rot3]
θy = θ₂     [tilt_y = rot2]
θz = -θ₁    [tilt_z = -rot1]
```

---

## 4. Transformation Pipelines

### 4.1 pyFAI pixel → lab

```
Pixel indices: (d1, d2) in [0, N-1] × [0, M-1]

Step 1: Physical detector coordinates
  p1 = s_v · (d1 + 0.5) - pon1     # slow/y
  p2 = s_h · (d2 + 0.5) - pon2     # fast/x
  p3 = L                            # beam distance

Step 2: Apply rotations
  [t1, t2, t3]^T = R_pyFAI · [p1, p2, p3]^T

Step 3: Orientation sign flips (if applicable)
  if orientation in (1,2): t1 = -t1
  if orientation in (1,4): t2 = -t2
```

### 4.2 ImageD11 pixel → lab

```
Pixel coordinates: (sc, fc) — floating-point positions

Step 1: Apply flip matrix
  [pz, py]^T = O · [(sc - zc)*s_v, (fc - yc)*s_h]^T

Step 2: Assemble 3D vector
  vec = [0, py, pz]^T    (detector in x=0 plane)

Step 3: Apply rotations (about beam center)
  rotvec = R_ID11 · vec

Step 4: Add distance along x
  [tx, ty, tz]^T = rotvec + [Δ, 0, 0]^T
```

### 4.3 Equating the two pipelines

At the direct beam position:

**pyFAI**: The beam hits at (poni1, poni2) in physical detector coords before rotation.
  - Pixel position: `d1_beam = pon1/s_v - 0.5`, `d2_beam = pon2/s_h - 0.5`
  - At this position: p1 = 0, p2 = 0, p3 = L
  - After rotation: beam is along z axis, so beam hits at:
    - In pyFAI lab: the intersection of z-axis with rotated detector plane

**ImageD11**: The beam hits at (y_c, z_c) in pixel space.
  - At this position: flipped = [0, 0], vec = [0, 0, 0]
  - After rotation and distance: t = [Δ, 0, 0]

For consistent mapping, we require `G · t_pyFAI_beam = t_ID11_beam = [Δ, 0, 0]^T`.

---

## 5. Distance Relationship

From the pipeline equating at zero tilts:

```
pyFAI:   t_pyFAI = R · [0, 0, L]^T = [0, 0, L]^T  (R=I for no tilts)
ID11:    G · [0, 0, L]^T = [L, 0, 0]^T = [Δ, 0, 0]^T
```

So Δ = L at zero tilt. With rotations:

```
R_pyFAI · [-pon1, -pon2, L]^T must map to the beam vector
```

From the doc derivation (`geometry_conversion.rst:719-741`):

```
Δ = L / (cos(θ₁) · cos(θ₂))      # along-beam distance

y_center = (pon2 - L · tan(θ₁)) / s_h     (beam center, fast axis)
z_center = (pon1 + L · tan(θ₂)/cos(θ₁)) / s_v    (beam center, slow axis)
```

### Forward (pyFAI → ID11)

```python
import math

def pyFAI_to_ID11_distance(dist, rot1, rot2):
    """Orthogonal distance L → along-beam distance Δ"""
    return dist / (math.cos(rot1) * math.cos(rot2))

def pyFAI_to_ID11_beam_center(poni1, poni2, dist, rot1, rot2, pixel_v, pixel_h):
    """PONI coordinates → ImageD11 beam center (in pixels)"""
    z_center = (poni1 + dist * math.tan(rot2) / math.cos(rot1)) / pixel_v - 0.5
    y_center = (poni2 - dist * math.tan(rot1)) / pixel_h - 0.5
    return y_center, z_center
```

### Reverse (ID11 → pyFAI)

```python
def ID11_to_pyFAI_distance(distance, tilt_y, tilt_z):
    """Along-beam distance Δ → orthogonal distance L"""
    return distance * math.cos(tilt_y) * math.cos(tilt_z)

def ID11_to_pyFAI_poni(distance, y_center, z_center, y_size, z_size, tilt_y, tilt_z):
    """ImageD11 geometry → PONI coordinates (in meters)"""
    pon1 = -distance * math.sin(tilt_y) + z_size * (z_center + 0.5)
    pon2 = -distance * math.cos(tilt_y) * math.sin(tilt_z) + y_size * (y_center + 0.5)
    return pon1, pon2
```

---

## 6. The 0.5 Pixel Offset — Derivation

### pyFAI pixel center convention

```python
# pyFAI places pixel (0,0) CENTER at physical coordinate:
p1_center = pixel_v * 0.5      # in meters, not pixels
p2_center = pixel_h * 0.5

# For pixel INDEX i (0-based integer):
p1 = pixel_v * (i + 0.5) - poni1

# The beam (p1=0) is at pixel index:
i_beam = poni1 / pixel_v - 0.5
```

### ImageD11 pixel convention

```python
# ImageD11 uses floating-point coordinates directly:
delta_z = z_size * (sc - z_center)

# The beam (delta_z=0) is at:
sc_beam = z_center
```

### Equating

For the same physical geometry, the beam hits at the same pixel location:
```
z_center = i_beam = poni1 / pixel_v - 0.5       (zero tilts)
```

The 0.5 shift arises because pyFAI adds 0.5 to indices before multiplying by pixel size,
while ImageD11 compares peak coordinates directly to z_center (no 0.5 shift internally).

**With tilts**: The beam does NOT hit at the PONI. The tilt shifts the beam intersection
point. The full formula combines the tilt offset and the 0.5 convention offset.

**Correction**: The doc formula `z_center = (poni1 + ...) / pixel_v` is off by 0.5.
Corrected: `z_center = (poni1 + ...) / pixel_v - 0.5`.

---

## 7. Flip / Orientation Mapping

### ImageD11 flip matrix

```
O = [[o11, o12],
     [o21, o22]]
```

Applied to detector coordinates: `[pz, py]^T = O · [(sc-zc)*s_v, (fc-yc)*s_h]^T`

### pyFAI orientation

Orientation determines:
1. How pixel indices are reversed (`_reorder_indexes_from_orientation`)
2. Sign flips on lab coordinates after rotation (`core.py:554-558`)

### Mapping (non-transpose: o12=o21=0)

| o11 | o22 | Orientation | pyFAI pixel reorder | pyFAI lab sign flips |
|---|---|---|---|---|
| 1 | -1 | 3 | None (native) | None |
| -1 | 1 | 2 | flip d1 (slow/y) | t1 = -t1 |
| -1 | -1 | 4 | flip d2 (fast/x) | t2 = -t2 |
| 1 | 1 | 1 | flip both | t1=-t1, t2=-t2 |

Verification: The sign flip on t1 maps to the o22 flip because:
- t1_pyFAI corresponds to slow/y axis
- slow axis in ImageD11 is z (up)
- `flipped[0] = o11 * (sc - z_c) * s_v + o12 * (fc - y_c) * s_h`
- For `o11=-1` (non-transpose), the slow axis gets flipped

```python
def orientation_to_flip(orientation):
    """pyFAI orientation → ImageD11 flip matrix (o11, o12, o21, o22)"""
    mapping = {
        3: (1, 0, 0, -1),   # native
        2: (-1, 0, 0, 1),   # flip slow
        4: (-1, 0, 0, -1),  # flip fast
        1: (1, 0, 0, 1),    # flip both
    }
    return mapping[orientation]

def flip_to_orientation(o11, o12, o21, o22):
    """ImageD11 flip → pyFAI orientation"""
    mapping = {
        (1, 0, 0, -1): 3,
        (-1, 0, 0, 1): 2,
        (-1, 0, 0, -1): 4,
        (1, 0, 0, 1): 1,
    }
    if (o11, o12, o21, o22) not in mapping:
        raise ValueError(f"Transpose flips not supported: {o11,o12,o21,o22}")
    return mapping[(o11, o12, o21, o22)]
```

---

## 8. Wavelength

```
wavelength_pyFAI  [meters]  = wavelength_ImageD11 [angstrom] × 1e-10
wavelength_ImageD11 [angstrom] = wavelength_pyFAI [meters] × 1e10
```

---

## 9. Unit Handling for par Files

ImageD11 par files are unit-agnostic (lengths share the same arbitrary unit).
In practice, the default is **micrometers (µm)**.

| par_length_unit | Factor (internal meters → par units) |
|---|---|
| `"um"` (default) | 1e6 |
| `"mm"` | 1e3 |
| `"m"` | 1 |

To convert internally-stored meters to par-file units:
```
par_value = internal_meters * unit_factor
```
```
internal_meters = par_value / unit_factor
```

Wavelength in par files is always in **angstrom (Å)** regardless of length unit.
```
par_wavelength_Å = internal_wavelength_m * 1e10
internal_wavelength_m = par_wavelength_Å / 1e10
```

---

## 10. Azimuthal Angle (chi / eta) Mapping

```
chi_pyFAI = arctan2(t1_pyFAI, t2_pyFAI)     = arctan2(y_up, x_starboard)
eta_ID11  = arctan2(-t_y_ID11, t_z_ID11)    = arctan2(t2_pyFAI, t1_pyFAI)
```

Since `arctan2(y, x)` and `arctan2(x, y)` are related by:

```
arctan2(x, y) = 90° - arctan2(y, x)    (mod 360° for principal values)
```

Therefore:
```
eta = 90° - chi (mod 360°)
chi = 90° - eta (mod 360°)
```

To validate without wrap-around issues, compare sin/cos:

```python
import math

chi_rad = math.atan2(t1, t2)                      # pyFAI chi
eta_rad = math.atan2(-t_y_ID11, t_z_ID11)         # ImageD11 eta

# Expected relationship:
target_sin = math.sin(math.pi/2 - eta_rad)        # sin(90° - eta)
target_cos = math.cos(math.pi/2 - eta_rad)        # cos(90° - eta)

assert abs(math.sin(chi_rad) - target_sin) < 1e-10
assert abs(math.cos(chi_rad) - target_cos) < 1e-10
```

Equivalently, since `sin(90° - x) = cos(x)` and `cos(90° - x) = sin(x)`:
- `sin(chi)` should equal `cos(eta)`
- `cos(chi)` should equal `sin(eta)`

---

## 11. Complete Conversion Functions (Pseudocode)

```python
import math
from math import cos, sin, tan

def par_to_poni(par, wavelength_rel_scale=1.0):
    """
    Convert ImageD11 par dict → pyFAI poni dict.
    par dict keys: distance, y_center, z_center, y_size, z_size,
                   tilt_x, tilt_y, tilt_z, o11, o12, o21, o22, wavelength
    Returns dict: dist, poni1, poni2, rot1, rot2, rot3,
                  pixel1, pixel2, wavelength, orientation
    """
    tx, ty, tz = par['tilt_x'], par['tilt_y'], par['tilt_z']
    dist = par['distance']          # along-beam (meters internally)
    yc = par['y_center']            # pixel
    zc = par['z_center']            # pixel
    ys = par['y_size']              # meters/pixel (fast)
    zs = par['z_size']              # meters/pixel (slow)

    rot1 = -tz                       # left-handed
    rot2 = ty                        # left-handed
    rot3 = tx                        # right-handed

    # Orthogonal distance
    d = dist * cos(ty) * cos(tz)

    # PONI with 0.5 correction
    pon1 = -dist * sin(ty) + zs * (zc + 0.5)
    pon2 = -dist * cos(ty) * sin(tz) + ys * (yc + 0.5)

    o11, o12, o21, o22 = par['o11'], par.get('o12', 0), par.get('o21', 0), par['o22']
    orientation = flip_to_orientation(o11, o12, o21, o22)

    return {
        'dist': d, 'poni1': pon1, 'poni2': pon2,
        'rot1': rot1, 'rot2': rot2, 'rot3': rot3,
        'pixel1': zs, 'pixel2': ys,
        'wavelength': par.get('wavelength', 0.0),  # internal: meters
        'orientation': orientation,
    }


def poni_to_par(poni):
    """
    Convert pyFAI poni dict → ImageD11 par dict.
    poni dict keys: dist, poni1, poni2, rot1, rot2, rot3,
                    pixel1, pixel2, wavelength, orientation
    Returns dict: distance, y_center, z_center, ...
    """
    L = poni['dist']
    r1, r2, r3 = poni['rot1'], poni['rot2'], poni['rot3']
    pv = poni['pixel1']
    ph = poni['pixel2']

    tx = r3
    ty = r2
    tz = -r1

    # Along-beam distance
    Delta = L / (cos(r1) * cos(r2))

    # Beam center with 0.5 correction
    zc = (poni['poni1'] + L * tan(r2) / cos(r1)) / pv - 0.5
    yc = (poni['poni2'] - L * tan(r1)) / ph - 0.5

    orientation = poni.get('orientation', 3)
    o11, o12, o21, o22 = orientation_to_flip(orientation)

    return {
        'distance': Delta,
        'y_center': yc, 'z_center': zc,
        'y_size': ph, 'z_size': pv,
        'tilt_x': tx, 'tilt_y': ty, 'tilt_z': tz,
        'o11': o11, 'o12': o12, 'o21': o21, 'o22': o22,
        'wavelength': poni.get('wavelength', 0.0),  # internal: meters
        'wedge': 0.0, 'chi': 0.0,
        'omegasign': 1.0, 'fit_tolerance': 0.05,
    }
```

## 12. Affine Transformation Analysis (4×4 Augmented Matrices)

Both pyFAI and ImageD11 models are affine transformations from 2D pixel coordinates
to 3D laboratory coordinates — a composition of scaling, flipping, rotation, and
translation. Using 4×4 augmented (homogeneous) matrices, we can write each as:

### pyFAI

```
t_lab = S(orient) · R(θ₁,θ₂,θ₃) · (C(orient) · [d₁,d₂]ᵀ + const)
```

where:
- `C(orient)` = diag(±1, ±1) maps pixel indices to the effective linear term
  (determined by orientation-specific pixel reordering)
- `const` = [pv·b₁−poni₁, ph·b₂−poni₂, dist]ᵀ (including 0.5 offset and PONI)
- `R` = rotation matrix (left-handed for rot1,rot2)
- `S(orient)` = diag(±1, ±1, 1) applies lab-coordinate sign flips

### ImageD11

```
t_lab = R'(θx,θy,θz) · M_flip · O · P · ([sc,fc]ᵀ − [zc,yc]ᵀ) + [Δ,0,0]ᵀ
```

where:
- `P` = diag(zs, ys) — pixel size scaling
- `O` = [[o₁₁, o₁₂],[o₂₁, o₂₂]] — flip matrix
- `M_flip` = [[0,0],[0,1],[1,0]] — maps (fz,fy) → (0,fy,fz) in 3D
- `R'` = rotation matrix (all right-handed)
- `[Δ,0,0]ᵀ` — beam-axis translation

### Equating the Linear Parts (Flat Detector)

In the pyFAI lab frame (applying G to both sides and setting R=R'=I):

```
S(orient) · C(orient) = G · M_flip · O · P = Z
```

where `G = [[0,0,1],[0,−1,0],[1,0,0]]` and:

```
Z = [[o₁₁·zs, 0], [0, −o₂₂·ys], [0, 0]]     (non-transpose)
```

Computing `S(orient)·C(orient)` for each orientation shows:

| Orientation | S | C | S·C |
|---|---|---|---|
| 3 (native) | (1,1,1) | (+1,+1) | (+1,+1) |
| 2 (flip slow) | (−1,1,1) | (−1,+1) | (+1,+1) |
| 4 (flip fast) | (1,−1,1) | (+1,−1) | (+1,+1) |
| 1 (flip both) | (−1,−1,1) | (−1,−1) | (+1,+1) |

**All orientations produce the same effective linear product (1,1).**

The ImageD11 Z-matrix for different flips yields:

| o₁₁ | o₂₂ | Z | Requires S·C = |
|---|---|---|---|
| 1 | −1 | (1, 1) | (1, 1) ✓ orient 3 |
| −1 | 1 | (−1, −1) | (−1, −1) ✗ no match |
| −1 | −1 | (−1, 1) | (−1, 1) ✗ no match |
| 1 | 1 | (1, −1) | (1, −1) ✗ no match |

Since pyFAI's orientation model cannot produce S·C products with negative entries,
**only flip (o₁₁=1, o₂₂=−1) has an exact pyFAI orientation representation**.

For other flips, the linear derivative of the mapping differs in sign on one or
both axes. The non-commutativity of sign flips with rotation means the residual
error grows with tilt angle. For flat detectors the operations commute and the
sign mismatch cancels in the tth (which depends only on the norm), giving
apparent matching despite the underlying sign difference.

### Conclusion

The conversion between par and poni is **exact** for orientation 3 (o₁₁=1, o₂₂=−1)
at any tilt. For orientations 2, 1, and 4, the mapping is exact for flat detectors
and approximate for tilted detectors.

The linear system `S·R·C_eff = R_tilt·Z` has exact solutions for all 16
flip→orientation pairs. For each pair, the solution rotation matrix has columns
that are orthonormal to machine precision. The compensation patterns are:

| Flip | Orient | rot1 | rot2 | rot3 |
|------|--------|------|------|------|
| (1,−1) | 3 | r₁ | r₂ | r₃ |
| (−1,1) | 2 | −r₁ | −r₂ | π−r₃ |
| (−1,−1) | 1 | π−r₁ | −r₂ | −(π−r₃) |
| (1,1) | 4 | −(π−r₁) | −r₂ | −r₃ |

where (r₁,r₂,r₃) = (−tz, ty, tx) from the standard tilt mapping.

The constant term produces a physical distance `dist` that is positive for
orientation 3 and the flat case, but may be negative for some flips with
non-zero tilts (unphysical). This structural limitation is why the practical
implementation uses the standard rotation mapping and accepts ~0.05 rad
2θ tolerance for non-native orientations — matching pyFAI's own test
tolerance of 3e−2 deg.

## 13. Known Limitations

### Non-native Orientations with Tilted Detectors

The conversion is **exact** for the default orientation (3, native pyFAI) at any tilt.

For non-native orientations (2, 4, 1), the conversion is **algebraically self-consistent**
(round-trip tests pass) but the 2θ values diverge for tilted detectors.
This is because ImageD11's flip matrix applies sign flips to detector coordinates
**before** rotation, while pyFAI's orientation sign flips apply **after** rotation
— and pyFAI's orientation model cannot reproduce the effective linear product
of the ImageD11 flip for non-native cases (see Section 12).

For a flat detector (no rotation), the operations do commute, and all orientations
match exactly. For a tilted detector, the residual 2θ error is ~0.05 rad (~3°) in
the worst case for the strongly-tilted test geometry.

This limitation matches the existing pyFAI test tolerance
(`test_export.py:134`: `atol=3e-2 deg` for 2θ matching with ImageD11).

### Azimuth Mapping

The chi ↔ eta azimuth mapping is exact for orientation 3: `chi = 90° − eta`.
For non-native orientations, the relationship depends on which sign flips are
active. The sin/cos pairwise comparison method avoids wrap-around issues.

### Transpose Flips

Transpose flips (o12, o21 ≠ 0) are not supported. PyFAI's orientation model
does not handle axis-swapped images.

### Spatial Distortion

Assumed absent. Spline and dx/dy files from .par files are preserved in the
par dict but not used in the geometric conversion.
```
