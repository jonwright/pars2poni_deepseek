# PLAN: All Solutions — Full 4×4 Cross-Mapping

## Goal

Enumerate ALL valid pyFAI poni solutions for a given ImageD11 par geometry,
covering all 16 (flip, orientation) pairings × 4 mirror matrices × 2 Euler
representations → 32 distinct solutions per flip.

## The Equation

For a given (flip F, orientation O, mirror M) triple:

    S(O) · R · C(O) = M · R_tilt · Z(F)

where:
- S(O) = pyFAI post-rotation sign flips (diag(s1, s2, 1))
- C(O) = pyFAI pixel reordering (diag(c1, c2))
- Z(F) = ImageD11 flip matrix extended to 3D (diag(o11, -o22, 1))
- M ∈ {M1, M2, M3=I, M4} = one of four mirror matrices

## Solution Structure

| Parameter | Count |
|-----------|-------|
| Orientations (O) | 4 |
| Mirrors (M) per orientation | 4 |
| Euler reps per (O, M) | 2 |
| **Total per flip** | **32** |

Across all 4 flips: 128 raw entries, deduplicating to 32 unique (rot, dist) tuples.

## Implementation

### `_azimuth_factors(mirror_orient)`
Returns the (sin_factor, cos_factor) for χ ↔ η mapping given a mirror.
The factors are the first two diagonal entries of M(mirror_orient):
```
sin(χ) = M[0,0] · cos(η)
cos(χ) = M[1,1] · sin(η)
```

### `find_all_poni_solutions(par, detector_shape, include_backscattering)`
Enumerates all 32 solutions for a par dict.
- Loops over all 4 orientations and all 4 mirror matrices
- For each, computes compensated rotation via `_compute_compensated_rotation()`
- Finds all ZYX Euler equivalents via `_find_all_rot_equivs()`
- Builds poni dicts with per-pairing azimuth factors
- Tags each solution with metadata: flip_label, orient_tried, mirror_source,
  chi_eta_exact, dist_positive, is_canonical, rot_magnitude
- Deduplicates by (rot1, rot2, rot3, dist)
- Sorts: positive distance first, then smallest magnitude, then canonical pairing

### Modified `par_to_poni()`
```
def par_to_poni(par, detector_shape=None,
                prefer_positive_distance=True, exact_chi=False)
```
- Default: canonical pairing with positive distance (unchanged from pre-refactor)
- `exact_chi=True`: picks solutions where mirror matches trial orientation (χ = 90° − η)
- Falls back to backscattering search when no positive-distance solution exists

### Modified `poni_to_par()`
- Reads `_mirror_orient`, `_forward_o11`, `_forward_o22` metadata from poni dict
- Uses stored flip values for reverse compensation instead of `orientation_to_flip()`
- Metadata persisted in `Detector_config` JSON for disk round-trips

## Test Coverage (TestAllSolutions)

- Verifies 32 solutions per flip, grouped into 16 (orient, mirror) pairs with 2 Euler reps each
- 2θ matches for ALL solutions (≤ 10⁻⁷ rad)
- Azimuth matches using correct mirror-dependent per-pairing factors
- Round-trip exact for all solutions (metadata preserves original flip)
- Euler pairs within each (orient, mirror) group produce identical rotation matrices
- `exact_chi=True` picks chi_eta_exact solutions
- Default API unchanged from pre-refactor
- Backscattering (6 additional tests): positive distance, π rotations, all 4 orientations

## Key Findings

1. **32 solutions per flip.** 4 orientations × 4 mirrors × 2 Euler reps.
   Each (orient, mirror) group has exactly 2 Euler representations
   (the ZYX β-solution pair). Both share the same R[2,2] = cos(rot1)·cos(rot2)
   and thus the same distance sign.

2. **Det(R) = +1 always.** The two columns r_c0, r_c1 are orthonormal,
   so det([r_c0, r_c1, cross(r_c0, r_c1)]) = |cross|² = 1 > 0. The
   det < 0 check in the code is unreachable dead logic — the cross
   product of two orthonormal vectors always gives a right-handed frame.

3. **Azimuth depends only on the mirror M**, not on orientation or flip.
   From t_pyFAI = M · t_id: χ = atan2(m0·t0, m1·t1), η = atan2(t1, t0).
   Therefore sin(χ) = m0·cos(η), cos(χ) = m1·sin(η).

4. **Backscattering (include_backscattering=True).** For negative ImageD11
   distance, π-offset seeds uncover positive-distance pyFAI representations
   with rot1≈π or rot2≈π. Works for zero tilts; approximate for tilted
   backscattering (physically unrealistic anyway).

5. **Solution-independence across flips.** The same 32 unique (rot, dist)
   tuples appear for all 4 flip matrices — the solution set is flip-independent.

## Versions Tested

- pyFAI: 2026.5.0 (latest PyPI stable release)
- ImageD11: 2.1.3 (latest PyPI stable release)
(Installed in test environment: pyFAI 2026.6.0a0, ImageD11 2.1.5)

## Code Quality

- All 37 tests pass (482 subtest variations)
- Module docstring describes the full solution space
- `mapping.md` documents mathematical derivation of all 32 solutions
- `README.md` API reference covers find_all_poni_solutions and options
- Metadata round-trips correctly through disk I/O
