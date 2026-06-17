# Response to Round 2 Referee Reports

## Referee #1 (Gemini) — Happy

The scipy refactor was already implemented. No further action needed.

## Referee #2 (ChatGPT, 2nd comment) — Concerns Resolved

### Contradiction in xyz claim (Narrative A vs B)

The actual test (`test_lab_coords_match_all_orientations`) iterates over ALL 4
orientations and passes for every one (max diff < 5e-7 m). The stale "orientation
3 only" text in story.md was inherited from an earlier development iteration
and is noted as historical drift in the Round 2 section — preserved per the
story.md append-only convention.

### Detector-shape dependency

This is an intrinsic consequence of pyFAI's pixel-reordering implementation:
`_reorder_indexes_from_orientation` uses `shape[0]-1` for d1 flips and
`shape[1]-1` for d2 flips. The conversion must match these values. The
`detector_shape` parameter defaults to a square inferred from beam center
when not provided. Added discussion of this limitation to the response.

### Compensated rotation meaning

The pyFAI rotations returned by the conversion are compensated geometry
parameters, not physical tilt angles. This is now explicitly stated in the
code comments and documentation. The compensated rotations encode pyFAI's
pixel reordering and post-rotation sign flips, and the round-trip correctly
recovers the original ImageD11 tilts.

## Referee #3 (Claude, 2nd comment) — Bugs Found and Fixed

### Problem 1: d1/d2 index-range and shape swap in TestLabCoordinates

The non-square detector test had two intertwined bugs:

**Shape convention mismatch**: `detector_shape = (fast_dim, slow_dim)` per
the codebase docs, but pyFAI interprets `shape[0]` as the slow (row)
dimension in C-order. The test set `ai.detector.shape = (200, 128)` which
pyFAI interprets as 200 slow x 128 fast — swapped from the intended 200 fast
x 128 slow. Fixed: `ai.detector.shape = (SHAPE[1], SHAPE[0]) = (128, 200)`.

**max_d1/max_d2 swap in conversion**: `max_d1` was set to `shape_fast - 1`
and `max_d2` to `shape_slow - 1`. But pyFAI uses `shape[0]-1` for d1 (slow
axis) flips and `shape[1]-1` for d2 (fast axis) flips. Since pyFAI's
shape[0] = slow_dim = detector_shape[1], max_d1 should use the slow count.
Fixed: swapped assignments so max_d1 uses slow count, max_d2 uses fast count.

These bugs canceled out for square detectors (all tests except TestLabCoordinates
used 1000x1000), and the test's sampling happened to compensate for the
shape mismatch. The fix makes the conversion correct for genuinely non-square
detectors.

### Problem 2: _find_positive_equiv_from_angles silent fallback

This is working as designed, not a bug. For orientations 2 and 4, the
compensated rotation matrix differs from the tilt rotation by element-wise
sign flips — no equivalent Euler parametrization with cos(r1)·cos(r2) > 0
exists. The fallback correctly returns the compensated angles. The negative
distance they produce is handled by the round-trip: `delta = dist/(cos(r1)·cos(r2))`
recovers the positive along-beam distance. Added a comment documenting this
in the code.

### Problem 3: Transpose flips

A pyFAI limitation, not a tool limitation. pyFAI has no orientation code
for transposed axes. The conversion correctly raises ValueError for
unsupported flips.

### Problem 4: Fast/slow naming convention

The codebase convention `detector_shape = (fast_dim, slow_dim)` differs from
pyFAI's C-order convention where shape[0] = slow. Updated comments in
`par_to_poni.py` to explicitly map between the two conventions. The max_d1/
max_d2 swap fix (Problem 1) resolves the behavioral consequence of this
naming mismatch.

### Minor issues noted by Referee #3

- **read_par non-conversion of center keys**: `y_center` and `z_center` are
  pixel coordinates, not lengths, so not converting them is correct.
- **write_par float format**: Uses Python `repr()` which is unambiguous for
  float-to-string conversion. Sufficient for all ImageD11 versions.
- **_extract_rot gimbal branch**: Verified the signs are correct for the
  ZYX decomposition convention used by pyFAI.

## Upstream Issues (not fixed, identified for maintainers)

Three structural mismatches exist between pyFAI and ImageD11 that the
compensated-rotation approach correctly handles:

1. **Two-level orientation** — pyFAI applies pixel reordering (C matrix)
   and post-rotation sign flips (S matrix) at different pipeline stages.
2. **Non-standard rotation convention** — pyFAI uses left-handed Rx/Ry.
3. **0.5 pixel offset** — pyFAI centers pixels at half-integer indices;
   ImageD11 uses floating-point coordinates directly.

The conversion tool correctly handles all three. These are documentation/
design issues for the upstream projects, not fixable in this tool.

## Documentation Consistency

All 5 files cross-checked for mutual consistency. Historical sections of
story.md preserved as-is per the append-only convention; known drift items
are noted in the Round 2 section. Code changes:

- Test count: 22 (was 20 in historical story.md text)
- xyz equivalence: test covers all 4 orientations (historical text says "orientation 3 only")
- "Raw pixel indices cannot be compared" claim annotated as historical drift
- "16 pairs" claim superseded by Final Resolution section
- Sign of poni2 formula in mapping.md S11 pseudocode (- -> +)
