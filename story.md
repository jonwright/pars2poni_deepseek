# Story: Deriving the Exact par ↔ poni Mapping

## The Problem

Convert calibration parameters between ImageD11 `.par` files and pyFAI `.poni` files.
Both describe the same physical detector geometry but use different conventions for
rotations, flips, distances, and pixel centers.

## Attempts

### Attempt 1: Direct algebraic mapping from geometry_conversion.rst

Followed the documented formulas from pyFAI's `geometry_conversion.rst` and the
existing `imaged11.py` code. The formulas give exact round-trip but produce
~0.04 rad 2θ errors for non-default orientations with tilted detectors.

**Why it didn't work**: The documented mapping doesn't account for the 0.5 pixel
offset, and the orientation→flip mapping was incorrect for 2 of 4 cases.

### Attempt 2: Account for 0.5 pixel offset

Added the 0.5 pixel correction to the conversion formulas. The round-trip
improved but 2θ errors remained for non-native orientations with tilts.

**Why it didn't work**: The 0.5 offset is absorbed into the round-trip formulas
but doesn't fix the linear part of the transformation.

### Attempt 3: 4×4 affine matrix analysis

Derived both transformations as 4×4 augmented matrices. Found that pyFAI's
orientation model cannot represent arbitrary sign flips because:
- S(orient) · C_eff(orient) = (+1, +1) for ALL orientations
- ImageD11 flips require different sign patterns

**Why it was wrong**: C_eff was defined as [[c1, 0], [0, c2], [0, 0]] where
c1, c2 are the orientation's sign flips. But the EFFECTIVE linear map includes
the pixel reordering, which has opposite signs. The pixel reordering signs
CANCEL the orientation signs: C_actual = diag(1, 1) for ALL orientations.

### Attempt 4: Compensating rotation Q (180° rotations)

Added 180° rotation matrices Q = diag(±1, ±1, ±1) to compensate the coordinate
sign flips. This worked for flat detectors but not for tilted ones because
R_tilt^T·S·R_tilt is non-diagonal.

**Why it didn't work**: Q commutes with S but not with R_tilt. The product
R_tilt^T·S·R_tilt·Q involves non-diagonal entries that prevent exact matching.

### Attempt 5: C_actual = diag(1, 1) — the breakthrough

The key insight: pyFAI's orientation applies BOTH pixel reordering AND lab sign
flips. For orientation 2 (flip slow):
- Pixel reorder: d1' = max-1-d1 → linear coefficient for d1 is -1
- Lab sign flip: t1 → -t1

The effective linear coefficient (pixel→pyFAI physical): c1_eff = (-1 for
pixel reorder) × (1 for no physical sign flip on this axis... wait, it's
actually: the pixel physical coordinate p1 comes from d1' not d1).

For orientation 4 (flip fast):
- Pixel reorder: d2' = max-1-d2 → linear coefficient for d2 is -1
- Lab sign flip: t2 → -t2
- C_eff physical: [[c1, 0], [0, c2], [0, 0]] where c1=1, c2=-1
- Pixel reorder: [[a1, 0], [0, a2]] where a1=1, a2=-1
- C_actual = C_eff · pixel_reorder = [[c1·a1, 0], [0, c2·a2]] = [[1, 0], [0, 1]]

So for ALL orientations, C_actual = diag(1, 1)! The pixel reordering signs
always cancel the orientation C signs.

With C_actual = diag(1, 1), the linear system becomes:
    S · R_comp · C_actual = R_tilt · Z
    S · R_comp = R_tilt · Z

The first two columns of R_comp are:
    R_comp[:, 0] = S · R_tilt · Z[:, 0]
    R_comp[:, 1] = S · R_tilt · Z[:, 1]

where Z[:, 0] = [o11, 0, 0]^T and Z[:, 1] = [0, -o22, 0]^T.

This gives different R_comp for each (orientation, flip) pair, and ALL 16 pairs
have exact solutions where the columns are orthonormal.

### Attempt 6: Positive distance via equivalent parametrizations

For some (orientation, flip) pairs, the extracted rotation parameters give
cos(r1)·cos(r2) < 0, producing negative distance (unphysical).

The user pointed out that adding 180° to two angles gives an equivalent
rotation with positive distance. For pyFAI's ZYX Euler convention:
    (r1, r2, r3) → (r1+π, -r2, r3+π)  [same rotation, cos(r1+π)=-cos(r1)]

By searching over all equivalent parametrizations (±π, ±2π shifts, sign flips
on r2), a positive-distance parametrization exists for all 16 pairs.

### Attempt 7: Correct PONI formulas per orientation

The PONI formulas must account for the orientation-specific pixel reordering:
- Orient 3 (native): poni1 = -Δ·sin(r2) + zs·(zc+0.5), poni2 = ...(yc+0.5)
- Orient 2 (flip slow): poni1 = -Δ·sin(r2) + zs·(max-0.5-zc), poni2 = ...(yc+0.5)
- Orient 4 (flip fast): poni1 = ...(zc+0.5), poni2 = -Δ·cos(r2)·sin(r1) + ys·(max-0.5-yc)
- Orient 1 (flip both): poni1 = ...(max-0.5-zc), poni2 = ...(max-0.5-yc)

This is because the pixel at the beam center maps to a different physical index
when the orientation flips the pixel coordinates.

### Final solution

Combining all corrections:
1. C_actual = diag(1, 1) — the effective linear map after pixel reordering
2. Exact R_comp from the linear system S·R_comp = R_tilt·Z
3. Equivalent parametrization with cos(r1)·cos(r2) > 0
4. Orientation-specific PONI formulas
5. Round-trip uses the same compensated rotations in both directions

The conversion is exact to machine precision for all 4 non-transpose flip→orientation pairs.

## Post-Mortem: Documentation Review and Cleanup

After the final solution was reached, a review of the documentation files (PLAN.md,
mapping.md, and the Python code) found **stale contradictory text and mismatches**
left over from earlier attempts. The LLM had correctly derived and implemented the
final compensated-rotation solution, but failed to revise the "Known Limitations"
section and §2 algorithmic formulas to reflect the new result.

### Issues Found and Fixed

**PLAN.md:**
- The "Known Limitations" section claimed ~0.05 rad 2θ errors for non-native
  orientations — directly contradicting the "Final Solution" section in the same
  file that says exact to machine precision. Removed the stale limitation.
- §1.4, §2.1, §2.2 still used the old (incorrect) flip→orientation mapping:
  `(-1,1)→2` and `(1,1)→1`. Corrected to `(-1,1)→1` and `(1,1)→2` to match
  the code's `_FLIP_TO_ORIENTATION` dict.
- §3.1 function signatures included `wavelength_m` and `par_length_unit`
  parameters that don't exist in the actual code — those parameters live in the
  IO layer (`read_par`/`write_par`), not in the conversion functions. Fixed.
- §0 claimed `omegasign` and `fit_tolerance` default to zero; code uses 1.0 and
  0.05 respectively. Corrected.
- Added a note to §2 clarifying that the formulas show the uncompensated mapping
  and the code applies rotation compensation per §5.

**mapping.md:**
- Section 13 ("Known Limitations") repeated the same stale ~0.05 rad claim.
  Removed entirely.
- §7 had the old flip→orientation mapping table and pseudocode. Corrected to match
  the code.
- §12 had contradictory conclusions (old limitation text + final solution text
  in the same section). Replaced with a single accurate conclusion.
- Added compensation context notes to §3 (rotation mapping), §5 (distance formulas),
  and §11 (pseudocode) since those sections present the uncompensated form.

**par_to_poni.py:**
- Module docstring claimed "All 16 flip→orientation pairs" — only 4 non-transpose
  pairs are implemented. Corrected.
- `par_to_poni()` and `poni_to_par()` docstrings had empty Parameters sections.
  Filled in.
- `_compute_id11_from_pyfai()` docstring referenced a non-existent variable Q
  and described the wrong algorithm. Rewrote to match the actual implementation.
- Comment referenced `Z·Z = I` for a 3×2 matrix (invalid dimensions). Corrected
  to `Z^T·Z = I₂`.

### Root Cause: LLM Output Drift

The LLM generated the documents **incrementally** over multiple reasoning steps.
When earlier reasoning was superseded (e.g., "only orient 3 is exact" was
replaced by "all 4 are exact"), the model failed to revisit and revise the
earlier text. This is a failure mode specific to large-generation tasks where
conclusions change mid-stream.

The pattern was:
1. Early sections written with a temporary conclusion
2. Later sections written with the corrected conclusion
3. Both co-exist, creating internal contradictions

### Mitigation Strategy for LLM Agent Programming

When using LLM agents for multi-file programming tasks with evolving understanding:

1. **Require a "final review" pass.** Before accepting output, instruct the
   agent to re-read all generated files and flag every factual claim that
   appears in more than one location. Any contradiction is a bug.

2. **Use single-source-of-truth patterns.** Express key facts (mapping tables,
   function signatures, known limitations) in exactly ONE file, referenced
   by others. Duplication invites drift. For this project: the code IS the
   spec; markdown files should cite it, not reinvent it.

3. **Tests that validate documentation.** Write assertions that parse markdown
   tables and compare them against the code's data structures. E.g., extract
   the flip→orientation table from PLAN.md and assert it matches
   `_FLIP_TO_ORIENTATION`.

4. **Separate derivation from specification.** Derivation documents (like
   mapping.md) should clearly label which sections are "working notes from
   exploration" vs "final validated result." When the LLM advances to a new
   conclusion, it should strike through or remove the working notes, not
   append the new conclusion below them.

5. **Round-trip the conclusions.** After generating all files, ask the LLM
   to re-read them and answer the question: "What is the known accuracy of
   the conversion for non-native orientations?" If it gives two different
   answers from different parts of the output, there's a bug.

6. **Versioned thinking.** Tag each section with the "attempt number" that
   produced it. When a new attempt supersedes an old one, deprecate or
   delete sections from the old attempt. The story.md attempt log helps
   with this, but the derivative documents (PLAN.md, mapping.md) were not
   cleaned up accordingly.

## LLM Attribution

Model used: DeepSeek V4 Pro on medium thinking (generation).
Review cleanup: DeepSeek V4 Pro, same session.
Total original cost: $0.74.
