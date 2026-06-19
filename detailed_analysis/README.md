# detailed_analysis/ — Archived pre-simplification material

This folder preserves the full solver, derivation, and development history
from before the Round 9 simplification (June 2026).

## Contents

| File | Description |
|------|-------------|
| `par_to_poni_full.py` | Original solver — `find_all_poni_solutions()`, mirror matrices, rotation compensation, 32-solution enumeration |
| `mapping.md` | Complete mathematical derivation of the affine pipeline equation |
| `story.md` | Development history, attempt log, and narrative |
| `task.md` | Original task description |

## Why archived

The conversion was found to reduce to a simple direct formula:
`rot = (−tz, ty, tx)` with an orientation from a corrected
flip→orientation lookup table.  See `../par_to_poni.py` for the
current (simplified) implementation and `../README.md` for usage.

All test cases that exercise the full solver are in `../test_conversion.py`
(classes that need `find_all_poni_solutions` import from `par_to_poni_full`).
- `demo_ceo2.ipynb` — full 32-solution demo notebook (imports find_all_poni_solutions from `par_to_poni_full.py`)
