# AGENTS.md — Guidelines for LLM Agents on This Project

## Single-Source-of-Truth

The **code is the spec** (`par_to_poni.py`). All tables, mappings, formulas,
and function signatures must match the code exactly. Documentation files
(mapping.md, README.md) should cite the code, not reinvent it.

Key facts that appear in multiple places (flip→orientation mapping,
azimuth factor tables, PONI formulas) must be checked against the code
after every change. Any contradiction is a bug.

## LLM Output Drift Prevention

This project was built incrementally over multiple reasoning rounds by an
LLM agent. A recurring failure mode was **output drift**: early sections
were written with a temporary conclusion, later sections superseded it,
but the agent failed to revise the earlier text. Both co-exist, creating
internal contradictions.

### Mitigation Strategy

1. **Require a final review pass.** Before accepting output, instruct the
   agent to re-read all generated files and flag every factual claim that
   appears in more than one location.

2. **Use single-source-of-truth patterns.** Express key facts (mapping tables,
   function signatures, known limitations) in exactly ONE file, referenced
   by others. Duplication invites drift. Code is the spec; markdown should
   cite it.

3. **Tests that validate documentation.** Write assertions that parse markdown
   tables and compare them against the code's data structures. E.g., extract
   the flip→orientation table from README.md and assert it matches
   `_FLIP_TO_ORIENTATION`.

4. **Separate derivation from specification.** Derivation documents (like
   mapping.md) should clearly label which sections are "final validated
   result" vs "working notes from exploration." Strike through or remove
   working notes when superseded.

5. **Round-trip the conclusions.** After generating all files, ask the LLM
   to re-read them and answer a key factual question. If it gives two
   different answers from different parts of the output, there's a bug.

6. **Versioned thinking.** Tag sections with the attempt number that
   produced them. When a new attempt supersedes an old one, deprecate or
   delete sections from the old attempt. The `story.md` attempt log helps
   with this.

## Naming Conventions

- **detector_shape** is `(slow, fast)` matching pyFAI's C-order convention
  `(shape[0], shape[1])` = `(rows, cols)`. Never swap these indices.
- **Always test on non-square detectors.** Square inputs hide axis-swap bugs.
- **Never invent a coordinate convention.** Look up the actual convention
  from the target library's source or API.

## Version History

- This file created during Round 6 replay cleanup (June 2026).
  Content extracted from `story.md` Post-Mortem section.
- **Round 9 (June 2026):** Major simplification — the solver was replaced
  with a direct formula.  `detailed_analysis/` is an archived snapshot
  of the pre-simplification code and docs.  **Ignore `detailed_analysis/`**
  when assessing correctness of the current codebase — it is historical
  reference only.  The authoritative code is `par_to_poni.py` (simplified)
  and `test_conversion.py` (core + large-tilt + force_orient3 coverage).
  `demo_ceo2_simple.ipynb` uses the current API; the old demo notebook
  is in `detailed_analysis/`.
