# par_to_poni.py

Convert ImageD11 `.par` ↔ pyFAI `.poni` geometry files.

## Get the file

```bash
curl -O https://raw.githubusercontent.com/jonwright/pars2poni_deepseek/main/par_to_poni.py
```

Requires `numpy`, `scipy`.  Put `par_to_poni.py` next to your script and import it.

## Convert a file

```python
import par_to_poni as pp

shape = (2162, 2068)   # Eiger 4M: (slow, fast)

# par → poni
par = pp.read_par("geometry.par")
poni = pp.par_to_poni(par, detector_shape=shape)
pp.write_poni(poni, "geometry.poni")

# poni → par
poni = pp.read_poni("geometry.poni")
par = pp.poni_to_par(poni, detector_shape=shape)
pp.write_par(par, "geometry.par")

## Azimuth mapping (chi ↔ eta)

PyFAI chi and ImageD11 eta are related by **orientation-dependent** formulas.
The mapping is **not** simply `chi = 90° − eta` for all orientations.

```python
chi_rad = ...                               # from pyFAI
eta = pp.chi_to_eta(chi_rad, orientation=3) # → ImageD11 eta (radians)
chi = pp.eta_to_chi(eta, orientation=3)     # → pyFAI chi (radians)
```

Orientation can be an `int` (1–4), a par dict, or a poni dict.
See the function docstrings for the full per-orientation mapping table.

| Orient | chi = | sin(chi) | cos(chi) |
|--------|-------|----------|----------|
| 3 | 90° − eta | +cos(eta) | +sin(eta) |
| 2 | eta − 90° | −cos(eta) | +sin(eta) |
| 4 | eta + 90° | +cos(eta) | −sin(eta) |
| 1 | 270° − eta | −cos(eta) | −sin(eta) |

## Status

All 4 non-transpose flip orientations (1, 2, 3, 4) match exactly:
- 2θ at machine precision (10⁻¹⁶ rad)
- Azimuth with orientation-dependent simple mapping (above)
- Lab coordinates with per-orientation mirror reflections
- Round-trip par ↔ poni exact

Transpose flips (`o12, o21 ≠ 0`) are not supported.
Spatial distortion is not handled.

---

*Author: DeepSeek V4 Pro (opencode), June 2026, for Jon Wright*
