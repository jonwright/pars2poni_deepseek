"""
par_to_poni.py — Convert between ImageD11 .par and pyFAI .poni geometry parameters.

Coordinate systems (from library source code):
- pyFAI: axis1=y_up(slow), axis2=x_starboard(fast), axis3=z_downstream
  (pyFAI/geometry/core.py:2668-2670).  Rotation: Rz(rot3)·Ry(−rot2)·Rx(−rot1)
  (core.py:2690-2702).  R1,R2 are right-handed with negated angles; pyFAI's own
  comments label them "left-handed" referring to the sign convention.
- ImageD11: X=downstream, Y=port, Z=up(slow) (transform.py:133-141).
  Rotation: Rx(tilt_x)·Ry(tilt_y)·Rz(tilt_z) (transform.py:57-81, all r.h.).

Frame transform: t_ID11 = G·t_pyFAI, G = [[0,0,1],[0,-1,0],[1,0,0]].

Orientations from pyFAI internals:
- Pixel reordering: _reorder_indexes_from_orientation  (_common.py:657)
- Sign flips: f_t1 / f_t2                            (_geometry.pyx:68-105)

Equating the full affine pipeline gives exact closed-form solutions for all
4×4 = 16 valid (ImageD11 flip, pyFAI orientation) mappings.  The equation

    S(orient) · R · C(orient) = M · R_tilt · Z(flip)

is solved with each of the 4 distinct mirror matrices {M1, M2, M3=I, M4}.
Per (orient, mirror) group there are two ZYX Euler-angle representations
(the β-solution pair), giving 32 unique solutions per flip.  Transpose
flips (o12,o21≠0) are not supported.

Azimuth mapping — the pyFAI chi and ImageD11 eta angles are related by
mirror-dependent formulas; see _azimuth_factors().

Dependencies: numpy, scipy (for Rotation). All internal units are meters for
lengths and meters for wavelength.

Usage:
    import par_to_poni as pp

    par = pp.read_par("geometry.par")
    poni = pp.par_to_poni(par, detector_shape=(2162, 2068))
    pp.write_poni(poni, "geometry.poni")

    poni = pp.read_poni("geometry.poni")
    par = pp.poni_to_par(poni, detector_shape=(2162, 2068))
    pp.write_par(par, "geometry.par")

    # Convert azimuth angles between the two programs:
    chi_rad = 1.2                     # from pyFAI
    eta_rad = pp.chi_to_eta(chi_rad, orientation=3)
    chi_rad = pp.eta_to_chi(eta_rad, orientation=3)

    # Orientation can come from a par or poni dict:
    eta_rad = pp.chi_to_eta(chi_rad, par)
    chi_rad = pp.eta_to_chi(eta_rad, poni)

    # Find all valid solutions (4 orientations × 4 mirrors × 2 Euler reps):
    solutions = pp.find_all_poni_solutions(par, detector_shape=(2162, 2068))
    for s in solutions:
        print(s["flip_label"], s["orient_tried"], s["mirror_source"],
              s["dist_positive"], s["chi_eta_exact"])
"""

import json
import math
from math import cos, sin, tan, atan2, pi, sqrt

import numpy as np
from scipy.spatial.transform import Rotation as ScipyRotation


# ---------------------------------------------------------------------------
# Unit handling
# ---------------------------------------------------------------------------

_LENGTH_UNIT_FACTORS = {"um": 1e6, "mm": 1e3, "m": 1.0}
_WAVELENGTH_A_PER_M = 1e10


def _parse_length_unit(unit):
    u = str(unit).lower().replace("µm", "um").replace("μ", "u")
    if u in _LENGTH_UNIT_FACTORS:
        return u, _LENGTH_UNIT_FACTORS[u]
    raise ValueError(f"Unsupported length unit: {unit!r}. Use um, mm, or m.")


# ---------------------------------------------------------------------------
# Flip / Orientation mapping (non-transpose only)
# ---------------------------------------------------------------------------

_FLIP_TO_ORIENTATION = {
    (1, 0, 0, -1): 3,
    (-1, 0, 0, 1): 1,
    (-1, 0, 0, -1): 4,
    (1, 0, 0, 1): 2,
}

_ORIENTATION_TO_FLIP = {v: k for k, v in _FLIP_TO_ORIENTATION.items()}


def flip_to_orientation(o11, o12, o21, o22):
    key = (int(o11), int(o12), int(o21), int(o22))
    if key not in _FLIP_TO_ORIENTATION:
        raise ValueError(
            f"Unsupported flip matrix [{o11},{o12},{o21},{o22}]. "
            f"Only non-transpose flips are supported."
        )
    return _FLIP_TO_ORIENTATION[key]


def orientation_to_flip(orientation):
    return _ORIENTATION_TO_FLIP[orientation]


# ---------------------------------------------------------------------------
# Azimuth angle conversion (chi ↔ eta)
# ---------------------------------------------------------------------------

# Per-orientation mapping: sin(chi)/cos(chi) expressed as (s0·cos(eta), s1·sin(eta))
# Derived from equating the pyFAI coordinate frame (with per-orientation mirror)
# against the ImageD11 coordinate frame.  Verified by test_conversion.py.
_CHI_ETA_SIN_COS_FACTORS = {
    3: (+1, +1),   # chi =  90° − eta    sin(chi) = +cos(eta)  cos(chi) = +sin(eta)
    2: (-1, +1),   # chi =  eta − 90°    sin(chi) = −cos(eta)  cos(chi) = +sin(eta)
    4: (+1, -1),   # chi =  eta + 90°    sin(chi) = +cos(eta)  cos(chi) = −sin(eta)
    1: (-1, -1),   # chi = 270° − eta    sin(chi) = −cos(eta)  cos(chi) = −sin(eta)
}


def _azimuth_factors(mirror_orient):
    """Return the (sin_factor, cos_factor) for the azimuth mapping
    chi ↔ eta given a mirror orientation.

    The relationship is:
        sin(chi) = sin_factor · cos(eta)
        cos(chi) = cos_factor · sin(eta)

    Derivation:  From t_pyFAI = M · t_id (the mirror-transformed
    ID11 coordinates in the pyFAI lab frame), chi = atan2(m0·t0, m1·t1)
    and eta = atan2(t1, t0), where M = diag(m0, m1, 1) = M(mirror_orient).
    Therefore sin(chi) = m0·cos(eta), cos(chi) = m1·sin(eta).
    The factors are the first two diagonal entries of the mirror matrix.
    """
    M_diag_axes = {
        3: (+1, +1),   # M(3) = I
        2: (-1, +1),   # M(2) = diag(-1, 1, 1)
        4: (+1, -1),   # M(4) = diag(1, -1, 1)
        1: (-1, -1),   # M(1) = diag(-1, -1, 1)
    }
    return M_diag_axes.get(mirror_orient, (+1, +1))


def chi_to_eta(chi_rad, orientation):
    """Convert pyFAI azimuthal angle chi → ImageD11 azimuthal angle eta.

    Parameters
    ----------
    chi_rad : float or array-like
        Azimuthal angle from pyFAI (radians).  Defined as
        ``atan2(t1, t2)`` in pyFAI lab coordinates; origin at +x
        (starboard), positive CCW toward +y (up).
    orientation : int, dict
        pyFAI orientation (1–4), or a par dict (with o11..o22 keys),
        or a poni dict (with an 'orientation' key).

    Returns
    -------
    float or ndarray
        ImageD11 eta angle (radians).  Defined as ``atan2(−t_y, t_z)``
        in ID11 lab coordinates; origin at +z (up), positive CW
        facing downstream.

    Notes
    -----
    The mapping depends on orientation because pyFAI's pixel reordering
    and sign flips change the effective azimuth origin.  All results are
    equivalent modulo 2π — use sin/cos comparisons to avoid
    wrap-around ambiguity.

    ┌─────────┬──────────────────────┬─────────────────────────┐
    │ orient  │ chi = f(eta)         │ sin(chi), cos(chi)      │
    ├─────────┼──────────────────────┼─────────────────────────┤
    │ 3       │ chi =  90° − eta     │ (+cos(eta), +sin(eta))  │
    │ 2       │ chi =  eta − 90°     │ (−cos(eta), +sin(eta))  │
    │ 4       │ chi =  eta + 90°     │ (+cos(eta), −sin(eta))  │
    │ 1       │ chi = 270° − eta     │ (−cos(eta), −sin(eta))  │
    └─────────┴──────────────────────┴─────────────────────────┘

    Examples
    --------
    >>> import numpy as np
    >>> eta = chi_to_eta(np.radians(45), orientation=3)
    >>> np.degrees(eta)
    45.0
    >>> eta = chi_to_eta(np.radians(120), orientation=2)
    >>> np.degrees(eta) % 360
    30.0
    """
    import numpy as np
    orientation = _extract_orientation_from_arg(orientation)
    s0, s1 = _CHI_ETA_SIN_COS_FACTORS[orientation]
    return np.arctan2(s1 * np.cos(chi_rad), s0 * np.sin(chi_rad))


def eta_to_chi(eta_rad, orientation):
    """Convert ImageD11 azimuthal angle eta → pyFAI azimuthal angle chi.

    Parameters
    ----------
    eta_rad : float or array-like
        Azimuthal angle from ImageD11 (radians).  Defined as
        ``atan2(−t_y, t_z)`` in ID11 lab coordinates; origin at +z
        (up), positive CW facing downstream.
    orientation : int, dict
        pyFAI orientation (1–4), or a par dict (with o11..o22 keys),
        or a poni dict (with an 'orientation' key).

    Returns
    -------
    float or ndarray
        pyFAI chi angle (radians).  Defined as ``atan2(t1, t2)`` in
        pyFAI lab coordinates; origin at +x (starboard), positive
        CCW toward +y (up).

    Notes
    -----
    See `chi_to_eta` for the per-orientation mapping table.
    Inverse of `chi_to_eta`: ``eta_to_chi(chi_to_eta(c, o), o)``
    recovers the original chi (modulo 2π).

    Examples
    --------
    >>> import numpy as np
    >>> chi = eta_to_chi(np.radians(45), orientation=3)
    >>> np.degrees(chi)
    45.0
    >>> chi = eta_to_chi(np.radians(30), orientation=4)
    >>> np.degrees(chi) % 360
    120.0
    """
    import numpy as np
    orientation = _extract_orientation_from_arg(orientation)
    s0, s1 = _CHI_ETA_SIN_COS_FACTORS[orientation]
    return np.arctan2(s0 * np.cos(eta_rad), s1 * np.sin(eta_rad))


def _extract_orientation_from_arg(arg):
    """Extract pyFAI orientation (1–4) from a par dict, poni dict, or int.

    Parameters
    ----------
    arg : dict or int
        If dict: must contain either 'orientation' key (poni) or
        'o11','o12','o21','o22' keys (par).
        If int: returned directly.

    Returns
    -------
    int
        pyFAI orientation (1, 2, 3, or 4).
    """
    if isinstance(arg, dict):
        if "orientation" in arg:
            return int(arg["orientation"])
        elif "o11" in arg and "o22" in arg:
            o11 = int(arg.get("o11", 1))
            o12 = int(arg.get("o12", 0))
            o21 = int(arg.get("o21", 0))
            o22 = int(arg.get("o22", -1))
            return flip_to_orientation(o11, o12, o21, o22)
        else:
            raise ValueError(
                "Dict must contain 'orientation' (poni) or "
                "'o11','o12','o21','o22' (par)")
    return int(arg)


# ---------------------------------------------------------------------------
# Rotation matrix utilities
# ---------------------------------------------------------------------------

def _pyfai_rotation_matrix(rot1, rot2, rot3):
    """pyFAI rotation matrix (pyFAI/geometry/core.py:2690-2702).

    pyFAI builds: rot1=Rx(−rot1), rot2=Ry(−rot2), rot3=Rz(rot3)
    then R = rot3·rot2·rot1 = Rz(rot3)·Ry(−rot2)·Rx(−rot1).

    pyFAI's source comments label R1,R2 "left-handed" but all three are
    standard right-handed rotation matrices — the angles for R1,R2 are
    simply negated.  Net result is intrinsic ZYX with [rot3, −rot2, −rot1].

    Matches pyFAI's actual rotation_matrix() to machine precision
    (verified by test_pyfai_rotation_matrix_matches_actual).

    Returns a tuple-of-tuples for backward compatibility.
    """
    R = ScipyRotation.from_euler('ZYX', [rot3, -rot2, -rot1]).as_matrix()
    return ((R[0, 0], R[0, 1], R[0, 2]),
            (R[1, 0], R[1, 1], R[1, 2]),
            (R[2, 0], R[2, 1], R[2, 2]))


def _extract_rot(R):
    """Extract rot1,rot2,rot3 from a pyFAI rotation matrix.

    R = Rz(rot3).Ry_left(rot2).Rx_left(rot1).
    Entries: R[2,0]=sin(rot2), R[2,1]=-cos(rot2).sin(rot1),
             R[2,2]=cos(rot2).cos(rot1),
             R[1,0]=sin(rot3).cos(rot2), R[0,0]=cos(rot3).cos(rot2).
    """
    r00, r01, r02 = R[0]
    r10, r11, r12 = R[1]
    r20, r21, r22 = R[2]

    if abs(r20) < 0.999999:
        rot2 = atan2(r20, sqrt(r21 * r21 + r22 * r22))
        rot1 = atan2(-r21, r22)
        rot3 = atan2(r10, r00)
    else:
        rot3 = 0.0
        if r20 > 0:
            rot2 = pi / 2
            rot1 = atan2(-r01, -r02)
        else:
            rot2 = -pi / 2
            rot1 = atan2(r01, r02)
    return rot1, rot2, rot3


def _find_positive_equiv_from_angles(rot1, rot2, rot3):
    """Find equivalent Euler angles with cos(rot1)*cos(rot2) > 0.

    For pyFAI's ZYX convention, equivalent parametrizations include:
      (rot1+π, -rot2, rot3+π)  and  (rot1-π, -rot2, rot3-π)
    Searches over ±π offsets on all three angles and sign flip on rot2.

    With the mirror-matrix compensation (see _get_mirror_matrix), all
    orientations now produce R[2,2] > 0 in the raw decomposition.  This
    function remains as a safety net for edge cases.
    """
    R_target = _pyfai_rotation_matrix(rot1, rot2, rot3)
    best = None
    for d1 in (0, pi, -pi):
        for d2 in (0, pi, -pi, 2 * pi, -2 * pi):
            for d3 in (0, pi, -pi):
                for s2 in (1, -1):
                    rt1, rt2, rt3 = rot1 + d1, s2 * rot2 + d2, rot3 + d3
                    if abs(rt1) > 10 or abs(rt2) > 10 or abs(rt3) > 10:
                        continue
                    Rt = _pyfai_rotation_matrix(rt1, rt2, rt3)
                    maxdiff = max(abs(Rt[i][j] - R_target[i][j])
                                  for i in range(3) for j in range(3))
                    if maxdiff < 1e-8:
                        dc = cos(rt1) * cos(rt2)
                        if dc > 0:
                            if best is None or abs(rt1)+abs(rt2) < abs(best[0])+abs(best[1]):
                                best = (rt1, rt2, rt3)
    return best


def _find_all_rot_equivs(rot1, rot2, rot3):
    """Find all equivalent (rot1,rot2,rot3) that produce the same rotation matrix.

    For pyFAI's ZYX convention the key equivalence is:
      Rz(rot3)·Ry(−rot2)·Rx(−rot1)  =  Rz(rot3+π)·Ry(rot2)·Rx(−rot1+π)
    which in terms of the pyFAI angle tuple is (rot1−π, −rot2, rot3+π).

    Searches the same grid as _find_positive_equiv_from_angles but collects
    *all* matches, both positive and negative distance.  Results are
    normalised to (−π, π] and deduplicated.
    """
    R_target = _pyfai_rotation_matrix(rot1, rot2, rot3)
    seen = set()
    results = []
    for d1 in (0, pi, -pi):
        for d2 in (0, pi, -pi, 2 * pi, -2 * pi):
            for d3 in (0, pi, -pi):
                for s2 in (1, -1):
                    rt1, rt2, rt3 = rot1 + d1, s2 * rot2 + d2, rot3 + d3
                    if abs(rt1) > 10 or abs(rt2) > 10 or abs(rt3) > 10:
                        continue
                    Rt = _pyfai_rotation_matrix(rt1, rt2, rt3)
                    maxdiff = max(abs(Rt[i][j] - R_target[i][j])
                                  for i in range(3) for j in range(3))
                    if maxdiff < 1e-8:
                        n1 = atan2(sin(rt1), cos(rt1))
                        n2 = atan2(sin(rt2), cos(rt2))
                        n3 = atan2(sin(rt3), cos(rt3))
                        key = (round(n1, 8), round(n2, 8), round(n3, 8))
                        if key not in seen:
                            seen.add(key)
                            results.append((n1, n2, n3))
    return results


def _get_mirror_matrix(orient):
    """Return the mirror matrix for coordinate-frame relaxation.

    For each non-native orientation, pyFAI flips specific pixel axes and
    lab-coordinate signs.  A matching mirror in the rotation constraint
    aligns the effective coordinate system with the detector's fast/slow
    axes while keeping distance positive:

      orient 2 (flip slow / pyFAI axis 1):  diag(-1,  1,  1)
      orient 4 (flip fast / pyFAI axis 2):  diag( 1, -1,  1)
      orient 1 (flip both):                 diag(-1, -1,  1)
      orient 3 (native, no flip):           identity

    Each mirror is self-inverse.  The mirror relaxes xyz coordinate
    matching in the ID11 frame but preserves 2θ and azimuth exactly.
    See chi_to_eta / eta_to_chi for the per-orientation azimuth mapping.
    """
    _m = {
        3: np.eye(3),
        2: np.diag([-1.0, 1.0, 1.0]),
        4: np.diag([1.0, -1.0, 1.0]),
        1: np.diag([-1.0, -1.0, 1.0]),
    }
    return _m[orient]


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------

def _compute_compensated_rotation(o11, o22, orient, r1_std, r2_std, r3_std,
                                  mirror_M=None, find_positive_equiv=True):
    """Compute compensated pyFAI rotation for a given (flip, orientation) pair.

    Derivation: equating the full pyFAI pipeline against the ID11 pipeline.
    For each orientation, pyFAI applies:
      - Pixel reordering: C = diag(c1, c2)  pre-rotation
      - Rotation: R
      - Sign flips:  S = diag(s1, s2, 1)  post-rotation

    ID11 applies the flip matrix Z = diag(o11, -o22) pre-rotation (in the
    pyFAI lab frame after G transformation).  The linear constraint is:

        S . R_comp . C = M . R_tilt . Z

    where M is the per-orientation mirror matrix (see _get_mirror_matrix)
    that relaxes xyz coordinate matching while preserving 2θ and azimuth.

    Solving:  R_comp[:,0] = S . M . R_tilt[:,0] . (o11 / c1)
              R_comp[:,1] = S . M . R_tilt[:,1] . (-o22 / c2)

    Parameters
    ----------
    find_positive_equiv : bool, optional
        If True (default) apply _find_positive_equiv_from_angles safety net.
        If False, return the raw Euler angles directly (caller is
        responsible for finding equivalent representations).
    """
    S_diag = {3: (1, 1, 1), 2: (-1, 1, 1), 4: (1, -1, 1), 1: (-1, -1, 1)}[orient]
    c1 = -1.0 if orient in (2, 1) else 1.0
    c2 = -1.0 if orient in (4, 1) else 1.0

    R_tilt = np.array(_pyfai_rotation_matrix(r1_std, r2_std, r3_std))
    if mirror_M is not None:
        R_tilt = mirror_M @ R_tilt

    r_c0 = np.array([S_diag[0] * R_tilt[0, 0] * (o11 / c1),
                     S_diag[1] * R_tilt[1, 0] * (o11 / c1),
                     S_diag[2] * R_tilt[2, 0] * (o11 / c1)])
    r_c1 = np.array([S_diag[0] * R_tilt[0, 1] * (-o22 / c2),
                     S_diag[1] * R_tilt[1, 1] * (-o22 / c2),
                     S_diag[2] * R_tilt[2, 1] * (-o22 / c2)])

    r_c2 = np.cross(r_c0, r_c1)
    if np.linalg.det(np.column_stack([r_c0, r_c1, r_c2])) < 0:
        r_c2 = -r_c2

    R_comp = np.column_stack([r_c0, r_c1, r_c2])
    rot_s = ScipyRotation.from_matrix(R_comp)
    angles = rot_s.as_euler('ZYX')
    rot3_c, rot2_c, rot1_c = angles[0], -angles[1], -angles[2]

    if find_positive_equiv:
        result = _find_positive_equiv_from_angles(rot1_c, rot2_c, rot3_c)
        if result is None:
            result = (rot1_c, rot2_c, rot3_c)
        return result
    return (rot1_c, rot2_c, rot3_c)


def _compute_id11_from_pyfai(rot1, rot2, rot3, orient, mirror_M=None,
                              o11=None, o22=None):
    """Recover ID11 tilt rotation from compensated pyFAI params.

    Reverse of _compute_compensated_rotation. From the forward equation:
      S . R_comp . C = M . R_tilt . Z
    reverse:
      R_tilt[:,0] = M^{-1} . S . R_comp[:,0] . (c1 / o11)
      R_tilt[:,1] = M^{-1} . S . R_comp[:,1] . (c2 / (-o22))

    Parameters
    ----------
    o11, o22 : int, optional
        Override the ImageD11 flip values (from the original par).
        If None, defaults to ``orientation_to_flip(orient)``.
    """
    S_diag = {3: (1, 1, 1), 2: (-1, 1, 1), 4: (1, -1, 1), 1: (-1, -1, 1)}[orient]
    c1 = -1.0 if orient in (2, 1) else 1.0
    c2 = -1.0 if orient in (4, 1) else 1.0

    if o11 is None or o22 is None:
        _o11, _o12, _o21, _o22 = orientation_to_flip(orient)
        if o11 is None:
            o11 = _o11
        if o22 is None:
            o22 = _o22

    R_comp = np.array(_pyfai_rotation_matrix(rot1, rot2, rot3))

    rt_c0 = np.array([S_diag[0] * R_comp[0, 0] * (c1 / o11),
                      S_diag[1] * R_comp[1, 0] * (c1 / o11),
                      S_diag[2] * R_comp[2, 0] * (c1 / o11)])
    rt_c1 = np.array([S_diag[0] * R_comp[0, 1] * (c2 / (-o22)),
                      S_diag[1] * R_comp[1, 1] * (c2 / (-o22)),
                      S_diag[2] * R_comp[2, 1] * (c2 / (-o22))])

    if mirror_M is not None:
        rt_c0 = mirror_M @ rt_c0
        rt_c1 = mirror_M @ rt_c1

    rt_c2 = np.cross(rt_c0, rt_c1)
    if np.linalg.det(np.column_stack([rt_c0, rt_c1, rt_c2])) < 0:
        rt_c2 = -rt_c2

    R_tilt = (
        (rt_c0[0], rt_c1[0], rt_c2[0]),
        (rt_c0[1], rt_c1[1], rt_c2[1]),
        (rt_c0[2], rt_c1[2], rt_c2[2]),
    )
    return _extract_rot(R_tilt)


def _build_poni_from_compensated_rots(par, trial_orient, rot1, rot2, rot3,
                                       detector_shape, use_mirror,
                                       forward_mirror_orient=None,
                                       forward_o11=None, forward_o22=None):
    """Build a full poni dict from compensated rotations and a par dict.

    This factors out the poni-building logic so that the solution finder
    can generate multiple candidate poni dicts from different (rot1,rot2,rot3)
    angle triples.

    Parameters
    ----------
    trial_orient : int
        The pyFAI orientation (1-4) used for PONI formulas and stored as
        ``orientation`` in the poni dict.  May differ from the canonical
        orientation for the par's flip.
    forward_o11, forward_o22 : int, optional
        The original ImageD11 flip values for this solution.  Used during
        reverse conversion (poni_to_par) to recover the correct flip.
    forward_mirror_orient : int, optional
        Which orientation's mirror matrix was used (1-4), or None for
        identity mirror.
    """
    distance = float(par["distance"])
    yc = float(par["y_center"])
    zc = float(par["z_center"])
    ys = float(par["y_size"])
    zs = float(par["z_size"])
    wl_m = float(par.get("wavelength", 0.0))
    delta = distance

    if detector_shape is None:
        shape_fast = max(int(2 * (yc if yc > 0 else -yc) + 1), 2)
        shape_slow = max(int(2 * (zc if zc > 0 else -zc) + 1), 2)
        if shape_fast < 2 or shape_fast > 100000:
            shape_fast = 2
        if shape_slow < 2 or shape_slow > 100000:
            shape_slow = 2
        detector_shape = (shape_slow, shape_fast)

    shape_slow, shape_fast = int(detector_shape[0]), int(detector_shape[1])
    max_d1 = shape_slow - 1.0
    max_d2 = shape_fast - 1.0

    dist = delta * cos(rot2) * cos(rot1)

    if trial_orient in (2, 1):
        poni1 = -delta * sin(rot2) + zs * (max_d1 - zc + 0.5)
    else:
        poni1 = -delta * sin(rot2) + zs * (zc + 0.5)

    if trial_orient in (4, 1):
        poni2 = delta * cos(rot2) * sin(rot1) + ys * (max_d2 - yc + 0.5)
    else:
        poni2 = delta * cos(rot2) * sin(rot1) + ys * (yc + 0.5)

    result = {
        "dist": dist,
        "poni1": poni1,
        "poni2": poni2,
        "rot1": rot1,
        "rot2": rot2,
        "rot3": rot3,
        "pixel1": zs,
        "pixel2": ys,
        "wavelength": wl_m,
        "orientation": trial_orient,
        "_mirror_used": use_mirror,
    }
    if forward_o11 is not None and forward_o22 is not None:
        result["_forward_o11"] = int(forward_o11)
        result["_forward_o22"] = int(forward_o22)
    if forward_mirror_orient is not None:
        result["_mirror_orient"] = int(forward_mirror_orient)
    return result


def _deduplicate_solutions(solutions):
    """Remove duplicate solutions that have the same rot angles and distance."""
    seen = set()
    unique = []
    for sol in solutions:
        p = sol["poni"]
        key = (round(p["rot1"], 8), round(p["rot2"], 8), round(p["rot3"], 8),
               round(p["dist"], 12))
        if key not in seen:
            seen.add(key)
            unique.append(sol)
    return unique


def find_all_poni_solutions(par, detector_shape=None,
                            include_backscattering=False):
    """Find all valid pyFAI poni solutions for a given ImageD11 par dict.

    For the flip carried by the ``par`` argument, enumerates **32 solutions**
    — 4 trial pyFAI orientations × 4 mirror matrices {M1, M2, M3=I, M4}
    × 2 ZYX Euler-angle representations per rotation matrix.  Across all 4
    possible flip matrices the solution set is flip-independent (the same
    32 unique (rot, dist) tuples appear for every flip).

    When ``include_backscattering=True``, seed rotations with ±π offsets
    on rot1 or rot2 are additionally explored to discover representations
    suitable for backscattering geometry (negative ImageD11 distance).

    Parameters
    ----------
    par : dict
        ImageD11 par dict.
    detector_shape : (slow_dim, fast_dim) tuple, optional
        Detector shape.
    include_backscattering : bool, optional
        If True, also search using π-offset seed rotations (for
        backscattering geometry).

    Returns
    -------
    list of dict
        Each element has keys:
        ``poni`` — the poni dict,
        ``use_mirror`` — whether a non-identity mirror was used,
        ``dist_positive`` — whether orthogonal distance is positive,
        ``chi_eta_exact`` — True when mirror matches trial orientation
            (canonical azimuth mapping: S·M=I in the coordinate frame),
        ``rot_magnitude`` — |rot1| + |rot2| + |rot3| (for ranking),
        ``flip_label`` — which ImageD11 flip (e.g. "F_o3"),
        ``orient_tried`` — which pyFAI orientation was tried,
        ``mirror_source`` — which mirror was used ("M3", "M1", "M2", "M4"),
        ``is_canonical`` — whether orientation matches canonical for this flip.
        List is sorted best-first (positive distance preferred, then
        smallest rotation magnitude).
    """
    tx = float(par.get("tilt_x", 0.0))
    ty = float(par.get("tilt_y", 0.0))
    tz = float(par.get("tilt_z", 0.0))
    o11 = int(par.get("o11", 1))
    o12 = int(par.get("o12", 0))
    o21 = int(par.get("o21", 0))
    o22 = int(par.get("o22", -1))

    r1 = -tz
    r2 = ty
    r3 = tx

    seeds = [(r1, r2, r3)]
    if include_backscattering:
        seeds.extend([
            (r1, r2 + pi, r3),
            (r1, r2 - pi, r3),
            (r1 + pi, r2, r3),
            (r1 - pi, r2, r3),
        ])

    if detector_shape is None:
        shape_fast = max(int(2 * float(par.get("y_center", 500)) + 1), 2)
        shape_slow = max(int(2 * float(par.get("z_center", 500)) + 1), 2)
        detector_shape = (shape_slow, shape_fast)

    _ALL_ORIENTS = [1, 2, 3, 4]
    _MIRROR_CHOICES = [
        (3, "M3", True),                # M(3) = I = diag(1,1,1)
        (1, "M1", True),                # M(1) = diag(-1,-1,1)
        (2, "M2", True),                # M(2) = diag(-1,1,1)
        (4, "M4", True),                # M(4) = diag(1,-1,1)
    ]

    flip_key = (o11, o12, o21, o22)
    if flip_key == (1, 0, 0, -1):
        flip_label = "F_o3"
    elif flip_key == (-1, 0, 0, 1):
        flip_label = "F_o1"
    elif flip_key == (-1, 0, 0, -1):
        flip_label = "F_o4"
    elif flip_key == (1, 0, 0, 1):
        flip_label = "F_o2"
    else:
        flip_label = f"F_({o11},{o12},{o21},{o22})"

    canonical_orient = flip_to_orientation(o11, o12, o21, o22)

    solutions = []
    for trial_orient in _ALL_ORIENTS:
        for mirror_orient, mirror_name, use_mirror in _MIRROR_CHOICES:
            mirror_M = _get_mirror_matrix(mirror_orient)

            for sr1, sr2, sr3 in seeds:
                rot1_base, rot2_base, rot3_base = _compute_compensated_rotation(
                    o11, o22, trial_orient, sr1, sr2, sr3,
                    mirror_M=mirror_M, find_positive_equiv=False)

                for erot1, erot2, erot3 in _find_all_rot_equivs(rot1_base, rot2_base, rot3_base):
                    chi_eta_exact = (mirror_orient == trial_orient)

                    poni = _build_poni_from_compensated_rots(
                        par, trial_orient, erot1, erot2, erot3,
                        detector_shape, use_mirror,
                        forward_mirror_orient=mirror_orient,
                        forward_o11=o11, forward_o22=o22)

                    solutions.append({
                        "poni": poni,
                        "use_mirror": use_mirror,
                        "dist_positive": poni["dist"] > 0,
                        "chi_eta_exact": chi_eta_exact,
                        "rot_magnitude": abs(erot1) + abs(erot2) + abs(erot3),
                        "flip_label": flip_label,
                        "orient_tried": trial_orient,
                        "mirror_source": mirror_name,
                        "is_canonical": trial_orient == canonical_orient,
                    })

    unique = _deduplicate_solutions(solutions)
    unique.sort(key=lambda s: (
        0 if s["dist_positive"] else 1,
        s["rot_magnitude"],
        0 if s["is_canonical"] else 1,
    ))
    return unique


def par_to_poni(par, detector_shape=None,
                prefer_positive_distance=True,
                exact_chi=False):
    """Convert ImageD11 .par parameters -> pyFAI .poni parameters.

    Parameters
    ----------
    par : dict
        Keys: distance, y_center, z_center, y_size, z_size,
        tilt_x, tilt_y, tilt_z, o11, o12, o21, o22, wavelength.
        All lengths in meters internally, wavelength in meters.
    detector_shape : (slow_dim, fast_dim) tuple, optional
        Detector pixel dimensions matching pyFAI's C-order shape
        convention: shape[0] = slow (height/rows), shape[1] = fast
        (width/columns).  Required for non-native orientations
        (2 and 4).  Defaults to square inferred from beam center.
    prefer_positive_distance : bool, optional
        If True (default) prefer solutions with positive orthogonal
        distance (the current behaviour).  If False, also consider
        solutions with negative distance.
    exact_chi : bool, optional
        If True, prefer the no-mirror solution where chi = 90°−eta
        for *all* orientations.  Within that family the solution with
        positive distance is preferred when available.

    Returns
    -------
    dict
        Keys: dist, poni1, poni2, rot1, rot2, rot3,
        pixel1, pixel2, wavelength, orientation, _mirror_used.
    """
    solutions = find_all_poni_solutions(par, detector_shape=detector_shape)

    if not solutions:
        raise RuntimeError("No valid poni solution found for the given par parameters.")

    if exact_chi:
        candidates = [s for s in solutions if s["chi_eta_exact"]]
        if not candidates:
            candidates = solutions
    else:
        candidates = solutions

    if prefer_positive_distance:
        pos = [s for s in candidates if s["dist_positive"]]
        if pos:
            candidates = pos
        else:
            bs_solutions = find_all_poni_solutions(
                par, detector_shape=detector_shape, include_backscattering=True)
            if exact_chi:
                bs_candidates = [s for s in bs_solutions if s["chi_eta_exact"]]
                if bs_candidates:
                    bs_solutions = bs_candidates
            pos = [s for s in bs_solutions if s["dist_positive"]]
            if pos:
                candidates = pos
            else:
                candidates = bs_solutions if bs_solutions else candidates

    best = min(candidates, key=lambda s: s["rot_magnitude"])
    return best["poni"]


def poni_to_par(poni, detector_shape=None):
    """Convert pyFAI .poni parameters -> ImageD11 .par parameters.

    Parameters
    ----------
    poni : dict
        Keys: dist, poni1, poni2, rot1, rot2, rot3,
        pixel1, pixel2, wavelength, orientation.
        May optionally contain ``_mirror_used`` (bool) indicating which
        solution family was used in the forward conversion.
        All lengths and wavelength in meters.
    detector_shape : (slow_dim, fast_dim) tuple, optional
        Detector pixel dimensions matching pyFAI's C-order shape
        convention: shape[0] = slow (height/rows), shape[1] = fast
        (width/columns).  Defaults to square inferred from poni.

    Returns
    -------
    dict
        Keys: distance, y_center, z_center, y_size,
        z_size, tilt_x, tilt_y, tilt_z, o11, o12, o21, o22,
        wavelength, wedge, chi, omegasign, fit_tolerance.
    """
    L = float(poni["dist"])
    rot1 = float(poni.get("rot1", 0.0))
    rot2 = float(poni.get("rot2", 0.0))
    rot3 = float(poni.get("rot3", 0.0))
    poni1 = float(poni["poni1"])
    poni2 = float(poni["poni2"])
    pv = float(poni["pixel1"])
    ph = float(poni["pixel2"])
    orientation = int(poni.get("orientation", 3))
    wl_m = float(poni.get("wavelength", 0.0))

    forward_o11 = poni.get("_forward_o11")
    forward_o22 = poni.get("_forward_o22")
    if forward_o11 is not None and forward_o22 is not None:
        o11, o22 = int(forward_o11), int(forward_o22)
        o12 = 0
        o21 = 0
    else:
        o11, o12, o21, o22 = orientation_to_flip(orientation)

    mirror_orient = poni.get("_mirror_orient")
    if mirror_orient is not None:
        mirror_M = _get_mirror_matrix(mirror_orient)
    elif poni.get("_mirror_used", True):
        mirror_M = _get_mirror_matrix(orientation)
    else:
        mirror_M = np.eye(3)

    tr1, tr2, tr3 = _compute_id11_from_pyfai(
        rot1, rot2, rot3, orientation, mirror_M=mirror_M,
        o11=o11, o22=o22)

    tx = tr3
    ty = tr2
    tz = -tr1

    delta = L / (cos(rot1) * cos(rot2))

    if detector_shape is None:
        shape_fast = shape_slow = max(int(2 * max(abs(poni1/pv), abs(poni2/ph)) + 2), 2)
        detector_shape = (shape_slow, shape_fast)
    else:
        shape_slow, shape_fast = int(detector_shape[0]), int(detector_shape[1])

    max_d1 = shape_slow - 1.0
    max_d2 = shape_fast - 1.0

    if orientation in (2, 1):
        zc = max_d1 + 0.5 - (poni1 + L * tan(rot2) / cos(rot1)) / pv
    else:
        zc = (poni1 + L * tan(rot2) / cos(rot1)) / pv - 0.5

    if orientation in (4, 1):
        yc = max_d2 + 0.5 - (poni2 - L * tan(rot1)) / ph
    else:
        yc = (poni2 - L * tan(rot1)) / ph - 0.5

    return {
        "distance": delta,
        "y_center": yc,
        "z_center": zc,
        "y_size": ph,
        "z_size": pv,
        "tilt_x": tx,
        "tilt_y": ty,
        "tilt_z": tz,
        "o11": o11,
        "o12": o12,
        "o21": o21,
        "o22": o22,
        "wavelength": wl_m,
        "wedge": 0.0,
        "chi": 0.0,
        "omegasign": 1.0,
        "fit_tolerance": 0.05,
    }


# ---------------------------------------------------------------------------
# File I/O -- .par
# ---------------------------------------------------------------------------

_PAR_GEOMETRY_KEYS = [
    "distance", "y_center", "z_center", "y_size", "z_size",
    "tilt_x", "tilt_y", "tilt_z",
    "o11", "o12", "o21", "o22",
    "wavelength", "wedge", "chi", "omegasign",
    "fit_tolerance", "min_bin_prob", "no_bins", "weight_hist_intensities",
    "t_x", "t_y", "t_z",
]

_PAR_KEY_ORDER = {name: idx for idx, name in enumerate(_PAR_GEOMETRY_KEYS)}


def read_par(filepath, par_length_unit="um"):
    """Read an ImageD11 .par file. Returns dict with all lengths in meters."""
    unit_name, unit_factor = _parse_length_unit(par_length_unit)
    par = {}
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split(None, 1)
            if len(parts) != 2:
                continue
            key = parts[0].replace("-", "_").strip()
            value = parts[1].strip()
            par[key] = value

    result = {}
    for key, value in par.items():
        try:
            result[key] = int(value)
        except (ValueError, TypeError):
            try:
                result[key] = float(value)
            except (ValueError, TypeError):
                result[key] = value

    length_keys = {"distance", "y_size", "z_size", "t_x", "t_y", "t_z"}
    for k in length_keys & result.keys():
        if isinstance(result[k], (int, float)):
            result[k] = float(result[k]) / unit_factor

    if "wavelength" in result and isinstance(result["wavelength"], (int, float)):
        result["wavelength"] = float(result["wavelength"]) / _WAVELENGTH_A_PER_M

    return result


def write_par(par, filepath, par_length_unit="um"):
    """Write an ImageD11 .par file. Lengths in meters internally."""
    unit_name, unit_factor = _parse_length_unit(par_length_unit)
    length_keys = {"distance", "y_size", "z_size", "t_x", "t_y", "t_z"}

    out = {}
    for key, value in par.items():
        if isinstance(value, float) and key in length_keys:
            out[key] = value * unit_factor
        elif key == "wavelength" and isinstance(value, (int, float)):
            out[key] = value * _WAVELENGTH_A_PER_M
        else:
            out[key] = value

    def _sort_key(k):
        return (_PAR_KEY_ORDER.get(k, 9999), k)

    lines = []
    for k in sorted(out.keys(), key=_sort_key):
        v = out[k]
        lines.append(f"{k} {v!r}" if isinstance(v, float) else f"{k} {v}")

    with open(filepath, "w") as fh:
        fh.write("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# File I/O -- .poni
# ---------------------------------------------------------------------------

def _detector_config_from_poni(poni):
    pv = float(poni["pixel1"])
    ph = float(poni["pixel2"])
    orientation = int(poni.get("orientation", 3))
    config = {"pixel1": pv, "pixel2": ph, "max_shape": None, "orientation": orientation}
    for meta_key in ("_mirror_used", "_forward_o11", "_forward_o22", "_mirror_orient"):
        if meta_key in poni:
            config[meta_key] = poni[meta_key]
    return config


def read_poni(filepath):
    """Read a pyFAI .poni file (v1 or v2/v3). Returns dict with lengths in m."""
    data = {}
    with open(filepath) as fh:
        for line in fh:
            line = line.strip()
            if not line or ":" not in line:
                continue
            if line.startswith("#"):
                continue
            words = line.split(":", 1)
            key = words[0].strip().lower()
            value = words[1].strip()
            data[key] = value

    version = float(data.get("poni_version", 1))
    orientation = 3
    pixel1 = pixel2 = None
    dc = {}
    if "detector_config" in data and version >= 2:
        try:
            dc = json.loads(data["detector_config"])
        except (json.JSONDecodeError, TypeError):
            dc = {}
        pixel1 = dc.get("pixel1")
        pixel2 = dc.get("pixel2")
        orientation = dc.get("orientation", 3)
    else:
        pixel1 = float(data.get("pixelsize1", 0))
        pixel2 = float(data.get("pixelsize2", 0))

    result = {
        "dist": float(data.get("distance", 0)),
        "poni1": float(data.get("poni1", 0)),
        "poni2": float(data.get("poni2", 0)),
        "rot1": float(data.get("rot1", 0)),
        "rot2": float(data.get("rot2", 0)),
        "rot3": float(data.get("rot3", 0)),
        "pixel1": pixel1,
        "pixel2": pixel2,
        "wavelength": float(data.get("wavelength", 0)),
        "orientation": int(orientation),
    }
    for meta_key in ("_mirror_used", "_forward_o11", "_forward_o22", "_mirror_orient"):
        if meta_key in dc:
            result[meta_key] = dc[meta_key]
    return result


def write_poni(poni, filepath):
    """Write a pyFAI v2.1 .poni file."""
    import time
    detector_config = _detector_config_from_poni(poni)
    lines = [
        "# Nota: C-Order, 1 refers to the Y axis, 2 to the X axis",
        f"# Calibration done with par_to_poni.py on {time.ctime()}",
        "poni_version: 2.1",
        "Detector: Detector",
        f"Detector_config: {json.dumps(detector_config)}",
        f"Distance: {float(poni['dist']):.12e}",
        f"Poni1: {float(poni['poni1']):.12e}",
        f"Poni2: {float(poni['poni2']):.12e}",
        f"Rot1: {float(poni['rot1']):.12e}",
        f"Rot2: {float(poni['rot2']):.12e}",
        f"Rot3: {float(poni['rot3']):.12e}",
    ]
    wl = poni.get("wavelength")
    if wl is not None:
        lines.append(f"Wavelength: {float(wl):.12e}")
    lines.append("")
    with open(filepath, "w") as fh:
        fh.write("\n".join(lines))
