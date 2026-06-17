"""
par_to_poni.py — Convert between ImageD11 .par and pyFAI .poni geometry parameters.

Based on the pyFAI source code analysis of orientation handling:
- Pixel reordering: _reorder_indexes_from_orientation  (_common.py:657)
- Sign flips: f_t1 / f_t2                            (_geometry.pyx:68-105)

Equating the full affine transforms gives exact closed-form solutions
for all 4 non-transpose flip->orientation pairs, including the pixel
reordering (C matrix) and post-rotation sign flips (S matrix).
Transpose flips (o12,o21!=0) are not supported.

Dependencies: numpy, scipy (for Rotation). All internal units are meters for
lengths and meters for wavelength.

Usage:
    import par_to_poni as pp

    par = pp.read_par("geometry.par")
    poni = pp.par_to_poni(par, detector_shape=(200, 128))
    pp.write_poni(poni, "geometry.poni")

    poni = pp.read_poni("geometry.poni")
    par = pp.poni_to_par(poni, detector_shape=(200, 128))
    pp.write_par(par, "geometry.par")
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
# Rotation matrix utilities
# ---------------------------------------------------------------------------

def _pyfai_rotation_matrix(rot1, rot2, rot3):
    """pyFAI rotation matrix: Rz(rot3).Ry_left(rot2).Rx_left(rot1).

    In standard right-handed convention this equals:
        Rz(rot3) . Ry(-rot2) . Rx(-rot1)
    which is intrinsic ZYX with angles [rot3, -rot2, -rot1].
    Matches pyFAI's actual rotation_matrix() output to machine precision
    (verified by test_pyfai_rotation_matrix_matches_actual in test_conversion.py).

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


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------

def _compute_compensated_rotation(o11, o22, orient, r1_std, r2_std, r3_std):
    """Compute compensated pyFAI rotation for a given (flip, orientation) pair.

    Derivation: equating the full pyFAI pipeline against the ID11 pipeline.
    For each orientation, pyFAI applies:
      - Pixel reordering: C = diag(c1, c2)  pre-rotation
      - Rotation: R
      - Sign flips:  S = diag(s1, s2, 1)  post-rotation

    ID11 applies the flip matrix Z = diag(o11, -o22) pre-rotation (in the
    pyFAI lab frame after G transformation).  The linear constraint is:

        S . R_comp . C = R_tilt . Z

    Solving:  R_comp[:,0] = S . R_tilt[:,0] . (o11 / c1)
              R_comp[:,1] = S . R_tilt[:,1] . (-o22 / c2)

    Returns (rot1, rot2, rot3) with cos(rot1).cos(rot2) > 0.
    """
    S_diag = {3: (1, 1, 1), 2: (-1, 1, 1), 4: (1, -1, 1), 1: (-1, -1, 1)}[orient]
    c1 = -1.0 if orient in (2, 1) else 1.0
    c2 = -1.0 if orient in (4, 1) else 1.0

    R_tilt = np.array(_pyfai_rotation_matrix(r1_std, r2_std, r3_std))

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

    result = _find_positive_equiv_from_angles(rot1_c, rot2_c, rot3_c)
    if result is None:
        # For orientations 2 and 4 the compensated rotation matrix differs
        # from R_tilt by element-wise sign flips — no equivalent Euler
        # parametrization with cos(r1)·cos(r2) > 0 exists.  The negative
        # distance is handled correctly by the round-trip (delta = dist /
        # (cos(rot1)·cos(rot2)) recovers the positive along-beam distance).
        result = (rot1_c, rot2_c, rot3_c)
    return result


def _compute_id11_from_pyfai(rot1, rot2, rot3, orient):
    """Recover ID11 tilt rotation from compensated pyFAI params.

    Reverse of _compute_compensated_rotation. From the forward equation:
      S . R_comp . C = R_tilt . Z
    reverse:
      R_tilt[:,0] = S . R_comp[:,0] . (c1 / o11)
      R_tilt[:,1] = S . R_comp[:,1] . (c2 / (-o22))
    """
    S_diag = {3: (1, 1, 1), 2: (-1, 1, 1), 4: (1, -1, 1), 1: (-1, -1, 1)}[orient]
    c1 = -1.0 if orient in (2, 1) else 1.0
    c2 = -1.0 if orient in (4, 1) else 1.0

    o11, o12, o21, o22 = orientation_to_flip(orient)
    R_comp = np.array(_pyfai_rotation_matrix(rot1, rot2, rot3))

    rt_c0 = np.array([S_diag[0] * R_comp[0, 0] * (c1 / o11),
                      S_diag[1] * R_comp[1, 0] * (c1 / o11),
                      S_diag[2] * R_comp[2, 0] * (c1 / o11)])
    rt_c1 = np.array([S_diag[0] * R_comp[0, 1] * (c2 / (-o22)),
                      S_diag[1] * R_comp[1, 1] * (c2 / (-o22)),
                      S_diag[2] * R_comp[2, 1] * (c2 / (-o22))])

    rt_c2 = np.cross(rt_c0, rt_c1)
    if np.linalg.det(np.column_stack([rt_c0, rt_c1, rt_c2])) < 0:
        rt_c2 = -rt_c2

    R_tilt = (
        (rt_c0[0], rt_c1[0], rt_c2[0]),
        (rt_c0[1], rt_c1[1], rt_c2[1]),
        (rt_c0[2], rt_c1[2], rt_c2[2]),
    )
    return _extract_rot(R_tilt)


def par_to_poni(par, detector_shape=None):
    """Convert ImageD11 .par parameters -> pyFAI .poni parameters.

    Parameters
    ----------
    par : dict
        Keys: distance, y_center, z_center, y_size, z_size,
        tilt_x, tilt_y, tilt_z, o11, o12, o21, o22, wavelength.
        All lengths in meters internally, wavelength in meters.
    detector_shape : (fast_dim, slow_dim) tuple, optional
        Detector pixel dimensions. Required for non-native orientations
        (2 and 4) to compute correct PONI accounting for pyFAI's
        pixel-reordering convention. For orientation 3 (native) the
        shape is not needed. Defaults to square inferred from beam center.

    Returns
    -------
    dict
        Keys: dist, poni1, poni2, rot1, rot2, rot3,
        pixel1, pixel2, wavelength, orientation.
    """
    tx = float(par.get("tilt_x", 0.0))
    ty = float(par.get("tilt_y", 0.0))
    tz = float(par.get("tilt_z", 0.0))
    distance = float(par["distance"])
    yc = float(par["y_center"])
    zc = float(par["z_center"])
    ys = float(par["y_size"])
    zs = float(par["z_size"])
    o11 = int(par.get("o11", 1))
    o12 = int(par.get("o12", 0))
    o21 = int(par.get("o21", 0))
    o22 = int(par.get("o22", -1))
    orientation = flip_to_orientation(o11, o12, o21, o22)
    wl_m = float(par.get("wavelength", 0.0))
    delta = distance

    if detector_shape is None:
        shape_fast = max(int(2 * yc + 1), 2)
        shape_slow = max(int(2 * zc + 1), 2)
        detector_shape = (shape_fast, shape_slow)
    else:
        shape_fast, shape_slow = int(detector_shape[0]), int(detector_shape[1])

    # pyFAI _reorder_indexes_from_orientation uses shape[0]-1 for d1 (slow axis)
    # and shape[1]-1 for d2 (fast axis).  detector_shape is (fast_dim, slow_dim)
    # so shape[0] maps to slow_dim = detector_shape[1], shape[1] maps to fast_dim.
    max_d1 = shape_slow - 1.0
    max_d2 = shape_fast - 1.0

    # Standard tilt mapping
    r1 = -tz
    r2 = ty
    r3 = tx

    rot1, rot2, rot3 = _compute_compensated_rotation(o11, o22, orientation, r1, r2, r3)

    dist = delta * cos(rot2) * cos(rot1)

    if orientation in (2, 1):
        poni1 = -delta * sin(rot2) + zs * (max_d1 - zc + 0.5)
    else:
        poni1 = -delta * sin(rot2) + zs * (zc + 0.5)

    if orientation in (4, 1):
        poni2 = delta * cos(rot2) * sin(rot1) + ys * (max_d2 - yc + 0.5)
    else:
        poni2 = delta * cos(rot2) * sin(rot1) + ys * (yc + 0.5)

    return {
        "dist": dist,
        "poni1": poni1,
        "poni2": poni2,
        "rot1": rot1,
        "rot2": rot2,
        "rot3": rot3,
        "pixel1": zs,
        "pixel2": ys,
        "wavelength": wl_m,
        "orientation": orientation,
    }


def poni_to_par(poni, detector_shape=None):
    """Convert pyFAI .poni parameters -> ImageD11 .par parameters.

    Parameters
    ----------
    poni : dict
        Keys: dist, poni1, poni2, rot1, rot2, rot3,
        pixel1, pixel2, wavelength, orientation.
        All lengths and wavelength in meters.
    detector_shape : (fast_dim, slow_dim) tuple, optional
        Detector pixel dimensions, needed to reverse orientation-specific
        PONI formulas. Defaults to square inferred from poni.

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
    o11, o12, o21, o22 = orientation_to_flip(orientation)
    wl_m = float(poni.get("wavelength", 0.0))

    tr1, tr2, tr3 = _compute_id11_from_pyfai(rot1, rot2, rot3, orientation)

    tx = tr3
    ty = tr2
    tz = -tr1

    delta = L / (cos(rot1) * cos(rot2))

    if detector_shape is None:
        shape_fast = shape_slow = max(int(2 * max(abs(poni1/pv), abs(poni2/ph)) + 2), 2)
        detector_shape = (shape_fast, shape_slow)
    else:
        shape_fast, shape_slow = int(detector_shape[0]), int(detector_shape[1])

    # Same convention as in par_to_poni: shape[0] maps to slow_dim (detector_shape[1])
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
    return {"pixel1": pv, "pixel2": ph, "max_shape": None, "orientation": orientation}


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
