"""
par_to_poni.py — Convert between ImageD11 .par and pyFAI .poni geometry parameters.

Coordinate systems:
- pyFAI: axis1=y_up(slow), axis2=x_starboard(fast), axis3=z_downstream
  Rotation: Rz(rot3)·Ry(−rot2)·Rx(−rot1)
- ImageD11: X=downstream, Y=port, Z=up(slow)
  Rotation: Rx(tilt_x)·Ry(tilt_y)·Rz(tilt_z)

Frame transform between systems:  G = [[0,0,1],[0,-1,0],[1,0,0]]

Conversion formula:
    rot1 = −tilt_z    rot2 = tilt_y    rot3 = tilt_x
    orient = flip_to_orientation(o11, o12, o21, o22)
    distance, PONI from orientation-dependent formulas below.

Each of the 4 non-transpose flip matrices maps to exactly one pyFAI
orientation (1–4).  Transpose flips (o12, o21 ≠ 0) are not supported.
Spatial distortion is not handled.

Dependencies: numpy.  All internal units are metres for lengths,
metres for wavelength.

Usage:
    import par_to_poni as pp

    par = pp.read_par("geometry.par")
    poni = pp.par_to_poni(par, detector_shape=(2162, 2068))
    pp.write_poni(poni, "geometry.poni")

    poni = pp.read_poni("geometry.poni")
    par = pp.poni_to_par(poni, detector_shape=(2162, 2068))
    pp.write_par(par, "geometry.par")

    # Convert azimuth angles between the two programs:
    eta_rad = pp.chi_to_eta(chi_rad, orientation=3)
    chi_rad = pp.eta_to_chi(eta_rad, orientation=3)
    # Orientation can come from a par or poni dict:
    eta_rad = pp.chi_to_eta(chi_rad, par)
    chi_rad = pp.eta_to_chi(chi_rad, poni)
"""

import json
import math
import time
from math import cos, sin, tan, atan2

import numpy as np


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
    (-1, 0, 0, -1): 2,
    (1, 0, 0, 1): 4,
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

_CHI_ETA_SIN_COS_FACTORS = {
    3: (+1, +1),   # chi =  90° − η    sin(chi) = +cos(eta)  cos(chi) = +sin(eta)
    2: (-1, +1),   # chi =  eta − 90°    sin(chi) = −cos(eta)  cos(chi) = +sin(eta)
    4: (+1, -1),   # chi =  eta + 90°    sin(chi) = +cos(eta)  cos(chi) = −sin(eta)
    1: (-1, -1),   # chi = 270° − eta    sin(chi) = −cos(eta)  cos(chi) = −sin(eta)
}


def _azimuth_factors(orientation):
    """Return (sin_factor, cos_factor) for chi ↔ eta given an orientation.

    sin(chi) = sin_factor · cos(eta)
    cos(chi) = cos_factor · sin(eta)
    """
    _m = {
        3: (+1, +1),   # identity — no flip
        2: (-1, +1),   # flip slow axis
        4: (+1, -1),   # flip fast axis
        1: (-1, -1),   # flip both
    }
    return _m.get(orientation, (+1, +1))


def chi_to_eta(chi_rad, orientation):
    """Convert pyFAI chi → ImageD11 eta (radians).

    orientation : int (1–4), par dict, or poni dict.

    ┌─────────┬──────────────────┬─────────────────────────┐
    │ orient  │ chi = f(eta)      │ sin(chi), cos(chi)      │
    ├─────────┼──────────────────┼─────────────────────────┤
    │ 3       │ chi =  90° − eta  │ (+cos eta, +sin eta)    │
    │ 2       │ chi =  eta − 90°  │ (−cos eta, +sin eta)    │
    │ 4       │ chi =  eta + 90°  │ (+cos eta, −sin eta)    │
    │ 1       │ chi = 270° − eta  │ (−cos eta, −sin eta)    │
    └─────────┴──────────────────┴─────────────────────────┘
    """
    orientation = _extract_orientation_from_arg(orientation)
    s0, s1 = _CHI_ETA_SIN_COS_FACTORS[orientation]
    return np.arctan2(s1 * np.cos(chi_rad), s0 * np.sin(chi_rad))


def eta_to_chi(eta_rad, orientation):
    """Convert ImageD11 eta → pyFAI chi (radians).

    Inverse of chi_to_eta modulo 2π.
    """
    orientation = _extract_orientation_from_arg(orientation)
    s0, s1 = _CHI_ETA_SIN_COS_FACTORS[orientation]
    return np.arctan2(s0 * np.cos(eta_rad), s1 * np.sin(eta_rad))


def _extract_orientation_from_arg(arg):
    """Extract pyFAI orientation (1–4) from a par dict, poni dict, or int."""
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
# Core conversion
# ---------------------------------------------------------------------------

def par_to_poni(par, detector_shape=None):
    """Convert ImageD11 .par parameters → pyFAI .poni parameters.

    Parameters
    ----------
    par : dict
        Keys: distance, y_center, z_center, y_size, z_size,
        tilt_x, tilt_y, tilt_z, o11, o12, o21, o22, wavelength.
        All lengths in metres, wavelength in metres.
    detector_shape : (slow_dim, fast_dim) tuple, optional
        pyFAI C-order shape (rows, cols).  Defaults to square
        inferred from beam centre.

    Returns
    -------
    dict
        Keys: dist, poni1, poni2, rot1, rot2, rot3,
        pixel1, pixel2, wavelength, orientation.
    """
    o11 = int(par.get("o11", 1))
    o12 = int(par.get("o12", 0))
    o21 = int(par.get("o21", 0))
    o22 = int(par.get("o22", -1))
    orientation = flip_to_orientation(o11, o12, o21, o22)

    # Direct rot → tilt mapping:  rot1 = −tilt_z,  rot2 = tilt_y,  rot3 = tilt_x
    rot1 = -float(par.get("tilt_z", 0.0))
    rot2 = float(par.get("tilt_y", 0.0))
    rot3 = float(par.get("tilt_x", 0.0))

    delta = float(par["distance"])
    yc = float(par["y_center"])
    zc = float(par["z_center"])
    ys = float(par["y_size"])
    zs = float(par["z_size"])
    wl_m = float(par.get("wavelength", 0.0))

    if detector_shape is None:
        shape_fast = max(int(2 * (abs(yc) if yc else 500) + 1), 2)
        shape_slow = max(int(2 * (abs(zc) if zc else 500) + 1), 2)
        shape_fast = min(shape_fast, 100000)
        shape_slow = min(shape_slow, 100000)
        detector_shape = (shape_slow, shape_fast)

    shape_slow, shape_fast = int(detector_shape[0]), int(detector_shape[1])
    max_d1 = shape_slow - 1.0
    max_d2 = shape_fast - 1.0

    dist = delta * cos(rot1) * cos(rot2)

    beam_z = max_d1 - zc if orientation in (2, 1) else zc
    beam_y = max_d2 - yc if orientation in (4, 1) else yc

    poni1 = -delta * sin(rot2) + zs * (beam_z + 0.5)
    poni2 = delta * cos(rot2) * sin(rot1) + ys * (beam_y + 0.5)

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
    """Convert pyFAI .poni parameters → ImageD11 .par parameters.

    Parameters
    ----------
    poni : dict
        Keys: dist, poni1, poni2, rot1, rot2, rot3,
        pixel1, pixel2, wavelength, orientation.
        All lengths and wavelength in metres.
    detector_shape : (slow_dim, fast_dim) tuple, optional
        pyFAI C-order shape (rows, cols).

    Returns
    -------
    dict
        Keys: distance, y_center, z_center, y_size, z_size,
        tilt_x, tilt_y, tilt_z, o11, o12, o21, o22,
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

    o11, o12, o21, o22 = orientation_to_flip(orientation)

    # Direct rot → tilt mapping:  tx = rot3,  ty = rot2,  tz = −rot1
    tx = rot3
    ty = rot2
    tz = -rot1

    delta = L / (cos(rot1) * cos(rot2))

    if detector_shape is None:
        shape_fast = shape_slow = max(
            int(2 * max(abs(poni1 / pv), abs(poni2 / ph)) + 2), 2)
        detector_shape = (shape_slow, shape_fast)

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
# File I/O — .par
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
    """Read an ImageD11 .par file.  Returns dict with all lengths in metres."""
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
    """Write an ImageD11 .par file.  Lengths in metres internally."""
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
# File I/O — .poni
# ---------------------------------------------------------------------------

def _detector_config_from_poni(poni):
    pv = float(poni["pixel1"])
    ph = float(poni["pixel2"])
    orientation = int(poni.get("orientation", 3))
    config = {"pixel1": pv, "pixel2": ph, "max_shape": None,
              "orientation": orientation}
    return config


def read_poni(filepath):
    """Read a pyFAI .poni file (v1 or v2/v3).  Returns dict with lengths in m."""
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

    return {
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


def write_poni(poni, filepath):
    """Write a pyFAI v2.1 .poni file."""
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
