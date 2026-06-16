"""
par_to_poni.py — Convert between ImageD11 .par and pyFAI .poni geometry parameters.

Based on the 4x4 affine matrix analysis of the two geometry frameworks.
All 4 non-transpose flip→orientation pairs have exact closed-form solutions
via rotation compensation. Transpose flips (o12,o21≠0) are not supported.

No dependencies beyond Python stdlib. All internal units are meters for lengths
and meters for wavelength.

Usage:
    import par_to_poni as pp

    par = pp.read_par("geometry.par")
    poni = pp.par_to_poni(par)
    pp.write_poni(poni, "geometry.poni")

    poni = pp.read_poni("geometry.poni")
    par = pp.poni_to_par(poni)
    pp.write_par(par, "geometry.par")
"""

import json
import math
from math import cos, sin, tan, atan2, pi, sqrt


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
    """pyFAI rotation matrix:  Rz(rot3)·Ry_left(rot2)·Rx_left(rot1).

    This equals Rz(rot3)·Ry(-rot2)·Rx(-rot1) in standard right-handed
    convention, matching pyFAI's rotation_matrix() output.
    """
    s1, c1 = sin(rot1), cos(rot1)
    s2, c2 = sin(rot2), cos(rot2)
    s3, c3 = sin(rot3), cos(rot3)

    Rx = ((1, 0, 0), (0, c1, s1), (0, -s1, c1))
    Ry = ((c2, 0, -s2), (0, 1, 0), (s2, 0, c2))
    Rz = ((c3, -s3, 0), (s3, c3, 0), (0, 0, 1))

    a00 = Ry[0][0] * Rx[0][0] + Ry[0][1] * Rx[1][0] + Ry[0][2] * Rx[2][0]
    a01 = Ry[0][0] * Rx[0][1] + Ry[0][1] * Rx[1][1] + Ry[0][2] * Rx[2][1]
    a02 = Ry[0][0] * Rx[0][2] + Ry[0][1] * Rx[1][2] + Ry[0][2] * Rx[2][2]
    a10 = Ry[1][0] * Rx[0][0] + Ry[1][1] * Rx[1][0] + Ry[1][2] * Rx[2][0]
    a11 = Ry[1][0] * Rx[0][1] + Ry[1][1] * Rx[1][1] + Ry[1][2] * Rx[2][1]
    a12 = Ry[1][0] * Rx[0][2] + Ry[1][1] * Rx[1][2] + Ry[1][2] * Rx[2][2]
    a20 = Ry[2][0] * Rx[0][0] + Ry[2][1] * Rx[1][0] + Ry[2][2] * Rx[2][0]
    a21 = Ry[2][0] * Rx[0][1] + Ry[2][1] * Rx[1][1] + Ry[2][2] * Rx[2][1]
    a22 = Ry[2][0] * Rx[0][2] + Ry[2][1] * Rx[1][2] + Ry[2][2] * Rx[2][2]

    return (
        (Rz[0][0] * a00 + Rz[0][1] * a10 + Rz[0][2] * a20,
         Rz[0][0] * a01 + Rz[0][1] * a11 + Rz[0][2] * a21,
         Rz[0][0] * a02 + Rz[0][1] * a12 + Rz[0][2] * a22),
        (Rz[1][0] * a00 + Rz[1][1] * a10 + Rz[1][2] * a20,
         Rz[1][0] * a01 + Rz[1][1] * a11 + Rz[1][2] * a21,
         Rz[1][0] * a02 + Rz[1][1] * a12 + Rz[1][2] * a22),
        (Rz[2][0] * a00 + Rz[2][1] * a10 + Rz[2][2] * a20,
         Rz[2][0] * a01 + Rz[2][1] * a11 + Rz[2][2] * a21,
         Rz[2][0] * a02 + Rz[2][1] * a12 + Rz[2][2] * a22),
    )


def _extract_rot(R):
    """Extract rot1,rot2,rot3 from a pyFAI rotation matrix.

    R = Rz(rot3)·Ry_left(rot2)·Rx_left(rot1).
    Entries: R[2,0]=sin(rot2), R[2,1]=-cos(rot2)·sin(rot1),
             R[2,2]=cos(rot2)·cos(rot1),
             R[1,0]=sin(rot3)·cos(rot2), R[0,0]=cos(rot3)·cos(rot2).
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


def _find_positive_equiv(R_comp):
    """Find an equivalent (rot1,rot2,rot3) for R_comp with cos(rot1)·cos(rot2)>0.

    Returns (rot1, rot2, rot3) or None if no positive-dist equivalent found.
    """
    a, b, c = _extract_rot(R_comp)
    best = None
    for d1 in (0, pi, -pi, 2 * pi):
        for d2 in (0, pi, -pi, 2 * pi):
            for s2 in (1, -1):
                rt1, rt2, rt3 = a + d1, s2 * b + d2, c
                if abs(rt1) > 10 or abs(rt2) > 10:
                    continue
                Rt = _pyfai_rotation_matrix(rt1, rt2, rt3)
                maxdiff = max(abs(Rt[i][j] - R_comp[i][j])
                             for i in range(3) for j in range(3))
                if maxdiff < 1e-8:
                    dc = cos(rt1) * cos(rt2)
                    if dc > 0:
                        if best is None or abs(rt1) + abs(rt2) < abs(best[0]) + abs(best[1]):
                            best = (rt1, rt2, rt3)
    return best


# ---------------------------------------------------------------------------
# Core conversion
# ---------------------------------------------------------------------------

def _compute_compensated_rotation(o11, o22, orient, r1_std, r2_std, r3_std):
    """Compute compensated pyFAI rotation for a given (flip, orientation) pair.

    Uses the linear system: S·R_comp·C_actual = R_tilt·Z where
    C_actual = diag(1,1) (the effective linear map after pixel reordering).

    Returns (rot1, rot2, rot3) with cos(rot1)·cos(rot2) > 0.
    """
    S_diag = {3: (1, 1, 1), 2: (-1, 1, 1), 4: (1, -1, 1), 1: (-1, -1, 1)}[orient]
    R_tilt = _pyfai_rotation_matrix(r1_std, r2_std, r3_std)

    # Z = [[o11, 0], [0, -o22], [0, 0]] (no pixel-size scaling needed)
    # target = S · R_tilt · Z  (3×2)
    # R_comp[:, 0:2] = target (since C_actual = diag(1,1))

    z_col0 = (o11, 0.0, 0.0)
    z_col1 = (0.0, float(-o22), 0.0)

    # S · R_tilt · z_col0
    Rt_z0 = (
        R_tilt[0][0] * z_col0[0] + R_tilt[0][1] * z_col0[1] + R_tilt[0][2] * z_col0[2],
        R_tilt[1][0] * z_col0[0] + R_tilt[1][1] * z_col0[1] + R_tilt[1][2] * z_col0[2],
        R_tilt[2][0] * z_col0[0] + R_tilt[2][1] * z_col0[1] + R_tilt[2][2] * z_col0[2],
    )
    target_c0 = (
        S_diag[0] * Rt_z0[0], S_diag[1] * Rt_z0[1], S_diag[2] * Rt_z0[2],
    )

    Rt_z1 = (
        R_tilt[0][0] * z_col1[0] + R_tilt[0][1] * z_col1[1] + R_tilt[0][2] * z_col1[2],
        R_tilt[1][0] * z_col1[0] + R_tilt[1][1] * z_col1[1] + R_tilt[1][2] * z_col1[2],
        R_tilt[2][0] * z_col1[0] + R_tilt[2][1] * z_col1[1] + R_tilt[2][2] * z_col1[2],
    )
    target_c1 = (
        S_diag[0] * Rt_z1[0], S_diag[1] * Rt_z1[1], S_diag[2] * Rt_z1[2],
    )

    # Cross product for third column
    r2c0 = target_c0[1] * target_c1[2] - target_c0[2] * target_c1[1]
    r2c1 = target_c0[2] * target_c1[0] - target_c0[0] * target_c1[2]
    r2c2 = target_c0[0] * target_c1[1] - target_c0[1] * target_c1[0]

    det = (target_c0[0] * (target_c1[1] * r2c2 - target_c1[2] * r2c1)
           - target_c1[0] * (target_c0[1] * r2c2 - target_c0[2] * r2c1)
           + r2c0 * (target_c0[1] * target_c1[2] - target_c0[2] * target_c1[1]))
    if det < 0:
        r2c0, r2c1, r2c2 = -r2c0, -r2c1, -r2c2

    R_comp = (
        (target_c0[0], target_c1[0], r2c0),
        (target_c0[1], target_c1[1], r2c1),
        (target_c0[2], target_c1[2], r2c2),
    )

    result = _find_positive_equiv(R_comp)
    if result is None:
        result = _extract_rot(R_comp)
    return result


def _compute_id11_from_pyfai(rot1, rot2, rot3, orient):
    """Recover original ID11 tilt rotation from compensated pyFAI params.

    From the forward equation: S(orient) · R_comp = R_tilt · Z(flip),
    we reverse to get R_tilt[:,0:2] = S · R_comp · Z and build the
    third column via cross product, ensuring right-handed orientation.
    Returns (tr1, tr2, tr3) as the tilt rotation angles (R3·R2·R1 convention).
    """
    S_diag = {3: (1, 1, 1), 2: (-1, 1, 1), 4: (1, -1, 1), 1: (-1, -1, 1)}[orient]
    R_total = _pyfai_rotation_matrix(rot1, rot2, rot3)

    # Build R_tilt from R_total.  Since R_comp[:,0:2] = S · R_tilt · Z
    # and Z^T·Z = I₂ (because o11²=o22²=1), we have
    # R_tilt[:,0:2] = S · R_comp · Z.  For the full R_tilt,

    o11, o12, o21, o22 = orientation_to_flip(orient)
    z_col0 = (o11, 0.0, 0.0)
    z_col1 = (0.0, float(-o22), 0.0)

    # R_total · Z gives R_total[:, 0]*o11 and R_total[:, 1]*(-o22)
    Rt0 = (
        R_total[0][0] * z_col0[0], R_total[1][0] * z_col0[0], R_total[2][0] * z_col0[0],
    )
    Rt1 = (
        R_total[0][1] * z_col1[1], R_total[1][1] * z_col1[1], R_total[2][1] * z_col1[1],
    )
    # S^{-1} = S:  Rt_col0 = S · (R_total · z0)
    st_col0 = (S_diag[0] * Rt0[0], S_diag[1] * Rt0[1], S_diag[2] * Rt0[2])
    st_col1 = (S_diag[0] * Rt1[0], S_diag[1] * Rt1[1], S_diag[2] * Rt1[2])

    # Cross product for third column of R_tilt
    # R_tilt[:, 2] should be the cross product of the first two columns
    r20 = st_col0[1] * st_col1[2] - st_col0[2] * st_col1[1]
    r21 = st_col0[2] * st_col1[0] - st_col0[0] * st_col1[2]
    r22 = st_col0[0] * st_col1[1] - st_col0[1] * st_col1[0]

    det = (st_col0[0] * (st_col1[1] * r22 - st_col1[2] * r21)
           - st_col1[0] * (st_col0[1] * r22 - st_col0[2] * r21)
           + r20 * (st_col0[1] * st_col1[2] - st_col0[2] * st_col1[1]))
    if det < 0:
        r20, r21, r22 = -r20, -r21, -r22

    R_tilt = (
        (st_col0[0], st_col1[0], r20),
        (st_col0[1], st_col1[1], r21),
        (st_col0[2], st_col1[2], r22),
    )

    return _extract_rot(R_tilt)


def par_to_poni(par):
    """Convert ImageD11 .par parameters -> pyFAI .poni parameters.

    Parameters
    ----------
    par : dict
        Keys: distance, y_center, z_center, y_size, z_size,
        tilt_x, tilt_y, tilt_z, o11, o12, o21, o22, wavelength.
        All lengths in meters internally, wavelength in meters.

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

    # Standard tilt mapping
    r1 = -tz
    r2 = ty
    r3 = tx

    # Compute compensated rotation (with positive distance)
    rot1, rot2, rot3 = _compute_compensated_rotation(o11, o22, orientation, r1, r2, r3)

    delta = distance  # ImageD11 along-beam distance
    dist = delta * cos(rot2) * cos(rot1)  # orthogonal distance to PONI

    # Standard PONI formulas (same for all orientations).
    # The orientation-specific pixel reordering is handled by the
    # orientation flag in pyFAI; the PONI parameters are in the
    # native un-flipped coordinate system.
    poni1 = -delta * sin(rot2) + zs * (zc + 0.5)
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


def poni_to_par(poni):
    """Convert pyFAI .poni parameters -> ImageD11 .par parameters.

    Parameters
    ----------
    poni : dict
        Keys: dist, poni1, poni2, rot1, rot2, rot3,
        pixel1, pixel2, wavelength, orientation.
        All lengths and wavelength in meters.

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

    # Recover ID11 tilts from the compensated pyFAI rotation.
    # tr1,tr2,tr3 = TILT rotation (uncompensated), i.e. the original tilts.
    tr1, tr2, tr3 = _compute_id11_from_pyfai(rot1, rot2, rot3, orientation)

    tx = tr3
    ty = tr2
    tz = -tr1

    # Recover along-beam distance and beam center using the COMPENSATED
    # rotation parameters (rot1,rot2,rot3 are COMPENSATED, matching
    # what was used to compute PONI in the forward direction).
    delta = L / (cos(rot1) * cos(rot2))

    zc = (poni1 + L * tan(rot2) / cos(rot1)) / pv - 0.5
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
# File I/O — .poni
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
        f"Distance: {poni['dist']!r}",
        f"Poni1: {poni['poni1']!r}",
        f"Poni2: {poni['poni2']!r}",
        f"Rot1: {poni['rot1']!r}",
        f"Rot2: {poni['rot2']!r}",
        f"Rot3: {poni['rot3']!r}",
    ]
    wl = poni.get("wavelength")
    if wl is not None:
        lines.append(f"Wavelength: {wl!r}")
    lines.append("")
    with open(filepath, "w") as fh:
        fh.write("\n".join(lines))
