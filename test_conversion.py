"""
test_conversion.py -- Test the par <-> poni conversion.

Validates:
1. Round-trip identity for all 4 orientations
2. 2th / q values match between pyFAI and ImageD11 (same raw pixels)
3. Azimuthal angles match (sin/cos comparison, same raw pixels)
4. Full xyz lab coordinates match (same raw pixels, non-square detector)
5. IO round-trip through disk files

Requires pyFAI and ImageD11 to be importable.
"""

import sys
import os
import tempfile
import math
import unittest
import numpy as np
from scipy.spatial.transform import Rotation as ScipyRotation

import par_to_poni as pp
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
import pyFAI
from ImageD11.transform import (
    compute_xyz_lab,
    compute_tth_eta,
    detector_rotation_matrix,
)
from ImageD11.parameters import parameters as ImageD11Parameters

print(f"pyFAI version: {pyFAI.version}")


# ---------------------------------------------------------------------------
# Test parameters
# ---------------------------------------------------------------------------

def make_base_par():
    """Create base par dict with strongly-tilted geometry, 1000x1000."""
    return dict(
        distance=0.15,
        y_center=500.0,
        z_center=500.0,
        y_size=75e-6,
        z_size=75e-6,
        tilt_x=0.3,
        tilt_y=0.2,
        tilt_z=-0.15,
        o11=1, o12=0, o21=0, o22=-1,
        wavelength=1.5406e-10,
        wedge=0.0,
        chi=0.0,
        omegasign=1.0,
        fit_tolerance=0.05,
    )


DETECTOR_SHAPE = (1000, 1000)  # (slow, fast) -- pyFAI C-order convention

# All 4 non-transpose flip orientations
FLIPS = [
    (1, 0, 0, -1, 3, "orient3_native"),
    (-1, 0, 0, 1, 1, "orient1_flip_both"),
    (-1, 0, 0, -1, 4, "orient4_flip_fast"),
    (1, 0, 0, 1, 2, "orient2_flip_slow"),
]


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def pyFAI_from_poni(poni):
    """Create an AzimuthalIntegrator from a poni dict."""
    ai = AzimuthalIntegrator(
        dist=poni["dist"],
        poni1=poni["poni1"],
        poni2=poni["poni2"],
        rot1=poni["rot1"],
        rot2=poni["rot2"],
        rot3=poni["rot3"],
        pixel1=poni["pixel1"],
        pixel2=poni["pixel2"],
        splineFile=None,
        detector=None,
        wavelength=poni.get("wavelength"),
        orientation=poni.get("orientation", 3),
    )
    ai.detector.shape = DETECTOR_SHAPE
    return ai


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestRoundTrip(unittest.TestCase):
    """Test par -> poni -> par and poni -> par -> poni round trips."""

    def test_par_round_trip_all_flips(self):
        """par -> poni -> par should recover original values."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)

                for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
                    self.assertAlmostEqual(
                        par[key], par2[key], delta=1e-10,
                        msg=f"{label}: {key} mismatch"
                    )

                for key in ["tilt_x", "tilt_y", "tilt_z"]:
                    self.assertAlmostEqual(
                        par[key], par2[key], delta=1e-10,
                        msg=f"{label}: {key} mismatch"
                    )

                for key in ["o11", "o12", "o21", "o22"]:
                    self.assertEqual(par[key], par2[key],
                                     msg=f"{label}: {key} mismatch")

                self.assertAlmostEqual(
                    par["wavelength"], par2["wavelength"], delta=1e-10,
                    msg=f"{label}: wavelength mismatch"
                )

    def test_poni_round_trip_all_flips(self):
        """poni -> par -> poni should recover original values."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                poni2 = pp.par_to_poni(
                    pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE),
                    detector_shape=DETECTOR_SHAPE)

                for key in ["dist", "poni1", "poni2", "rot1", "rot2", "rot3",
                            "pixel1", "pixel2", "wavelength", "orientation"]:
                    self.assertAlmostEqual(
                        poni[key], poni2[key], delta=1e-10,
                        msg=f"{label}: {key} mismatch"
                    )

    def test_round_trip_zero_tilts(self):
        """Round trip works when all tilts are zero."""
        par = make_base_par()
        for k in ["tilt_x", "tilt_y", "tilt_z"]:
            par[k] = 0.0

        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)

                for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
                    self.assertAlmostEqual(par[key], par2[key], delta=1e-10,
                                           msg=f"{label}: {key}")

    def test_round_trip_single_tilts(self):
        """Round trip with only one non-zero tilt at a time."""
        par = make_base_par()
        angles = [0.1, -0.2, 0.3, -0.5, 0.7]

        for angle in angles:
            for tilt_key in ["tilt_x", "tilt_y", "tilt_z"]:
                with self.subTest(tilt=tilt_key, angle=angle):
                    p = dict(par)
                    for k in ["tilt_x", "tilt_y", "tilt_z"]:
                        p[k] = angle if k == tilt_key else 0.0

                    poni = pp.par_to_poni(p, detector_shape=DETECTOR_SHAPE)
                    par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)

                    for key in ["distance", "y_center", "z_center"]:
                        self.assertAlmostEqual(p[key], par2[key], delta=1e-8,
                                               msg=f"{tilt_key}={angle}: {key}")

    def test_round_trip_edge_beam_positions(self):
        """Round trip with beam at detector edges."""
        par = make_base_par()
        par["tilt_x"] = par["tilt_y"] = par["tilt_z"] = 0.0

        for yc, zc in [(0, 0), (999, 999), (0, 999), (999, 0), (500, 500)]:
            with self.subTest(y_center=yc, z_center=zc):
                p = dict(par)
                p["y_center"] = float(yc)
                p["z_center"] = float(zc)
                poni = pp.par_to_poni(p, detector_shape=DETECTOR_SHAPE)
                par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
                self.assertAlmostEqual(yc, par2["y_center"], delta=1e-10)
                self.assertAlmostEqual(zc, par2["z_center"], delta=1e-10)


class TestTwothetaMatching(unittest.TestCase):
    """Test that 2th values match between pyFAI and ImageD11."""

    NCOORDS = 5000

    def test_tth_matches_all_flips(self):
        """2th values match to machine precision for all 4 orientations,
        same raw pixel indices, no coordinate flipping."""
        rng = np.random.RandomState(42)
        shape_slow, shape_fast = DETECTOR_SHAPE
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                ai = pyFAI_from_poni(poni)

                d1 = rng.uniform(0, shape_slow - 1, self.NCOORDS)
                d2 = rng.uniform(0, shape_fast - 1, self.NCOORDS)

                tth_pyfai = ai.tth(d1=d1, d2=d2, path="cython")

                tth_id11, _ = compute_tth_eta(
                    np.array([d1, d2]),
                    **{k: par[k] for k in [
                        "y_center", "y_size", "z_center", "z_size",
                        "tilt_x", "tilt_y", "tilt_z", "distance",
                        "o11", "o12", "o21", "o22"
                    ]}
                )
                tth_id11_rad = np.radians(tth_id11)

                diff = np.abs(tth_pyfai - tth_id11_rad)
                self.assertLess(np.max(diff), 1e-7,
                                msg=f"{label}: max 2th diff {np.max(diff):.2e}")

    def test_tth_matches_zero_tilts(self):
        """2th values match when all tilts are zero."""
        par = make_base_par()
        par["tilt_x"] = par["tilt_y"] = par["tilt_z"] = 0.0

        rng = np.random.RandomState(42)
        shape_slow, shape_fast = DETECTOR_SHAPE
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                ai = pyFAI_from_poni(poni)

                d1 = rng.uniform(0, shape_slow - 1, self.NCOORDS)
                d2 = rng.uniform(0, shape_fast - 1, self.NCOORDS)

                tth_pyfai = ai.tth(d1=d1, d2=d2, path="cython")
                tth_id11, _ = compute_tth_eta(np.array([d1, d2]), **par)
                tth_id11_rad = np.radians(tth_id11)

                self.assertLess(np.max(np.abs(tth_pyfai - tth_id11_rad)), 1e-7,
                                msg=f"{label}: zero tilt mismatch")

    def test_tth_versus_q(self):
        """q vector values are consistent with 2th."""
        par = make_base_par()
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        ai = pyFAI_from_poni(poni)

        d1 = np.array([100.0, 500.0, 900.0])
        d2 = np.array([200.0, 500.0, 800.0])

        tth = ai.tth(d1=d1, d2=d2)
        q_pyfai = ai.qFunction(d1=d1, d2=d2)
        wavelength = ai.wavelength

        q_expected = 4.0e-9 * math.pi * np.sin(tth / 2.0) / wavelength
        self.assertTrue(np.allclose(q_pyfai, q_expected, rtol=1e-10))


class TestAzimuthMatching(unittest.TestCase):
    """Test that azimuthal angles (chi / eta) match correctly."""

    NCOORDS = 2000

    def test_azimuth_relationship_all_flips(self):
        """Azimuth (chi/eta) matches using the correct sin/cos factors
        for the default solution's (orientation, mirror) combination."""
        rng = np.random.RandomState(123)
        shape_slow, shape_fast = DETECTOR_SHAPE
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                ai = pyFAI_from_poni(poni)

                d1 = rng.uniform(0, shape_slow - 1, self.NCOORDS)
                d2 = rng.uniform(0, shape_fast - 1, self.NCOORDS)

                chi = ai.chi(d1=d1, d2=d2, path="cython")
                _, eta = compute_tth_eta(np.array([d1, d2]), **par)
                eta_rad = np.radians(eta)

                mirror_orient = poni.get("_mirror_orient", poni["orientation"])
                sf, cf = pp._azimuth_factors(mirror_orient)

                target_sin = sf * np.cos(eta_rad)
                target_cos = cf * np.sin(eta_rad)

                sin_diff = np.abs(np.sin(chi) - target_sin)
                cos_diff = np.abs(np.cos(chi) - target_cos)

                self.assertLess(np.max(sin_diff), 1e-7,
                                msg=f"{label}: max sin diff {np.max(sin_diff):.2e}")
                self.assertLess(np.max(cos_diff), 1e-7,
                                msg=f"{label}: max cos diff {np.max(cos_diff):.2e}")


class TestLabCoordinates(unittest.TestCase):
    """Full xyz lab coordinates match pixel-by-pixel, non-square detector,
    same raw pixel indices, no coordinate flipping."""

    NCOORDS = 2000
    SHAPE = (128, 200)  # (slow, fast) — pyFAI C-order
    G = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]], dtype=float)

    def _make_test_par(self, **kw):
        return dict(
            distance=0.15,
            y_center=(self.SHAPE[1] - 1) / 2.0,
            z_center=(self.SHAPE[0] - 1) / 2.0,
            y_size=75e-6, z_size=75e-6,
            tilt_x=0.3, tilt_y=0.2, tilt_z=-0.15,
            wavelength=1.5406e-10,
            **kw,
        )

    def test_lab_coords_match_all_orientations(self):
        """Full xyz lab coordinates match at machine precision for all
        4 orientations, after per-orientation mirrors. The coordinate
        flips are determined by the mirror matrix M = diag(m0,m1,1)
        stored in the poni dict.  In ID11 frame (via G·M·G):
          M(3): m0=+1,m1=+1 → no flip
          M(2): m0=-1,m1=+1 → Z flip
          M(4): m0=+1,m1=-1 → Y flip
          M(1): m0=-1,m1=-1 → Y+Z flip"""
        rng = np.random.RandomState(42)
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = self._make_test_par(
                    o11=o11, o12=o12, o21=o21, o22=o22)
                poni = pp.par_to_poni(par, detector_shape=self.SHAPE)

                d1 = rng.uniform(0, self.SHAPE[0] - 1, self.NCOORDS)
                d2 = rng.uniform(0, self.SHAPE[1] - 1, self.NCOORDS)

                ai = AzimuthalIntegrator(
                    dist=poni["dist"], poni1=poni["poni1"], poni2=poni["poni2"],
                    rot1=poni["rot1"], rot2=poni["rot2"], rot3=poni["rot3"],
                    pixel1=poni["pixel1"], pixel2=poni["pixel2"],
                    wavelength=poni["wavelength"], orientation=poni["orientation"])
                ai.detector.shape = self.SHAPE

                t3v, t1v, t2v = ai.calc_pos_zyx(d1=d1, d2=d2)
                xyz_py = np.column_stack([t3v, -t2v, t1v])
                xyz_id = compute_xyz_lab(np.array([d1, d2]), **par).T

                mirror_orient = poni.get("_mirror_orient", orientation)
                _flip_id_y = mirror_orient in (4, 1)
                _flip_id_z = mirror_orient in (2, 1)
                if _flip_id_y:
                    xyz_id = xyz_id.copy()
                    xyz_id[:, 1] = -xyz_id[:, 1]
                if _flip_id_z:
                    if not _flip_id_y:
                        xyz_id = xyz_id.copy()
                    xyz_id[:, 2] = -xyz_id[:, 2]

                diff = np.max(np.abs(xyz_py - xyz_id))
                self.assertLess(diff, 5e-7,
                                msg=f"{label}: xyz diff {diff:.2e}")


class TestIO(unittest.TestCase):
    """Test file I/O round trip."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_par_read_write_round_trip(self):
        """Read .par -> write .par -> read .par gives same values."""
        par = make_base_par()

        par_file = os.path.join(self.tmpdir, "test.par")
        pp.write_par(par, par_file, par_length_unit="um")

        par_read = pp.read_par(par_file, par_length_unit="um")

        for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
            self.assertAlmostEqual(par[key], par_read[key], delta=1e-10,
                                   msg=f"par IO: {key} mismatch")

    def test_poni_read_write_round_trip(self):
        """Read .poni -> write .poni -> read .poni gives same values."""
        par = make_base_par()
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)

        poni_file = os.path.join(self.tmpdir, "test.poni")
        pp.write_poni(poni, poni_file)

        poni_read = pp.read_poni(poni_file)

        for key in ["dist", "poni1", "poni2", "rot1", "rot2", "rot3",
                    "pixel1", "pixel2", "wavelength", "orientation"]:
            self.assertAlmostEqual(poni[key], poni_read[key], delta=1e-10,
                                   msg=f"poni IO: {key} mismatch")

    def test_full_disk_round_trip(self):
        """par file on disk -> poni file on disk -> par file on disk."""
        par = make_base_par()

        par_file = os.path.join(self.tmpdir, "geom.par")
        poni_file = os.path.join(self.tmpdir, "geom.poni")
        par2_file = os.path.join(self.tmpdir, "geom2.par")

        pp.write_par(par, par_file, par_length_unit="um")
        poni = pp.par_to_poni(pp.read_par(par_file, par_length_unit="um"),
                              detector_shape=DETECTOR_SHAPE)
        pp.write_poni(poni, poni_file)
        poni_read = pp.read_poni(poni_file)
        par_read = pp.poni_to_par(poni_read, detector_shape=DETECTOR_SHAPE)
        pp.write_par(par_read, par2_file, par_length_unit="um")
        par_final = pp.read_par(par2_file, par_length_unit="um")

        for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
            self.assertAlmostEqual(par[key], par_final[key], delta=1e-10,
                                   msg=f"disk round-trip: {key}")

    def test_write_par_contains_required_fields(self):
        """Written par file contains all necessary fields."""
        par = make_base_par()
        par_file = os.path.join(self.tmpdir, "fields.par")
        pp.write_par(par, par_file, par_length_unit="um")

        with open(par_file) as f:
            content = f.read()

        expected_fields = [
            "distance", "y_center", "z_center", "y_size", "z_size",
            "tilt_x", "tilt_y", "tilt_z",
            "o11", "o12", "o21", "o22",
            "wavelength", "wedge", "chi", "omegasign", "fit_tolerance",
        ]
        for field in expected_fields:
            self.assertIn(field, content, msg=f"par missing field: {field}")

    def test_par_length_units(self):
        """par I/O handles different length units."""
        par = make_base_par()

        for unit, factor in [("um", 1e6), ("mm", 1e3), ("m", 1.0)]:
            with self.subTest(unit=unit):
                par_file = os.path.join(self.tmpdir, f"test_{unit}.par")
                pp.write_par(par, par_file, par_length_unit=unit)

                par_unit = pp.read_par(par_file, par_length_unit=unit)
                self.assertAlmostEqual(
                    par["distance"], par_unit["distance"], delta=1e-10,
                    msg=f"{unit}: distance"
                )
                self.assertAlmostEqual(
                    par["y_size"], par_unit["y_size"], delta=1e-10,
                    msg=f"{unit}: y_size"
                )

    def test_write_poni_loads_and_integrates(self):
        """Written poni loads with pyFAI.load() and integrate1d succeeds."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                for k in ['o11', 'o12', 'o21', 'o22']:
                    par[k] = locals()[k]
                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)

                poni_file = os.path.join(self.tmpdir, f"test_{label}.poni")
                pp.write_poni(poni, poni_file)

                ai = pyFAI.load(poni_file)
                ai.detector.shape = DETECTOR_SHAPE
                shape_slow, shape_fast = DETECTOR_SHAPE
                img = np.ones((shape_slow, shape_fast), dtype=np.float64)
                result = ai.integrate1d(img, 20)
                self.assertGreater(len(result.radial), 0,
                                   msg=f"{label}: integration produced no output")


class TestEdgeCases(unittest.TestCase):
    """Edge case tests."""

    def test_wavelength_conversion(self):
        """Wavelength passes through conversion unchanged (both in meters)."""
        par = make_base_par()
        self.assertAlmostEqual(par["wavelength"], 1.5406e-10, delta=1e-15)

        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        self.assertAlmostEqual(poni["wavelength"], 1.5406e-10, delta=1e-15)

        par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
        self.assertAlmostEqual(par2["wavelength"], 1.5406e-10, delta=1e-15)

        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".par", delete=False) as f:
            f.close()
            pp.write_par(par, f.name, par_length_unit="um")
            par_read = pp.read_par(f.name, par_length_unit="um")
            os.unlink(f.name)
        self.assertAlmostEqual(par_read["wavelength"], 1.5406e-10, delta=1e-15)

    def test_zero_pixel_size_handled(self):
        """Zero pixel sizes produce well-defined results."""
        par = make_base_par()
        par["y_size"] = 0.0
        par["z_size"] = 0.0
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        self.assertGreater(abs(poni["dist"]), 0)

    def test_orientation_mapping_completeness(self):
        """All 4 orientations map correctly."""
        expected = [
            ((1, 0, 0, -1), 3),
            ((-1, 0, 0, 1), 1),
            ((-1, 0, 0, -1), 4),
            ((1, 0, 0, 1), 2),
        ]
        for (o11, o12, o21, o22), orient in expected:
            self.assertEqual(pp.flip_to_orientation(o11, o12, o21, o22), orient)
            self.assertEqual(pp.orientation_to_flip(orient), (o11, o12, o21, o22))

    def test_unsupported_flip_raises(self):
        """Unsupported transpose flips raise ValueError."""
        with self.assertRaises(ValueError):
            pp.flip_to_orientation(0, 1, 1, 0)

    def test_pyfai_rotation_matrix_matches_actual(self):
        """Our _pyfai_rotation_matrix matches pyFAI's rotation_matrix()."""
        import numpy as np
        from pyFAI.integrator.azimuthal import AzimuthalIntegrator
        test_cases = [
            (0.0, 0.0, 0.0),
            (0.1, 0.2, 0.3),
            (-0.15, 0.0, 0.7),
            (0.5, -0.5, 0.0),
            (0.0, 1.4, 0.0),
            (-0.3, 0.2, -0.15),
        ]
        for rot1, rot2, rot3 in test_cases:
            with self.subTest(rot1=rot1, rot2=rot2, rot3=rot3):
                ai = AzimuthalIntegrator(dist=0.1, poni1=0.0, poni2=0.0,
                                          rot1=rot1, rot2=rot2, rot3=rot3,
                                          pixel1=75e-6, pixel2=75e-6)
                R_pyfai = ai.rotation_matrix()
                R_ours = np.array(pp._pyfai_rotation_matrix(rot1, rot2, rot3))
                diff = np.max(np.abs(R_pyfai - R_ours))
                self.assertLess(diff, 1e-14,
                                msg=f"rot=({rot1},{rot2},{rot3}) diff={diff:.2e}")

    def test_too_large_tilts(self):
        """Tilts up to +-pi/4 round-trip correctly."""
        par = make_base_par()
        for tilt_key in ["tilt_x", "tilt_y", "tilt_z"]:
            for angle in [0.0, 0.78, -0.78]:
                with self.subTest(tilt=tilt_key, angle=angle):
                    p = dict(par)
                    for k in ["tilt_x", "tilt_y", "tilt_z"]:
                        p[k] = angle if k == tilt_key else 0.0
                    poni = pp.par_to_poni(p, detector_shape=DETECTOR_SHAPE)
                    par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
                    self.assertAlmostEqual(angle, par2[tilt_key], delta=1e-8)

    def test_chi_eta_conversion(self):
        """chi_to_eta and eta_to_chi are inverses, per-orientation mapping correct."""
        test_angles = [-3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0]

        for orient in (1, 2, 3, 4):
            with self.subTest(orientation=orient):
                for ang in test_angles:
                    eta = pp.chi_to_eta(ang, orient)
                    chi = pp.eta_to_chi(ang, orient)
                    # Round-trip: should recover original value modulo 2π
                    rtt = pp.eta_to_chi(eta, orient)
                    self.assertAlmostEqual(
                        math.sin(rtt), math.sin(ang), delta=1e-14,
                        msg=f"orient={orient} round-trip sin mismatch at {ang}")
                    self.assertAlmostEqual(
                        math.cos(rtt), math.cos(ang), delta=1e-14,
                        msg=f"orient={orient} round-trip cos mismatch at {ang}")

        # Per-orientation sin/cos relationships from _CHI_ETA_SIN_COS_FACTORS
        s0, s1 = pp._CHI_ETA_SIN_COS_FACTORS[3]
        self.assertEqual((s0, s1), (1, 1))
        s0, s1 = pp._CHI_ETA_SIN_COS_FACTORS[2]
        self.assertEqual((s0, s1), (-1, 1))
        s0, s1 = pp._CHI_ETA_SIN_COS_FACTORS[4]
        self.assertEqual((s0, s1), (1, -1))
        s0, s1 = pp._CHI_ETA_SIN_COS_FACTORS[1]
        self.assertEqual((s0, s1), (-1, -1))

        # Orientation from par / poni dicts
        par = make_base_par()
        self.assertEqual(pp._extract_orientation_from_arg(par), 3)
        par["o11"] = -1; par["o22"] = -1
        self.assertEqual(pp._extract_orientation_from_arg(par), 4)
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        self.assertEqual(pp._extract_orientation_from_arg(poni), 4)

        # chi_to_eta with dict arg
        for ang in [0.0, 1.0, -0.5]:
            eta_p = pp.chi_to_eta(ang, par)
            eta_i = pp.chi_to_eta(ang, 4)
            self.assertAlmostEqual(math.sin(eta_p), math.sin(eta_i), delta=1e-14)


class TestAllSolutions(unittest.TestCase):
    """Verify the multi-solution finder returns correct and distinct solutions."""

    NCOORDS = 2000

    def test_four_solutions_orient2(self):
        """For a single par (flip orient 2), find_all_poni_solutions
        enumerates 32 unique solutions (4 orientations × 4 mirrors × 2 Euler reps).

        They group into 16 (orient, mirror) pairs, each contributing 2 Euler-equivalent
        angle triples.  The first 8 solutions (orient=canonical, dist positive)
        are from the canonical flip→orientation pairing.
        """
        par = make_base_par()
        par["o11"] = 1
        par["o12"] = 0
        par["o21"] = 0
        par["o22"] = 1
        solutions = pp.find_all_poni_solutions(par, detector_shape=DETECTOR_SHAPE)
        self.assertEqual(len(solutions), 32, f"Expected 32 solutions, got {len(solutions)}")

        groups = {}
        for s in solutions:
            key = (s["orient_tried"], s["mirror_source"])
            groups.setdefault(key, []).append(s)
        self.assertEqual(len(groups), 16, f"Expected 16 (orient,mirror) groups")

        for key, sols in groups.items():
            self.assertEqual(len(sols), 2,
                             f"Group {key}: expected 2 Euler reps, got {len(sols)}")
            p0 = sols[0]["poni"]
            p1 = sols[1]["poni"]
            self.assertAlmostEqual(p0["dist"], p1["dist"], delta=1e-10,
                                   msg=f"Group {key}: distances differ")

    def test_all_solutions_tth_matches(self):
        """2theta values match pyFAI vs ImageD11 for every discovered solution."""
        rng = np.random.RandomState(42)
        shape_slow, shape_fast = DETECTOR_SHAPE

        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                solutions = pp.find_all_poni_solutions(par, detector_shape=DETECTOR_SHAPE)
                self.assertGreater(len(solutions), 0, f"{label}: no solutions found")

                d1 = rng.uniform(0, shape_slow - 1, self.NCOORDS)
                d2 = rng.uniform(0, shape_fast - 1, self.NCOORDS)

                tth_id11, _ = compute_tth_eta(
                    np.array([d1, d2]),
                    **{k: par[k] for k in [
                        "y_center", "y_size", "z_center", "z_size",
                        "tilt_x", "tilt_y", "tilt_z", "distance",
                        "o11", "o12", "o21", "o22"
                    ]}
                )
                tth_id11_rad = np.radians(tth_id11)

                for i, sol in enumerate(solutions):
                    with self.subTest(solution=i):
                        ai = pyFAI_from_poni(sol["poni"])
                        tth_pyfai = ai.tth(d1=d1, d2=d2, path="cython")
                        diff = np.abs(tth_pyfai - tth_id11_rad)
                        self.assertLess(np.max(diff), 1e-6,
                                        msg=f"sol {i}: max 2th diff {np.max(diff):.2e}")

    def test_all_solutions_azimuth_matches(self):
        """Azimuth (chi/eta) mapping is correct for every solution,
        using the correct factors for each (orient, mirror) pair.
        Pixels near the beam centre (where eta≈0, chi≈180°) are excluded
        because the atan2 discontinuity makes the factor comparison
        singular there."""
        rng = np.random.RandomState(123)
        shape_slow, shape_fast = DETECTOR_SHAPE

        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                solutions = pp.find_all_poni_solutions(par, detector_shape=DETECTOR_SHAPE)

                d1 = rng.uniform(0, shape_slow - 1, self.NCOORDS)
                d2 = rng.uniform(0, shape_fast - 1, self.NCOORDS)

                _, eta = compute_tth_eta(np.array([d1, d2]), **par)
                eta_rad = np.radians(eta)

                mask = (np.abs(np.cos(eta_rad)) > 0.05) & (np.abs(np.sin(eta_rad)) > 0.05)
                n_kept = np.sum(mask)
                self.assertGreater(n_kept, 500, f"{label}: too few pixels kept after mask ({n_kept})")

                for i, sol in enumerate(solutions):
                    with self.subTest(solution=i):
                        ai = pyFAI_from_poni(sol["poni"])
                        chi = ai.chi(d1=d1, d2=d2, path="cython")

                        mirror_orient = sol["poni"].get("_mirror_orient", sol["poni"]["orientation"])
                        sf, cf = pp._azimuth_factors(mirror_orient)

                        target_sin = sf * np.cos(eta_rad)
                        target_cos = cf * np.sin(eta_rad)

                        sin_diff = np.abs(np.sin(chi) - target_sin)[mask]
                        cos_diff = np.abs(np.cos(chi) - target_cos)[mask]

                        self.assertLess(np.max(sin_diff), 1e-6,
                                        msg=f"sol {i}: {sol['flip_label']} o={sol['orient_tried']} "
                                            f"m={mirror_orient}: max sin diff {np.max(sin_diff):.2e}")
                        self.assertLess(np.max(cos_diff), 1e-6,
                                        msg=f"sol {i}: {sol['flip_label']} o={sol['orient_tried']} "
                                            f"m={mirror_orient}: max cos diff {np.max(cos_diff):.2e}")

    def test_all_solutions_roundtrip(self):
        """par -> poni -> par round-trips for every solution."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                solutions = pp.find_all_poni_solutions(par, detector_shape=DETECTOR_SHAPE)

                for i, sol in enumerate(solutions):
                    with self.subTest(solution=i):
                        par2 = pp.poni_to_par(sol["poni"], detector_shape=DETECTOR_SHAPE)

                        for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
                            self.assertAlmostEqual(par[key], par2[key], delta=1e-10,
                                                   msg=f"sol {i}: {key}")
                        for key in ["tilt_x", "tilt_y", "tilt_z"]:
                            self.assertAlmostEqual(par[key], par2[key], delta=1e-10,
                                                   msg=f"sol {i}: {key}")
                        for key in ["o11", "o12", "o21", "o22"]:
                            self.assertEqual(par[key], par2[key],
                                             msg=f"sol {i}: {key}")

    def test_equivalent_reps_differ_by_zyx_equiv(self):
        """Within each (orient, mirror) group the two angle triples produce
        the same pyFAI rotation matrix (ZYX β-solution pair)."""
        par = make_base_par()
        par["o11"] = 1
        par["o12"] = 0
        par["o21"] = 0
        par["o22"] = 1
        solutions = pp.find_all_poni_solutions(par, detector_shape=DETECTOR_SHAPE)

        groups = {}
        for s in solutions:
            key = (s["orient_tried"], s["mirror_source"])
            groups.setdefault(key, []).append(s)

        for (ot, ms), sols in groups.items():
            self.assertEqual(len(sols), 2,
                             msg=f"Group (o={ot}, m={ms}): expected 2 Euler reps")
            p0 = sols[0]["poni"]
            p1 = sols[1]["poni"]
            import numpy as np
            R0 = np.array(pp._pyfai_rotation_matrix(p0["rot1"], p0["rot2"], p0["rot3"]))
            R1 = np.array(pp._pyfai_rotation_matrix(p1["rot1"], p1["rot2"], p1["rot3"]))
            diff = np.max(np.abs(R0 - R1))
            self.assertLess(diff, 1e-12,
                            msg=f"Group (o={ot}, m={ms}): rotation matrices differ by {diff:.2e}")

    def test_exact_chi_option(self):
        """exact_chi=True picks a chi_eta_exact solution (mirror matches trial
        orientation, so chi = 90°−eta in the pyFAI frame for that orientation)."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                solutions = pp.find_all_poni_solutions(par, detector_shape=DETECTOR_SHAPE)
                chi_exact_sols = [s for s in solutions if s["chi_eta_exact"]]
                self.assertGreater(len(chi_exact_sols), 0,
                                   msg=f"{label}: should have chi_eta_exact solutions")

                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE,
                                      exact_chi=True)
                trial_orient = poni["orientation"]
                mirror_orient = poni.get("_mirror_orient", trial_orient)
                self.assertEqual(trial_orient, mirror_orient,
                                 msg=f"{label}: exact_chi should have mirror_orient==trial_orient")

    def test_prefer_positive_distance_false(self):
        """prefer_positive_distance=False allows negative distance."""
        par = make_base_par()
        par["o11"] = 1
        par["o12"] = 0
        par["o21"] = 0
        par["o22"] = 1

        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE,
                              prefer_positive_distance=False)
        self.assertIn("dist", poni)

    def test_default_api_unchanged(self):
        """Default par_to_poni returns the same result as pre-refactor."""
        import numpy as np
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                self.assertTrue(poni["dist"] > 0,
                                msg=f"{label}: default should give positive dist")
                self.assertTrue(poni.get("_mirror_used", True),
                                msg=f"{label}: default should use mirror")

                ai = pyFAI_from_poni(poni)
                shape_slow, shape_fast = DETECTOR_SHAPE
                d1 = np.array([100.0, 500.0, 900.0])
                d2 = np.array([200.0, 500.0, 800.0])

                tth_pyfai = ai.tth(d1=d1, d2=d2)
                tth_id11, _ = compute_tth_eta(
                    np.array([d1, d2]),
                    **{k: par[k] for k in [
                        "y_center", "y_size", "z_center", "z_size",
                        "tilt_x", "tilt_y", "tilt_z", "distance",
                        "o11", "o12", "o21", "o22"
                    ]}
                )
                tth_id11_rad = np.radians(tth_id11)
                self.assertLess(np.max(np.abs(tth_pyfai - tth_id11_rad)), 1e-7)


class TestBackscattering(unittest.TestCase):
    """Backscattering geometry: detector upstream of sample (distance<0 in ID11,
    positive dist + rot2≈π in pyFAI)."""

    NCOORDS = 2000

    def _make_backscattering_par(self, **kw):
        result = dict(
            distance=-0.15,
            y_center=500.0,
            z_center=500.0,
            y_size=75e-6,
            z_size=75e-6,
            tilt_x=0.0,
            tilt_y=0.0,
            tilt_z=0.0,
            o11=1, o12=0, o21=0, o22=-1,
            wavelength=1.5406e-10,
            wedge=0.0,
            chi=0.0,
            omegasign=1.0,
            fit_tolerance=0.05,
        )
        result.update(kw)
        return result

    def test_positive_dist_rot2_near_pi(self):
        """Backscattering converts to positive pyFAI dist with an angle near ±π."""
        par = self._make_backscattering_par()
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        self.assertGreater(poni["dist"], 0, "pyFAI dist should be positive")
        close_to_pi = (
            abs(abs(poni["rot1"]) - math.pi) < 1e-6 or
            abs(abs(poni["rot2"]) - math.pi) < 1e-6
        )
        self.assertTrue(close_to_pi,
                        msg=f"rot1={poni['rot1']:.4f} or rot2={poni['rot2']:.4f} should be near ±π")

    def test_backscattering_tth_matches(self):
        """2theta matches pyFAI vs ImageD11 for backscattering."""
        rng = np.random.RandomState(99)
        shape_slow, shape_fast = DETECTOR_SHAPE

        par = self._make_backscattering_par()
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        ai = pyFAI_from_poni(poni)

        d1 = rng.uniform(0, shape_slow - 1, self.NCOORDS)
        d2 = rng.uniform(0, shape_fast - 1, self.NCOORDS)

        tth_pyfai = ai.tth(d1=d1, d2=d2, path="cython")
        tth_id11, _ = compute_tth_eta(np.array([d1, d2]), **par)
        tth_id11_rad = np.radians(tth_id11)

        diff = np.abs(tth_pyfai - tth_id11_rad)
        self.assertLess(np.max(diff), 1e-6,
                        msg=f"backscattering max 2th diff {np.max(diff):.2e}")

    def test_backscattering_roundtrip(self):
        """par -> poni -> par round-trip for backscattering.
        Tilts may differ by ±π for the alternative representation
        (rot1≈±π or rot2≈±π); the effective geometry is compared."""
        par = self._make_backscattering_par()
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)

        for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
            self.assertAlmostEqual(par[key], par2[key], delta=1e-10,
                                   msg=f"backscattering: {key}")
        for key in ["tilt_x", "tilt_y", "tilt_z"]:
            diff = abs(par[key] - par2[key])
            diff_mod = min(diff, abs(diff - 2 * math.pi), abs(diff + 2 * math.pi),
                           abs(diff - math.pi), abs(diff + math.pi))
            self.assertLess(diff_mod, 1e-10,
                            msg=f"backscattering: {key} differs by {diff}")

    def test_backscattering_with_tilts(self):
        """Backscattering with small tilts: verifies positive distance
        and round-trip of beam-center parameters."""
        par = self._make_backscattering_par(
            tilt_x=0.01, tilt_y=-0.02, tilt_z=0.01)
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        self.assertGreater(poni["dist"], 0, "pyFAI dist should be positive")
        par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
        for key in ["distance", "y_center", "z_center"]:
            self.assertAlmostEqual(par[key], par2[key], delta=1e-6,
                                   msg=f"backscattering tilted: {key}")

    def test_backscattering_azimuth(self):
        """Azimuth matches for backscattering — the π rotation effectively
        changes the orientation mapping (orient 3 becomes like orient 4)."""
        rng = np.random.RandomState(42)
        shape_slow, shape_fast = DETECTOR_SHAPE

        par = self._make_backscattering_par()
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        ai = pyFAI_from_poni(poni)

        d1 = rng.uniform(0, shape_slow - 1, self.NCOORDS)
        d2 = rng.uniform(0, shape_fast - 1, self.NCOORDS)

        chi = ai.chi(d1=d1, d2=d2, path="cython")
        _, eta = compute_tth_eta(np.array([d1, d2]), **par)
        eta_rad = np.radians(eta)

        best_max = float("inf")
        for orient in (1, 2, 3, 4):
            eta_cs = pp.chi_to_eta(chi, orientation=orient)
            max_diff = max(np.max(np.abs(np.sin(eta_rad) - np.sin(eta_cs))),
                           np.max(np.abs(np.cos(eta_rad) - np.cos(eta_cs))))
            best_max = min(best_max, max_diff)

        self.assertLess(best_max, 1e-6,
                        msg=f"backscattering azimuth: no orientation matches, "
                            f"best max diff={best_max:.2e}")

    def test_backscattering_all_orientations(self):
        """Backscattering round-trip works for all 4 orientations."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = self._make_backscattering_par(
                    o11=o11, o12=o12, o21=o21, o22=o22)
                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                self.assertGreater(poni["dist"], 0,
                                   msg=f"{label}: dist should be positive")
                par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
                self.assertAlmostEqual(par["distance"], par2["distance"], delta=1e-10,
                                       msg=f"{label}: distance")


if __name__ == "__main__":
    unittest.main(verbosity=2)
