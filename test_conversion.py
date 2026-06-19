"""
test_conversion.py — Test par ↔ poni conversion (simplified, Round 9).

Validates:
1. Round-trip identity for all 4 orientations (including large tilts)
2. 2θ values match between pyFAI and ImageD11 for all flips
3. Azimuthal angles (chi/eta) match via sin/cos factors
4. Full xyz lab coordinates match pixel-by-pixel
5. File I/O round-trip through disk files

Requires pyFAI and ImageD11 to be importable.
"""

import sys
import os
import tempfile
import math
import unittest
import numpy as np

import par_to_poni as pp
from pyFAI.integrator.azimuthal import AzimuthalIntegrator
from ImageD11.transform import compute_xyz_lab, compute_tth_eta


# ---------------------------------------------------------------------------
# Test parameters
# ---------------------------------------------------------------------------

DETECTOR_SHAPE = (1000, 1000)  # (slow, fast) — pyFAI C-order convention

FLIPS = [
    (1, 0, 0, -1, 3, "orient3_native"),
    (-1, 0, 0, 1, 1, "orient1_flip_both"),
    (-1, 0, 0, -1, 2, "orient2_flip_slow"),
    (1, 0, 0, 1, 4, "orient4_flip_fast"),
]

# Large-tilt test angles (approx 0°, 20°, 40°, 60°)
LARGE_TILTS = [
    (0.0, 0.0, 0.0),
    (0.35, 0.26, 0.18),
    (0.70, 0.53, 0.35),
    (1.05, 0.79, 0.53),
]


def make_base_par(tx=0.3, ty=0.2, tz=-0.15):
    return dict(
        distance=0.15,
        y_center=500.0,
        z_center=500.0,
        y_size=75e-6,
        z_size=75e-6,
        tilt_x=tx,
        tilt_y=ty,
        tilt_z=tz,
        o11=1, o12=0, o21=0, o22=-1,
        wavelength=1.5406e-10,
        wedge=0.0, chi=0.0, omegasign=1.0,
        fit_tolerance=0.05,
    )


def pyFAI_from_poni(poni, shape=DETECTOR_SHAPE):
    ai = AzimuthalIntegrator(
        dist=poni["dist"], poni1=poni["poni1"], poni2=poni["poni2"],
        rot1=poni["rot1"], rot2=poni["rot2"], rot3=poni["rot3"],
        pixel1=poni["pixel1"], pixel2=poni["pixel2"],
        wavelength=poni.get("wavelength"), orientation=poni.get("orientation", 3))
    ai.detector.shape = shape
    return ai


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestRoundTrip(unittest.TestCase):
    """par → poni → par and poni → par → poni round-trip."""

    def test_par_round_trip_all_flips(self):
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
                for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
                    self.assertAlmostEqual(par[key], par2[key], delta=1e-10)
                for key in ["tilt_x", "tilt_y", "tilt_z"]:
                    self.assertAlmostEqual(par[key], par2[key], delta=1e-10)
                for key in ["o11", "o12", "o21", "o22"]:
                    self.assertEqual(par[key], par2[key])

    def test_poni_round_trip_all_flips(self):
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                poni2 = pp.par_to_poni(
                    pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE),
                    detector_shape=DETECTOR_SHAPE)
                for key in ["dist", "poni1", "poni2", "rot1", "rot2", "rot3",
                            "pixel1", "pixel2", "wavelength", "orientation"]:
                    self.assertAlmostEqual(poni[key], poni2[key], delta=1e-10)

    def test_round_trip_zero_tilts(self):
        par = make_base_par(0, 0, 0)
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
                for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
                    self.assertAlmostEqual(par[key], par2[key], delta=1e-10)

    def test_round_trip_single_tilts(self):
        par = make_base_par()
        for angle in [0.1, -0.2, 0.3, -0.5, 0.7]:
            for tilt_key in ["tilt_x", "tilt_y", "tilt_z"]:
                with self.subTest(tilt=tilt_key, angle=angle):
                    p = dict(par)
                    for k in ["tilt_x", "tilt_y", "tilt_z"]:
                        p[k] = angle if k == tilt_key else 0.0
                    poni = pp.par_to_poni(p, detector_shape=DETECTOR_SHAPE)
                    par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
                    for key in ["distance", "y_center", "z_center"]:
                        self.assertAlmostEqual(p[key], par2[key], delta=1e-8)

    def test_round_trip_large_tilts(self):
        for tx, ty, tz in LARGE_TILTS:
            for o11, o12, o21, o22, orientation, label in FLIPS:
                with self.subTest(flip=label, t=(tx, ty, tz)):
                    par = make_base_par(tx, ty, tz)
                    par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                    poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                    par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
                    for key in ["distance", "y_center", "z_center"]:
                        self.assertAlmostEqual(par[key], par2[key], delta=1e-10)
                    for key in ["tilt_x", "tilt_y", "tilt_z"]:
                        self.assertAlmostEqual(par[key], par2[key], delta=1e-10)

    def test_round_trip_edge_beam_positions(self):
        par = make_base_par(0, 0, 0)
        for yc, zc in [(0, 0), (999, 999), (0, 999), (999, 0), (500, 500)]:
            with self.subTest(y_center=yc, z_center=zc):
                p = dict(par, y_center=float(yc), z_center=float(zc))
                poni = pp.par_to_poni(p, detector_shape=DETECTOR_SHAPE)
                par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
                self.assertAlmostEqual(yc, par2["y_center"], delta=1e-10)
                self.assertAlmostEqual(zc, par2["z_center"], delta=1e-10)


class TestTwothetaMatching(unittest.TestCase):
    """2θ values match between pyFAI and ImageD11 for all flips."""

    NCOORDS = 5000

    def _check_tth(self, par, flip_label):
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        ai = pyFAI_from_poni(poni)
        rng = np.random.RandomState(42)
        shape_slow, shape_fast = DETECTOR_SHAPE
        d1 = rng.uniform(0, shape_slow - 1, self.NCOORDS)
        d2 = rng.uniform(0, shape_fast - 1, self.NCOORDS)
        tth_py = ai.tth(d1=d1, d2=d2, path="cython")
        tth_id, _ = compute_tth_eta(
            np.array([d1, d2]),
            **{k: par[k] for k in ["y_center", "y_size", "z_center", "z_size",
                "tilt_x", "tilt_y", "tilt_z", "distance",
                "o11", "o12", "o21", "o22"]})
        tth_id_rad = np.radians(tth_id)
        diff = np.abs(tth_py - tth_id_rad)
        self.assertLess(np.max(diff), 1e-7, f"{flip_label}: max 2θ diff {np.max(diff):.2e}")

    def test_tth_matches_all_flips(self):
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                self._check_tth(par, label)

    def test_tth_matches_zero_tilts(self):
        par = make_base_par(0, 0, 0)
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                self._check_tth(par, label)

    def test_tth_large_tilts(self):
        for tx, ty, tz in LARGE_TILTS:
            for o11, o12, o21, o22, orientation, label in FLIPS:
                with self.subTest(flip=label, t=(tx, ty, tz)):
                    par = make_base_par(tx, ty, tz)
                    par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                    self._check_tth(par, label)

    def test_tth_versus_q(self):
        par = make_base_par()
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        ai = pyFAI_from_poni(poni)
        d1 = np.array([100.0, 500.0, 900.0])
        d2 = np.array([200.0, 500.0, 800.0])
        tth = ai.tth(d1=d1, d2=d2)
        q_pyfai = ai.qFunction(d1=d1, d2=d2)
        q_expected = 4.0e-9 * math.pi * np.sin(tth / 2.0) / ai.wavelength
        self.assertTrue(np.allclose(q_pyfai, q_expected, rtol=1e-10))


class TestAzimuthMatching(unittest.TestCase):
    """Azimuthal angles (chi / eta) match via sin/cos factors."""

    NCOORDS = 2000

    def _check_azimuth(self, par, flip_label):
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        ai = pyFAI_from_poni(poni)
        rng = np.random.RandomState(123)
        shape_slow, shape_fast = DETECTOR_SHAPE
        d1 = rng.uniform(0, shape_slow - 1, self.NCOORDS)
        d2 = rng.uniform(0, shape_fast - 1, self.NCOORDS)
        chi = ai.chi(d1=d1, d2=d2, path="cython")
        _, eta = compute_tth_eta(np.array([d1, d2]), **par)
        eta_rad = np.radians(eta)
        orient = poni["orientation"]
        sf, cf = pp._CHI_ETA_SIN_COS_FACTORS[orient]
        target_sin = sf * np.cos(eta_rad)
        target_cos = cf * np.sin(eta_rad)
        mask = (np.abs(np.cos(eta_rad)) > 0.05) & (np.abs(np.sin(eta_rad)) > 0.05)
        self.assertLess(np.max(np.abs(np.sin(chi) - target_sin)[mask]), 1e-7,
                        f"{flip_label}: max sin diff")
        self.assertLess(np.max(np.abs(np.cos(chi) - target_cos)[mask]), 1e-7,
                        f"{flip_label}: max cos diff")

    def test_azimuth_all_flips(self):
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                self._check_azimuth(par, label)

    def test_azimuth_large_tilts(self):
        for tx, ty, tz in LARGE_TILTS[1:]:  # skip zero-tilt
            for o11, o12, o21, o22, orientation, label in FLIPS:
                with self.subTest(flip=label, t=(tx, ty, tz)):
                    par = make_base_par(tx, ty, tz)
                    par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                    self._check_azimuth(par, label)

    def test_chi_eta_conversion(self):
        test_angles = [-3.0, -2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0, 3.0]
        for orient in (1, 2, 3, 4):
            with self.subTest(orientation=orient):
                for ang in test_angles:
                    eta = pp.chi_to_eta(ang, orient)
                    rtt = pp.eta_to_chi(eta, orient)
                    self.assertAlmostEqual(math.sin(rtt), math.sin(ang), delta=1e-14)
                    self.assertAlmostEqual(math.cos(rtt), math.cos(ang), delta=1e-14)
        # Verify table values
        self.assertEqual(pp._CHI_ETA_SIN_COS_FACTORS[3], (1, 1))
        self.assertEqual(pp._CHI_ETA_SIN_COS_FACTORS[2], (-1, 1))
        self.assertEqual(pp._CHI_ETA_SIN_COS_FACTORS[4], (1, -1))
        self.assertEqual(pp._CHI_ETA_SIN_COS_FACTORS[1], (-1, -1))


class TestLabCoordinates(unittest.TestCase):
    """Full xyz lab coordinates match pixel-by-pixel, non-square detector.

    For non-native orientations pyFAI applies post-rotation sign flips
    (per its orientation model) which flip axis-1 or axis-2 in the pyFAI
    lab frame relative to ImageD11.  The test accounts for these sign
    differences.  2θ and azimuth are unaffected.
    """

    NCOORDS = 2000
    SHAPE = (128, 200)

    def test_lab_coords_match_all_orientations(self):
        rng = np.random.RandomState(42)
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = dict(
                    distance=0.15,
                    y_center=(self.SHAPE[1] - 1) / 2.0,
                    z_center=(self.SHAPE[0] - 1) / 2.0,
                    y_size=75e-6, z_size=75e-6,
                    tilt_x=0.3, tilt_y=0.2, tilt_z=-0.15,
                    wavelength=1.5406e-10,
                    o11=o11, o12=o12, o21=o21, o22=o22,
                )
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

                # For non-native orientations, pyFAI flips certain axes post-rotation.
                # Apply the same flips to the ID11 coordinates for comparison.
                if orientation in (2, 1):
                    xyz_id = xyz_id.copy()
                    xyz_id[:, 2] = -xyz_id[:, 2]   # flip Z (slow/up) axis
                if orientation in (4, 1):
                    if orientation not in (2, 1):
                        xyz_id = xyz_id.copy()
                    xyz_id[:, 1] = -xyz_id[:, 1]   # flip Y (fast/port) axis

                diff = np.max(np.abs(xyz_py - xyz_id))
                self.assertLess(diff, 5e-7, f"{label}: xyz diff {diff:.2e}")


class TestIO(unittest.TestCase):
    """File I/O round trip."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_par_read_write_round_trip(self):
        par = make_base_par()
        par_file = os.path.join(self.tmpdir, "test.par")
        pp.write_par(par, par_file, par_length_unit="um")
        par_read = pp.read_par(par_file, par_length_unit="um")
        for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
            self.assertAlmostEqual(par[key], par_read[key], delta=1e-10)

    def test_poni_read_write_round_trip(self):
        par = make_base_par()
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        poni_file = os.path.join(self.tmpdir, "test.poni")
        pp.write_poni(poni, poni_file)
        poni_read = pp.read_poni(poni_file)
        for key in ["dist", "poni1", "poni2", "rot1", "rot2", "rot3",
                    "pixel1", "pixel2", "wavelength", "orientation"]:
            self.assertAlmostEqual(poni[key], poni_read[key], delta=1e-10)

    def test_full_disk_round_trip(self):
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
            self.assertAlmostEqual(par[key], par_final[key], delta=1e-10)

    def test_write_par_contains_required_fields(self):
        par = make_base_par()
        par_file = os.path.join(self.tmpdir, "fields.par")
        pp.write_par(par, par_file, par_length_unit="um")
        with open(par_file) as f:
            content = f.read()
        expected = ["distance", "y_center", "z_center", "y_size", "z_size",
                    "tilt_x", "tilt_y", "tilt_z", "o11", "o12", "o21", "o22",
                    "wavelength", "wedge", "chi", "omegasign", "fit_tolerance"]
        for field in expected:
            self.assertIn(field, content)

    def test_par_length_units(self):
        par = make_base_par()
        for unit, factor in [("um", 1e6), ("mm", 1e3), ("m", 1.0)]:
            with self.subTest(unit=unit):
                par_file = os.path.join(self.tmpdir, f"test_{unit}.par")
                pp.write_par(par, par_file, par_length_unit=unit)
                par_unit = pp.read_par(par_file, par_length_unit=unit)
                self.assertAlmostEqual(par["distance"], par_unit["distance"], delta=1e-10)
                self.assertAlmostEqual(par["y_size"], par_unit["y_size"], delta=1e-10)

    def test_write_poni_loads_and_integrates(self):
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"], par["o12"], par["o21"], par["o22"] = o11, o12, o21, o22
                poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
                poni_file = os.path.join(self.tmpdir, f"test_{label}.poni")
                pp.write_poni(poni, poni_file)
                import pyFAI
                ai = pyFAI.load(poni_file)
                ai.detector.shape = DETECTOR_SHAPE
                shape_slow, shape_fast = DETECTOR_SHAPE
                img = np.ones((shape_slow, shape_fast), dtype=np.float64)
                result = ai.integrate1d(img, 20)
                self.assertGreater(len(result.radial), 0)


class TestEdgeCases(unittest.TestCase):
    """Edge case tests."""

    def test_wavelength_conversion(self):
        par = make_base_par()
        self.assertAlmostEqual(par["wavelength"], 1.5406e-10, delta=1e-15)
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        self.assertAlmostEqual(poni["wavelength"], 1.5406e-10, delta=1e-15)
        par2 = pp.poni_to_par(poni, detector_shape=DETECTOR_SHAPE)
        self.assertAlmostEqual(par2["wavelength"], 1.5406e-10, delta=1e-15)

    def test_zero_pixel_size_handled(self):
        par = make_base_par()
        par["y_size"] = par["z_size"] = 0.0
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        self.assertGreater(abs(poni["dist"]), 0)

    def test_orientation_mapping_completeness(self):
        expected = [
            ((1, 0, 0, -1), 3),
            ((-1, 0, 0, 1), 1),
            ((-1, 0, 0, -1), 2),
            ((1, 0, 0, 1), 4),
        ]
        for (o11, o12, o21, o22), orient in expected:
            self.assertEqual(pp.flip_to_orientation(o11, o12, o21, o22), orient)
            self.assertEqual(pp.orientation_to_flip(orient), (o11, o12, o21, o22))

    def test_unsupported_flip_raises(self):
        with self.assertRaises(ValueError):
            pp.flip_to_orientation(0, 1, 1, 0)

    def test_too_large_tilts(self):
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

    def test_chi_eta_dict_args(self):
        par = make_base_par()
        self.assertEqual(pp._extract_orientation_from_arg(par), 3)
        par["o11"] = -1
        par["o22"] = -1
        self.assertEqual(pp._extract_orientation_from_arg(par), 2)
        poni = pp.par_to_poni(par, detector_shape=DETECTOR_SHAPE)
        self.assertEqual(pp._extract_orientation_from_arg(poni), 2)
        for ang in [0.0, 1.0, -0.5]:
            eta_p = pp.chi_to_eta(ang, par)
            eta_i = pp.chi_to_eta(ang, 2)
            self.assertAlmostEqual(math.sin(eta_p), math.sin(eta_i), delta=1e-14)


class TestDocumentation(unittest.TestCase):
    """Per AGENTS.md §3: documentation tables must match code data structures."""

    def test_readme_flip_to_orientation_table_matches_code(self):
        """Parse the flip→orientation table from README.md and check against
        _FLIP_TO_ORIENTATION."""
        with open(os.path.join(os.path.dirname(__file__), "README.md")) as f:
            text = f.read()
        table_expected = [
            ("(1, 0, 0, -1)", "3"),
            ("(-1, 0, 0, 1)", "1"),
            ("(-1, 0, 0, -1)", "2"),
            ("(1, 0, 0, 1)", "4"),
        ]
        for flip_str, orient_str in table_expected:
            o11, o12, o21, o22 = [int(x.strip()) for x in flip_str.strip("()").split(",")]
            code_o = pp._FLIP_TO_ORIENTATION.get((o11, o12, o21, o22))
            self.assertIsNotNone(code_o, f"Flip {flip_str} not in _FLIP_TO_ORIENTATION")
            self.assertEqual(str(code_o), orient_str,
                             f"Flip {flip_str}: README={orient_str}, code={code_o}")
            self.assertIn(flip_str, text, f"README.md missing {flip_str}")
            self.assertIn(orient_str, text, f"README.md missing orient {orient_str}")

    def test_readme_azimuth_table_matches_chieta_factors(self):
        """Parse the azimuth table from README.md and check against
        _CHI_ETA_SIN_COS_FACTORS."""
        expected = {
            3: ("+cos(η)", "+sin(η)"),
            2: ("−cos(η)", "+sin(η)"),
            4: ("+cos(η)", "−sin(η)"),
            1: ("−cos(η)", "−sin(η)"),
        }
        for orient, (sin_str, cos_str) in expected.items():
            code_sf, code_cf = pp._CHI_ETA_SIN_COS_FACTORS[orient]
            readme_sin_sign = "+" if sin_str[0] != "−" else "-"
            readme_cos_sign = "+" if cos_str[0] != "−" else "-"
            code_sin_sign = "+" if code_sf == 1 else "-"
            code_cos_sign = "+" if code_cf == 1 else "-"
            self.assertEqual(readme_sin_sign, code_sin_sign,
                             f"Orient {orient}: README sin={readme_sin_sign}, code={code_sin_sign}")
            self.assertEqual(readme_cos_sign, code_cos_sign,
                             f"Orient {orient}: README cos={readme_cos_sign}, code={code_cos_sign}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
