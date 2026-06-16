"""
test_conversion.py — Test the par ↔ poni conversion.

Validates:
1. Round-trip identity for all 4 orientations
2. 2θ / q values match between pyFAI and ImageD11
3. Azimuthal angles match (sin/cos comparison)
4. IO round-trip through disk files

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
from ImageD11.transform import (
    compute_xyz_lab,
    compute_tth_eta,
    detector_rotation_matrix,
)
from ImageD11.parameters import parameters as ImageD11Parameters


# ---------------------------------------------------------------------------
# Test parameters
# ---------------------------------------------------------------------------

def make_base_par():
    """Create base par dict with the specified strongly-tilted geometry."""
    return dict(
        distance=0.15,           # m (internal)
        y_center=500.0,          # px
        z_center=500.0,          # px
        y_size=75e-6,            # m/px (75 µm)
        z_size=75e-6,            # m/px (75 µm)
        tilt_x=0.3,              # rad (~17°)
        tilt_y=0.2,              # rad (~11°)
        tilt_z=-0.15,            # rad (~9°)
        o11=1, o12=0, o21=0, o22=-1,    # orientation 3
        wavelength=1.5406e-10,   # m (Cu Kα)
        wedge=0.0,
        chi=0.0,
        omegasign=1.0,
        fit_tolerance=0.05,
    )


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
    ai.detector.shape = (1000, 1000)
    return ai


def pixel_lut_from_par(par, shape=(50, 50)):
    """Create an ImageD11 PixelLUT from a par dict."""
    from ImageD11.transform import PixelLUT
    p = dict(par)
    p["shape"] = shape
    return PixelLUT(p)


def arctan2_sincos_pair(angle_rad):
    """Return (sin, cos) for angle comparison, avoiding wrap issues."""
    return math.sin(angle_rad), math.cos(angle_rad)


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

class TestRoundTrip(unittest.TestCase):
    """Test par → poni → par and poni → par → poni round trips."""

    def test_par_round_trip_all_flips(self):
        """par → poni → par should recover original values."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                poni = pp.par_to_poni(par)
                par2 = pp.poni_to_par(poni)

                # Check geometry fields match within tolerance
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
        """poni → par → poni should recover original values."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                poni = pp.par_to_poni(par)
                poni2 = pp.par_to_poni(pp.poni_to_par(poni))

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

                poni = pp.par_to_poni(par)
                par2 = pp.poni_to_par(poni)

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

                    poni = pp.par_to_poni(p)
                    par2 = pp.poni_to_par(poni)

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
                poni = pp.par_to_poni(p)
                par2 = pp.poni_to_par(poni)
                self.assertAlmostEqual(yc, par2["y_center"], delta=1e-10)
                self.assertAlmostEqual(zc, par2["z_center"], delta=1e-10)


class TestTwothetaMatching(unittest.TestCase):
    """Test that 2θ values match between pyFAI and ImageD11."""

    NCOORDS = 5000  # number of random test points

    def _compute_both_tth(self, par, orientation):
        """Compute 2θ from both codes using matching coordinates.

        PyFAI's orientation reorders pixel indices internally. ImageD11 uses
        raw pixel coords without reordering. So to compare at the same
        physical pixel, we pre-flip coordinates for ImageD11 based on
        pyFAI's orientation.
        """
        poni = pp.par_to_poni(par)
        self.assertEqual(poni["orientation"], orientation)

        ai = pyFAI_from_poni(poni)

        rng = np.random.RandomState(42)
        d1 = rng.uniform(0, 999, self.NCOORDS)  # slow (vertical)
        d2 = rng.uniform(0, 999, self.NCOORDS)  # fast (horizontal)

        # pyFAI: adds 0.5 internally, applies orientation-specific pixel reorder
        tth_pyfai = ai.tth(d1=d1, d2=d2, path="cython")

        # ImageD11: sc/fc used directly. Pre-flip to match pyFAI's reorder.
        # orientation 3 (native): sc=d1, fc=d2  (no reorder)
        # orientation 2 (flip slow): sc=max-1-d1, fc=d2
        # orientation 4 (flip fast): sc=d1, fc=max-1-d2
        # orientation 1 (flip both): sc=max-1-d1, fc=max-1-d2
        max_idx = 999.0
        if orientation in (2,):
            sc = max_idx - d1
            fc = d2.copy()
        elif orientation in (4,):
            sc = d1.copy()
            fc = max_idx - d2
        elif orientation in (1,):
            sc = max_idx - d1
            fc = max_idx - d2
        else:  # orientation 3
            sc = d1.copy()
            fc = d2.copy()

        tth_id11, _ = compute_tth_eta(
            np.array([sc, fc]),
            **{k: par[k] for k in [
                "y_center", "y_size", "z_center", "z_size",
                "tilt_x", "tilt_y", "tilt_z", "distance",
                "o11", "o12", "o21", "o22"
            ]}
        )
        # compute_tth_eta returns degrees, convert to radians
        tth_id11_rad = np.radians(tth_id11)

        return tth_pyfai, tth_id11_rad, d1, d2

    def test_tth_matches_all_flips(self):
        """2θ values match to machine precision for all 4 orientations."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                tth_pyfai, tth_id11, _, _ = self._compute_both_tth(par, orientation)

                diff = np.abs(tth_pyfai - tth_id11)
                max_diff = np.max(diff)

                self.assertLess(
                    max_diff, 1e-7,
                    msg=f"{label}: max 2θ diff {max_diff:.2e} exceeds 1e-7 rad"
                )

    def test_tth_matches_zero_tilts(self):
        """2θ values match when all tilts are zero."""
        par = make_base_par()
        par["tilt_x"] = par["tilt_y"] = par["tilt_z"] = 0.0

        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                tth_pyfai, tth_id11, _, _ = self._compute_both_tth(par, orientation)
                diff = np.abs(tth_pyfai - tth_id11)
                self.assertLess(np.max(diff), 1e-7, msg=f"{label}: zero tilt mismatch")

    def test_tth_versus_q(self):
        """q vector values are consistent with 2θ."""
        par = make_base_par()
        poni = pp.par_to_poni(par)
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

    def _compute_both_azimuth(self, par, orientation):
        """Compute azimuth from both codes using matching coordinates."""
        poni = pp.par_to_poni(par)
        self.assertEqual(poni["orientation"], orientation)

        ai = pyFAI_from_poni(poni)

        rng = np.random.RandomState(123)
        d1 = rng.uniform(0, 999, self.NCOORDS)
        d2 = rng.uniform(0, 999, self.NCOORDS)

        # pyFAI chi (radians)
        chi = ai.chi(d1=d1, d2=d2, path="cython")

        # Pre-flip coordinates for ImageD11 to match pyFAI's orientation reordering
        max_idx = 999.0
        if orientation in (2,):
            sc = max_idx - d1
            fc = d2.copy()
        elif orientation in (4,):
            sc = d1.copy()
            fc = max_idx - d2
        elif orientation in (1,):
            sc = max_idx - d1
            fc = max_idx - d2
        else:  # orientation 3
            sc = d1.copy()
            fc = d2.copy()

        _, eta = compute_tth_eta(
            np.array([sc, fc]),
            **{k: par[k] for k in [
                "y_center", "y_size", "z_center", "z_size",
                "tilt_x", "tilt_y", "tilt_z", "distance",
                "o11", "o12", "o21", "o22"
            ]}
        )
        eta_rad = np.radians(eta)  # compute_tth_eta returns degrees

        return chi, eta_rad

    def test_azimuth_relationship_all_flips(self):
        """Test chi = 90° - eta using sin/cos comparison — exact."""
        for o11, o12, o21, o22, orientation, label in FLIPS:
            with self.subTest(flip=label):
                par = make_base_par()
                par["o11"] = o11
                par["o12"] = o12
                par["o21"] = o21
                par["o22"] = o22

                chi, eta = self._compute_both_azimuth(par, orientation)

                target_sin = np.cos(eta)
                target_cos = np.sin(eta)

                sin_diff = np.abs(np.sin(chi) - target_sin)
                cos_diff = np.abs(np.cos(chi) - target_cos)

                self.assertLess(np.max(sin_diff), 1e-7,
                                msg=f"{label}: max sin diff {np.max(sin_diff):.2e}")
                self.assertLess(np.max(cos_diff), 1e-7,
                                msg=f"{label}: max cos diff {np.max(cos_diff):.2e}")


class TestIO(unittest.TestCase):
    """Test file I/O round trip."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp()

    def tearDown(self):
        import shutil
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_par_read_write_round_trip(self):
        """Read .par → write .par → read .par gives same values."""
        par = make_base_par()

        par_um = dict(par)
        for k in ["distance", "y_size", "z_size"]:
            par_um[k] = par[k] * 1e6   # m → µm

        par_file = os.path.join(self.tmpdir, "test.par")
        pp.write_par(par, par_file, par_length_unit="um")

        par_read = pp.read_par(par_file, par_length_unit="um")

        for key in ["distance", "y_center", "z_center", "y_size", "z_size"]:
            self.assertAlmostEqual(par[key], par_read[key], delta=1e-10,
                                   msg=f"par IO: {key} mismatch")

    def test_poni_read_write_round_trip(self):
        """Read .poni → write .poni → read .poni gives same values."""
        par = make_base_par()
        poni = pp.par_to_poni(par)

        poni_file = os.path.join(self.tmpdir, "test.poni")
        pp.write_poni(poni, poni_file)

        poni_read = pp.read_poni(poni_file)

        for key in ["dist", "poni1", "poni2", "rot1", "rot2", "rot3",
                    "pixel1", "pixel2", "wavelength", "orientation"]:
            self.assertAlmostEqual(poni[key], poni_read[key], delta=1e-10,
                                   msg=f"poni IO: {key} mismatch")

    def test_full_disk_round_trip(self):
        """par file on disk → poni file on disk → par file on disk."""
        par = make_base_par()

        par_file = os.path.join(self.tmpdir, "geom.par")
        poni_file = os.path.join(self.tmpdir, "geom.poni")
        par2_file = os.path.join(self.tmpdir, "geom2.par")

        pp.write_par(par, par_file, par_length_unit="um")
        poni = pp.par_to_poni(pp.read_par(par_file, par_length_unit="um"))
        pp.write_poni(poni, poni_file)
        poni_read = pp.read_poni(poni_file)
        par_read = pp.poni_to_par(poni_read)
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

                with open(par_file) as f:
                    content = f.read()

                # Verify the written values are in the right units
                par_unit = pp.read_par(par_file, par_length_unit=unit)
                self.assertAlmostEqual(
                    par["distance"], par_unit["distance"], delta=1e-10,
                    msg=f"{unit}: distance"
                )
                self.assertAlmostEqual(
                    par["y_size"], par_unit["y_size"], delta=1e-10,
                    msg=f"{unit}: y_size"
                )

                # Spot check: the written file should have values ≈ par * factor
                # distance was 0.15 m
                expected_distance = par["distance"] * factor
                self.assertIn(str(expected_distance)[:6], content,
                              msg=f"{unit}: expected distance value not found")


class TestEdgeCases(unittest.TestCase):
    """Edge case tests."""

    def test_wavelength_conversion(self):
        """Wavelength passes through conversion unchanged (both in meters)."""
        par = make_base_par()
        # 1.5406e-10 m = 1.5406 Å (Cu Kα) — stored internally in meters
        self.assertAlmostEqual(par["wavelength"], 1.5406e-10, delta=1e-15)

        poni = pp.par_to_poni(par)
        # Same value in meters
        self.assertAlmostEqual(poni["wavelength"], 1.5406e-10, delta=1e-15)

        par2 = pp.poni_to_par(poni)
        # Same value in meters (Å <-> m conversion happens in IO layer only)
        self.assertAlmostEqual(par2["wavelength"], 1.5406e-10, delta=1e-15)

        # Verify IO layer converts correctly: write .par in Å, read back in m
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".par", delete=False) as f:
            f.close()
            pp.write_par(par, f.name, par_length_unit="um")
            par_read = pp.read_par(f.name, par_length_unit="um")
            os.unlink(f.name)
        self.assertAlmostEqual(par_read["wavelength"], 1.5406e-10, delta=1e-15)

    def test_zero_pixel_size_handled(self):
        """Zero pixel sizes produce well-defined results (distance, tilts unaffected)."""
        par = make_base_par()
        par["y_size"] = 0.0
        par["z_size"] = 0.0
        # Conversion should not crash
        poni = pp.par_to_poni(par)
        # Distance and tilts still valid
        self.assertGreater(abs(poni["dist"]), 0)

    def test_orientation_mapping_completeness(self):
        """All 4 orientations map correctly (corrected from 4x4 affine analysis)."""
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

    def test_too_large_tilts(self):
        """Tilts up to ±π/4 round-trip correctly."""
        par = make_base_par()
        for tilt_key in ["tilt_x", "tilt_y", "tilt_z"]:
            for angle in [0.0, 0.78, -0.78]:  # ±45°
                with self.subTest(tilt=tilt_key, angle=angle):
                    p = dict(par)
                    for k in ["tilt_x", "tilt_y", "tilt_z"]:
                        p[k] = angle if k == tilt_key else 0.0
                    poni = pp.par_to_poni(p)
                    par2 = pp.poni_to_par(poni)
                    self.assertAlmostEqual(angle, par2[tilt_key], delta=1e-8)


if __name__ == "__main__":
    unittest.main(verbosity=2)
