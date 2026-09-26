"""Input/analysis checks; passing these does not establish a valid engine run."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from scipy.io import netcdf_file

from build_cases import build, mesh_arrays, vegetation_mask
from analyse import compute_metrics, water_balance

ROOT = Path(__file__).resolve().parent


class DeltaProjectTests(unittest.TestCase):
    def setUp(self):
        self.c = json.loads((ROOT / "experiment.json").read_text())

    def test_mesh_area_orientation_and_resolution(self):
        x, y, faces, edges = mesh_arrays(self.c)
        xf, yf = x[faces], y[faces]
        area = 0.5 * np.sum(xf * np.roll(yf, -1, axis=1) - yf * np.roll(xf, -1, axis=1), axis=1)
        self.assertTrue(np.all(area > 0))
        self.assertAlmostEqual(area.sum(), self.c["length_m"] * self.c["width_m"], places=12)
        self.assertTrue(np.allclose(area, self.c["grid_spacing_m"] ** 2))
        self.assertTrue(np.all((edges >= 0) & (edges < len(x))))

    def test_mirrored_masks(self):
        x, y, faces, _ = mesh_arrays(self.c)
        x, y = x[faces].mean(axis=1), y[faces].mean(axis=1)
        left = vegetation_mask(x, y, "half_delta", self.c)
        right = vegetation_mask(x, self.c["width_m"] - y, "half_delta_mirrored", self.c)
        np.testing.assert_array_equal(left, right)
        self.assertGreater(left.sum(), 0)
        self.assertLess(left.sum(), len(left))
        self.assertFalse(vegetation_mask(x, y, "bare", self.c).any())

    def test_shared_initial_state_and_forced_discharge(self):
        from hydrolib.core.dflowfm import FMModel
        from hydrolib.core.dflowfm.mdu.models import Output
        from pydantic import Field
        # HYDROLIB1.0.1 predates these native keys. Extend only their schema,
        # retaining actual files and all other strict/referenced validation.
        class NativeBudgetOutput(Output):
            mbainterval: float = Field(default=0, alias="MbaInterval", ge=0)
            mbawritenetcdf: int = Field(default=0, alias="MbaWriteNetCDF", ge=0, le=1)
            mbawritecsv: int = Field(default=0, alias="MbaWriteCsv", ge=0, le=1)
        class NativeBudgetFMModel(FMModel):
            output: NativeBudgetOutput = Field(default_factory=NativeBudgetOutput)
        c = dict(self.c, discharge_m3_s=0.000071, grid_spacing_m=0.1)
        with tempfile.TemporaryDirectory(prefix="silurian_inputs_") as temp:
            temp = Path(temp)
            config = temp / "experiment.json"
            config.write_text(json.dumps(c))
            cases = temp / "cases"
            build(config, cases)
            first_bed = None
            for name in c["cases"]:
                folder = cases / name
                bed = np.loadtxt(folder / "bedlevel.xyz")
                spatial = (folder / "spatial.ext").read_text()
                self.assertNotIn("METHOD=4", spatial)
                self.assertEqual(spatial.count("METHOD=5"), 5)
                self.assertIn("QUANTITY=massbalanceareaflume", spatial)
                morphology = (folder / "morphology.mor").read_text()
                self.assertIn("NeuBcSand=false", morphology)
                self.assertIn("IBedCond=5", morphology)
                self.assertIn("transport excl pores plastic", (folder / "bedload_boundary.bcm").read_text())
                if first_bed is None:
                    first_bed = bed
                np.testing.assert_array_equal(bed, first_bed)
                with netcdf_file(folder / "flume_net.nc", "r", mmap=False) as ds:
                    self.assertEqual(int(ds.variables["mesh2d_face_nodes"].start_index), 1)
                    self.assertEqual(ds.variables["mesh2d_face_nodes"].shape[1], 4)
                # Read actual boundary data, not just its presence in the configuration.
                lines = (folder / "boundaries.bc").read_text().splitlines()
                values = [float(lines[i+1]) for i, line in enumerate(lines) if line == "Unit=m3/s"]
                self.assertEqual(values, [c["discharge_m3_s"], c["discharge_m3_s"]])
                # Recurse validates the referenced external forcing and boundary files.
                model = NativeBudgetFMModel(folder / "flow.mdu", recurse=True)
                self.assertEqual(model.geometry.bedlevtype, 1)
                self.assertEqual(model.geometry.kmx, 0)
                self.assertEqual(model.sediment.sedimentmodelnr, 4)
                self.assertEqual(model.output.mbainterval, c["output_interval_s"])
                self.assertEqual(model.output.mbawritenetcdf, 1)

    def test_signed_bed_volume_metrics(self):
        data = {"x": np.array([2.5, 3.0]), "y": np.array([.3, .7]),
                "area": np.array([2., 3.]), "initial": np.zeros(2),
                "final": np.array([.1, -.2]), "times_s": np.array([0., 30.]),
                "bed_variable": "analytic_test"}
        r = compute_metrics(data, {"basin_start_m": 2., "bed_porosity": .4}, .01)
        self.assertAlmostEqual(r["positive_bed_change_bulk_volume_m3"], .2)
        self.assertAlmostEqual(r["negative_bed_change_bulk_volume_m3"], .6)
        self.assertAlmostEqual(r["net_bed_solid_volume_change_m3"], -.24)
        self.assertEqual(r["downstream_deposition_footprint_m2"], 2)
        self.assertEqual(r["downstream_deposition_centroid_x_m"], 2.5)
        data["final"] = np.array([-.1, -.2])
        r = compute_metrics(data, {"basin_start_m": 2., "bed_porosity": .4}, .01)
        self.assertIsNone(r["downstream_deposition_centroid_x_m"])
        self.assertEqual(r["maximum_aggradation_m"], 0)

    def test_native_water_ledger_diagnostic(self):
        from netCDF4 import Dataset
        with tempfile.TemporaryDirectory(prefix="silurian_ledger_") as temp:
            path = Path(temp) / "synthetic_his.nc"
            with Dataset(path, "w") as ds:
                ds.createDimension("time", 3)
                volume = ds.createVariable("water_balance_total_volume", "f8", ("time",))
                error = ds.createVariable("water_balance_volume_error", "f8", ("time",))
                volume.units = error.units = "m3"
                volume[:] = [1, 2, 4]
                error[:] = [0, -2e-6, 1e-6]
            result = water_balance(path)
            self.assertAlmostEqual(result["max_absolute_volume_error_m3"], 2e-6)
            self.assertAlmostEqual(result["max_error_relative_to_volume"], 5e-7)
            self.assertTrue(result["below_diagnostic_threshold"])
            with Dataset(path, "a") as ds:
                ds["water_balance_volume_error"][:] = [0, -.1, .01]
            self.assertFalse(water_balance(path)["below_diagnostic_threshold"])
            with Dataset(path, "a") as ds:
                ds["water_balance_total_volume"][:] = 0
            self.assertEqual(water_balance(path)["status"], "unavailable")
            with Dataset(path, "a") as ds:
                ds["water_balance_total_volume"].units = "m3/s"
            with self.assertRaisesRegex(ValueError, "cubic metres"):
                water_balance(path)


if __name__ == "__main__":
    unittest.main(verbosity=2)
