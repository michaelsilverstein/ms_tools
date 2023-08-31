from unittest import TestCase
import numpy as np

from ms_tools.platereader import Plate, CUEexperiment, CUEexperiments, od2biomassC, culture_volume, deepwell_volume

class testPlate(TestCase):
    def setUp(self):
        self.test_file_path = 'test/test_data/DE.Microresp1.xlsx'
        self.test_sheet_name = 0
        self.measurement_name = 'measure1'
        self.plate = Plate(self.test_file_path, self.test_sheet_name, self.measurement_name)
        
    def test_plate_data(self):
        test_data_text = """0.884	0.879	0.969	0.948	0.901	0.93	0.96	0.8	0.892	0.856	0.905	1.039
0.831	0.821	0.847	0.863	0.894	0.855	0.913	0.835	0.826	0.795	0.822	0.85
0.958	0.963	0.919	0.983	0.909	0.919	0.897	0.967	0.935	0.991	1.029	0.926
0.897	0.863	0.859	0.829	0.915	0.843	0.847	0.814	0.899	0.873	0.893	0.95
1.005	0.955	0.95	0.87	1.023	0.87	0.891	0.931	0.93	0.895	0.949	0.89
0.918	0.974	0.94	0.876	0.959	0.831	0.905	0.975	0.959	0.879	0.973	0.996
0.929	0.947	0.92	0.864	0.882	0.894	0.918	0.97	0.975	0.932	0.926	0.951
0.928	0.986	0.917	0.911	0.963	0.903	0.922	0.992	1.094	1.021	1.268	1.473"""
        test_data = np.array([list(map(float, line.split('\t'))) for line in test_data_text.split('\n')])
        self.assertTrue(np.array_equal(self.plate.df.values, test_data))
        
    def test_sheet_name(self):
        with self.assertRaises(AssertionError) as context:
            Plate('', None)
        self.assertEqual('"sheet_name" must be specified', str(context.exception))
        
    def test_measurement_name(self):
        self.assertTrue(self.measurement_name in self.plate.stacked.columns)
        
    def test_removeAndRestoreWells(self):
        wells_to_remove = [('A', 1), ('B', 5)]
        self.plate.removeWells(wells_to_remove)
        
        for well_idx in wells_to_remove:
            well_value = self.plate.df.at[well_idx]
            self.assertTrue(np.isnan(well_value))
            self.assertIn(well_idx, self.plate.removed_wells) 

        self.plate.restoreWells()
        
        for well_idx in wells_to_remove:
            well_value = self.plate.df.loc[well_idx]
            self.assertTrue(np.isreal(well_value))
            self.assertFalse(self.plate.removed_wells)
            
    def test_wellApply(self):
        wells = [('H', i) for i in range(1, 13)]
        last_row_mean = self.plate.df.loc['H'].mean()
        self.assertEqual(self.plate.wellApply(wells, np.mean), last_row_mean)

    def test_o2biomassC(self):
        # This should be able to accept a value and table
        volume = 500
        od2biomassC(1, volume)
        od2biomassC(self.plate.df, volume)

class testCUEexperiment(TestCase):
    def setUp(self):
        microresp_path = 'test/test_data/DE.Microresp1.xlsx'
        od_path = 'test/test_data/DE.Transfer1.xlsx'
        self.pre_od_plate, self.post_od_plate = [Plate(od_path, sheet) for sheet in [0, 1]]
        self.pre_microresp_plate, self.post_microresp_plate = [Plate(microresp_path, sheet) for sheet in [0, 1]]
        
    def test_backgroundOD(self):
        cue = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5)
        self.assertEqual(cue._assumed_background_od, cue._pre_od_background)
        self.assertEqual(cue._assumed_background_od, cue._post_od_background)
        
        control_wells = [('H', i) for i in range(1, 13)]
        
        cue = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5, control_wells=control_wells)
        self.assertAlmostEqual(0.09358333333333334, cue._pre_od_background)
        self.assertAlmostEqual(0.09641666666666666, cue._post_od_background)
    
    def test_removeBadWells(self):
        od_bad_wells = [('A', 1), ('A', 2)]
        microresp_bad_wells = [('B', 1), ('B', 2)]

        cue = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5, bad_wells_od=od_bad_wells, bad_wells_microresp=microresp_bad_wells)

        self.assertCountEqual(cue._pre_od.removed_wells, od_bad_wells)
        self.assertCountEqual(cue._post_od.removed_wells, od_bad_wells)

        self.assertCountEqual(cue._pre_microresp.removed_wells, microresp_bad_wells)
        self.assertCountEqual(cue._post_microresp.removed_wells, microresp_bad_wells)
        
    def test_negative_delta_biomass(self):
        # Without control wells
        cue_without_control_wells = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5)
        expected_negative_biomass_wells_without_control = [('A', 2), ('C', 10), ('D', 8), ('D', 10), ('E', 2), ('E', 3), ('F', 7), ('H', 5)]
        self.assertCountEqual(expected_negative_biomass_wells_without_control, cue_without_control_wells._negative_delta_biomass_wells)
        
        
        # With control wells
        control_wells = [('H', i) for i in range(1, 13)]
        cue_with_control_wells = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5, control_wells)
        expected_negative_biomass_wells_with_control = [('A', 2), ('C', 2), ('C', 10), ('D', 8), ('D', 10), ('E', 2), ('E', 3), ('F', 7), ('F', 11), ('G', 2)]
        self.assertCountEqual(expected_negative_biomass_wells_with_control, cue_with_control_wells._negative_delta_biomass_wells)
        
        self.assertTrue(all([np.isnan(cue_with_control_wells.delta_biomassC.at[i]) for i in cue_with_control_wells._negative_delta_biomass_wells]))

        # Underscored delta should keep actual delta (no nulls)
        self.assertFalse(cue_with_control_wells._delta_biomassC.isnull().any().any())
        
    def test_negative_microresp(self):
        # All CUE values as expected with original data
        cue_unedited = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5)
        self.assertEqual([], cue_unedited._negative_delta_microresp_wells)

        # Put unexpected relationship
        post_microresp_edited = self.post_microresp_plate
        well_to_edit = ('A', 1)
        post_microresp_edited.df.at[well_to_edit] = self.pre_microresp_plate.df.at[well_to_edit] + .1
        cue_edited = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5)
        
        self.assertCountEqual([well_to_edit], cue_edited._negative_delta_microresp_wells)
        
        # No nulls in underscore
        self.assertFalse(cue_edited._respirationC.isnull().any().any())
        self.assertTrue(np.isnan(cue_edited.respirationC.at[well_to_edit]))
        
    def test_null_values_from_negatives(self):
        post_microresp_edited = self.post_microresp_plate
        well_to_edit = ('A', 1)
        post_microresp_edited.df.at[well_to_edit] = self.pre_microresp_plate.df.at[well_to_edit] + .1
        
        control_wells = [('H', i) for i in range(1, 13)]
        cue_edited = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5, control_wells)
        
        expected_null_wells = [('A', 1), ('A', 2), ('C', 2), ('C', 10), ('D', 8), ('D', 10), ('E', 2), ('E', 3), ('F', 7), ('F', 11), ('G', 2)]
        
        self.assertCountEqual(expected_null_wells, cue_edited._negative_cue_wells)
        
        for r in cue_edited.cue.index:
            for c in cue_edited.cue.columns:
                well = (r, c)
                value = cue_edited.cue.at[well]
                if well in expected_null_wells:
                    self.assertTrue(np.isnan(value))
                else:
                    self.assertTrue(np.isfinite(value))
                    
class testCUEexperiments(TestCase):
    def setUp(self) -> None:
        self.od_filepaths = ['test/test_data/DE.Transfer1.xlsx', 'test/test_data/DE.Transfer3.xlsx']
        self.microresp_filepaths = ['test/test_data/DE.Microresp1.xlsx', 'test/test_data/DE.Microresp3.xlsx']
        self.n_experiments = len(self.od_filepaths)
        self.dilutions = 5
        self.control_wells = [('H', i) for i in range(1, 13)]
        self.cues = CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, self.control_wells)
        
        
    def test_same_number_filepaths(self):
        # Unequal number of paths
        od_filepaths = self.od_filepaths + ['test/test_data/DE.Transfer1.xlsx']
        with self.assertRaises(ValueError) as context:
            CUEexperiments(od_filepaths, self.microresp_filepaths, self.dilutions)
        self.assertEqual('The same number of OD and MicroResp filepaths must be provided', str(context.exception))
        
        # Equal number of paths
        self.assertEqual(2, self.cues.n_experiments)

    def test_dilutions(self):
        # Single value becoming n
        expected_dilutions = tuple([self.dilutions] * self.n_experiments)
        self.assertEqual(expected_dilutions, self.cues._dilutions)
        
        # Explicitly providing n
        dilutions = (5, 5)
        explicit_cues = CUEexperiments(self.od_filepaths, self.microresp_filepaths, dilutions, self.control_wells)
        self.assertTrue(dilutions, explicit_cues._dilutions)
        
        # Raise error
        dilutions = [5, 5, 5]
        with self.assertRaises(ValueError) as context:
            CUEexperiments(self.od_filepaths, self.microresp_filepaths, dilutions, self.control_wells)
        self.assertEqual(f'"_dilutions" must be of length {self.n_experiments}', str(context.exception))
        
    def test_control_wells(self):
        # None case
        none_control_wells = None
        none_wells_cue = CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, none_control_wells)
        expected_control_wells = tuple([None] * self.n_experiments)
        self.assertEqual(expected_control_wells, none_wells_cue._control_wells)
        
        # One list
        expected_control_wells = tuple([self.control_wells] * self.n_experiments)
        self.assertEqual(expected_control_wells, self.cues._control_wells)
        
        # Multiple lists correct
        correct_multiple_control_wells = (self.control_wells, [('A', 1)])
        multiple_control_wells_cues = CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, correct_multiple_control_wells)
        self.assertEqual(correct_multiple_control_wells, multiple_control_wells_cues._control_wells)
        
        # Wrong number of lists
        incorrect_multiple_control_wells = (self.control_wells, [('A', 1)], [('B', 1)])
        with self.assertRaises(ValueError) as context:
            CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, incorrect_multiple_control_wells)
        self.assertEqual('If different `control_wells` are provided for each experiment there must be only one for each experiment', str(context.exception))
        
    def test_culture_volumes(self):
        # Default: No volume provided
        expected_volumes = tuple([culture_volume] * self.n_experiments)
        self.assertEqual(expected_volumes, self.cues._culture_volumes)
        
        # Providing one volume
        provided_vol = 400
        provided_cues = CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, self.control_wells, provided_vol)
        expected_vols = tuple([provided_vol] * self.n_experiments)
        self.assertEqual(expected_vols, provided_cues._culture_volumes)
        
        # Providing multiple expected_volumes
        provided_vols = (100, 200)
        provided_cues = CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, self.control_wells, provided_vols)
        self.assertEqual(provided_vols, provided_cues._culture_volumes)
        
        # Wrong number of volumes
        provided_vols = [100, 200, 300]
        with self.assertRaises(ValueError) as context:
            CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, self.control_wells, provided_vols)
        self.assertEqual(f'"_culture_volumes" must be of length {self.n_experiments}', str(context.exception))

    def test_deepwell_volumes(self):
        # Default: No volume provided
        expected_volumes = tuple([deepwell_volume] * self.n_experiments)
        self.assertEqual(expected_volumes, self.cues._deepwell_volumes)
        
        # Providing one volume
        provided_vol = 400
        provided_cues = CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, self.control_wells, deepwell_volumes=provided_vol)
        expected_vols = tuple([provided_vol] * self.n_experiments)
        self.assertEqual(expected_vols, provided_cues._deepwell_volumes)
        
        # Providing multiple expected_volumes
        provided_vols = (100, 200)
        provided_cues = CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, self.control_wells, deepwell_volumes=provided_vols)
        self.assertEqual(provided_vols, provided_cues._deepwell_volumes)
        
        # Wrong number of volumes
        provided_vols = [100, 200, 300]
        with self.assertRaises(ValueError) as context:
            CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, self.control_wells, deepwell_volumes=provided_vols)
        self.assertEqual(f'"_deepwell_volumes" must be of length {self.n_experiments}', str(context.exception))
    
    def test_bad_wells(self):
        bad_od_wells = {1: [('A', 1), ('B', 2)]}
        bad_microresp_wells = {0: [('C', 3)]}
        
        cues = CUEexperiments(self.od_filepaths, self.microresp_filepaths, self.dilutions, self.control_wells, bad_wells_od=bad_od_wells, bad_wells_microresp=bad_microresp_wells)
        
        expected_bad_od_wells = (None, [('A', 1), ('B', 2)])
        self.assertEqual(expected_bad_od_wells, cues._bad_wells_od)
        
        expected_bad_microresp_wells = ([('C', 3)], None)
        self.assertEqual(expected_bad_microresp_wells, cues._bad_wells_microresp)
        
    def test_read_excel_files(self):
        pre_ods = [Plate(filepath, 0, 'od') for filepath in self.od_filepaths]
        post_ods = [Plate(filepath, 1, 'od') for filepath in self.od_filepaths]
        
        pre_microresps = [Plate(filepath, 0, 'od') for filepath in self.microresp_filepaths]
        post_microresps = [Plate(filepath, 1, 'od') for filepath in self.microresp_filepaths]
        
        for experiment in range(self.n_experiments):
            pre_od_df = pre_ods[experiment].df
            post_od_df = post_ods[experiment].df
            
            pre_microresp_df = pre_microresps[experiment].df
            post_microresp_df = post_microresps[experiment].df
            
            # Only test equality of DataFrames bc doing it with whole `Plate` object is obnoxious
            pre_od_eval = pre_od_df.eq(self.cues.pre_ods[experiment].df).all().all()
            self.assertTrue(pre_od_eval)
            
            post_od_eval = post_od_df.eq(self.cues.post_ods[experiment].df).all().all()
            self.assertTrue(post_od_eval)
            
            pre_microresp_eval = pre_microresp_df.eq(self.cues.pre_microresps[experiment].df).all().all()
            self.assertTrue(pre_microresp_eval)
            
            post_microresp_eval = post_microresp_df.eq(self.cues.post_microresps[experiment].df).all().all()
            self.assertTrue(post_microresp_eval)
            
    
    