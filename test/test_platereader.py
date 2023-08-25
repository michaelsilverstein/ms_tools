from unittest import TestCase
import numpy as np

from ms_tools.platereader import Plate, CUEexperiment, od2biomassC

class testPlate(TestCase):
    def setUp(self):
        self.test_file_path = 'test/test_data/DE.Microresp1.xlsx'
        self.test_sheet_name = 0
        self.measurement_name = 'measure1'
        test_data_text = """0.884	0.879	0.969	0.948	0.901	0.93	0.96	0.8	0.892	0.856	0.905	1.039
0.831	0.821	0.847	0.863	0.894	0.855	0.913	0.835	0.826	0.795	0.822	0.85
0.958	0.963	0.919	0.983	0.909	0.919	0.897	0.967	0.935	0.991	1.029	0.926
0.897	0.863	0.859	0.829	0.915	0.843	0.847	0.814	0.899	0.873	0.893	0.95
1.005	0.955	0.95	0.87	1.023	0.87	0.891	0.931	0.93	0.895	0.949	0.89
0.918	0.974	0.94	0.876	0.959	0.831	0.905	0.975	0.959	0.879	0.973	0.996
0.929	0.947	0.92	0.864	0.882	0.894	0.918	0.97	0.975	0.932	0.926	0.951
0.928	0.986	0.917	0.911	0.963	0.903	0.922	0.992	1.094	1.021	1.268	1.473"""
        self.test_data = np.array([list(map(float, line.split('\t'))) for line in test_data_text.split('\n')])
        self.plate = Plate(self.test_file_path, self.test_sheet_name, self.measurement_name)
        
    def test_plate_data(self):
        self.assertTrue(np.array_equal(self.plate.df.values, self.test_data))
        
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
        self.cue = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5)
        
    def test_backgroundOD(self):
        self.assertEqual(self.cue._assumed_background_od, self.cue._pre_od_background)
        self.assertEqual(self.cue._assumed_background_od, self.cue._post_od_background)
        
        control_wells = [('H', i) for i in range(1, 13)]
        
        cue = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5, control_wells=control_wells)
        self.assertAlmostEquals(0.09358333333333334, cue._pre_od_background)
        self.assertAlmostEquals(0.09641666666666666, cue._post_od_background)
    
    def test_removeBadWells(self):
        od_bad_wells = [('A', 1), ('A', 2)]
        microresp_bad_wells = [('B', 1), ('B', 2)]

        cue = CUEexperiment(self.pre_od_plate, self.post_od_plate, self.pre_microresp_plate, self.post_microresp_plate, 5, bad_wells_od=od_bad_wells, bad_wells_microresp=microresp_bad_wells)

        self.assertCountEqual(cue._pre_od.removed_wells, od_bad_wells)
        self.assertCountEqual(cue._post_od.removed_wells, od_bad_wells)

        self.assertCountEqual(cue._pre_microresp.removed_wells, microresp_bad_wells)
        self.assertCountEqual(cue._post_microresp.removed_wells, microresp_bad_wells)
        
