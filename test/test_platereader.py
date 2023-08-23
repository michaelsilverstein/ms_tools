from unittest import TestCase
import numpy as np

from ms_tools.platereader import Plate, CUEexperiment

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
        self.Plate = Plate(self.test_file_path, self.test_sheet_name, self.measurement_name)
        
    def test_plate_data(self):
        self.assertTrue(np.array_equal(self.Plate.df.values, self.test_data))
        
    def test_plate_df_idx(self):
        pass
        # self.assertTrue(self.Plate)
        
    def test_sheet_name(self):
        with self.assertRaises(AssertionError) as context:
            Plate('', None)
        self.assertEquals('"sheet_name" must be specified', str(context.exception))
        
    def test_measurement_name(self):
        self.assertTrue(self.measurement_name in self.Plate.stacked.columns)
        
    def test_removeAndRestoreWells(self):
        wells_to_remove = [('A', 1), ('B', 5)]
        self.Plate.removeWells(wells_to_remove)
        
        for well_idx in wells_to_remove:
            well_value = self.Plate.df.loc[well_idx]
            self.assertTrue(np.isnan(well_value))
        
        self.Plate.restoreWells()
        
        for well_idx in wells_to_remove:
            well_value = self.Plate.df.loc[well_idx]
            self.assertTrue(np.isreal(well_value))
            
class testCUEexperiment(TestCase):
    def setUp(self):
        microresp_path = 'test/test_data/DE.Microresp1.xlsx'
        biomass_path = 'test/test_data/DE.Transfer1.xlsx'
        self.pre_biomass_plate, self.post_biomass_plate = [Plate(biomass_path, sheet) for sheet in [0, 1]]
        self.pre_microresp_plate, self.post_microresp_plate = [Plate(microresp_path, sheet) for sheet in [0, 1]]
        self.cue = CUEexperiment(self.pre_biomass_plate, self.post_biomass_plate, self.pre_microresp_plate, self.post_microresp_plate, 5)
        
    def test_backgroundOD(self):
        self.assertEqual(0, self.cue._pre_od_background)
        self.assertEqual(0, self.cue._post_od_background)
        
        control_wells = [('H', i) for i in range(1, 13)]
        
        cue = CUEexperiment(self.pre_biomass_plate, self.post_biomass_plate, self.pre_microresp_plate, self.post_microresp_plate, 5, control_wells)
        self.assertAlmostEquals(0.09358333333333334, cue._pre_od_background)
        self.assertAlmostEquals(0.09641666666666666, cue._post_od_background)
        