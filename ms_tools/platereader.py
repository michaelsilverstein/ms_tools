"""
Plate reader utility functions
"""

import pandas as pd
import numpy as np
from typing import List, Tuple
import warnings

from ms_tools.utils import check_len_n, check_n_sheets

# class NegativeBiomassWarning(UserWarning):
#     pass

# class NegativeMicrorespWarning(UserWarning):
#     pass

# warnings.simplefilter('module', UserWarning)
# warnings.simplefilter('module', NegativeBiomassWarning)
# warnings.simplefilter('module', NegativeMicrorespWarning)
        
_assumed_background_OD = 0.04
_deepwell_volume = 1200
_culture_volume = 500

def od2biomassC(od, volume: float):
    """
    Convert an OD value to biomass of carbon (ug) assuming E.coli
    Inputs:
    | od: An OD reading (value or array)
    | volume: Culture volume (uL)
    Outputs:
    | biomass_c: Carbon biomass (ug)
    """
    # Change in biomass
    # https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=3&id=108127
    # in ug (1e-6 is uL -> L, then 1e6 is g -> ug)
    biomass = .56 * od * volume

    # Change in biomass -> C
    # https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100649&ver=8&trm=percent+carbon+e+coli&org=
    biomass_c = biomass * .47

    return biomass_c

def single2list(obj, n):
    "Return an object repeated in a list"

class Plate:
    """
    Excel plate reader data for single plate spectrophotometer reading
    """

    def __init__(self, filepath, sheet_name=0, measurement_name='measurement', columns=None, rows=None, column_type='column', row_type='row', plate_format=96, find_idx=True):
        """
        filepath: Path to plate reader Excel file
        sheet_name: Name of sheet to reader from `filepath`. Must be specified.
        measurement_name: Name of measurement
        columns: Name of columns
        rows: Name of rows
        {column, row}_type: The type of data contained in the {columns, rows} (ex. Sample or condition)
        plate_format: Format of plate. Current models: 96
        find_idx: Automatically find plate reader coordinates in spreadsheet and load it with loadPlateData(). (default: True)
        """
        assert sheet_name is not None, '"sheet_name" must be specified'
        self.filepath = filepath
        self.df = None
        self._original_df = None
        self.sheet_name = sheet_name
        self._original_sheet = pd.read_excel(self.filepath, self.sheet_name)
        self.measurement_name = measurement_name
        self.removed_wells = []

        # Setup plate format
        self._plate_formats = {96: {'data_columns': list(range(1, 13)), 'data_rows': list('ABCDEFGH')}}
        self._plate_format_data = self._plate_formats[plate_format]
        self._n_data_columns, self._n_data_rows = map(len, (self._plate_format_data[x] for x in ('data_columns', 'data_rows')))

        # Setup column and row names
        for idx_name, idx, n in zip(['columns', 'rows'], [columns, rows], [self._n_data_columns, self._n_data_rows]):
            if idx is not None:
                assert len(idx) == n, f'{n} {idx_name} must be provided.'
        self.columns = columns
        self.rows = rows
        self.column_type = column_type
        self.row_type = row_type

        # Plate reader location info
        self._starting_row_i, self._starting_column_i = None, None
        self._header = None
        self._usecols = None
        self._nrows = None
        self._index_col = None
        

        # Load plate data
        if find_idx:
            self.loadPlateData()

    def _findIdx(self):
        "Find plate reader table indices within Excel sheet"
        df = self._original_sheet

        # Find row where data begins
        data_columns = self._plate_format_data['data_columns']
        starting_row = df.isin(data_columns).sum(1).eq(
            self._n_data_columns).pipe(lambda x: x[x]).index
        starting_row_i = list(df.index).index(starting_row)

        # Find column where data begins
        data_rows = self._plate_format_data['data_rows']
        starting_column = df.isin(data_rows).sum().eq(
            self._n_data_rows).pipe(lambda x: x[x]).index
        starting_column_i = list(df.columns).index(starting_column)

        self._starting_row_i, self._starting_column_i = starting_row_i, starting_column_i
    
        # Assign load parameters
        self._header = self._starting_row_i + 1
        self._usecols = range(self._starting_column_i,
                        self._starting_column_i + self._n_data_columns + 1)
        self._nrows = self._n_data_rows
        self._index_col = 0

    def loadPlateData(self, find_idx=True, read_excel_kwargs={}):
        """
        Load plate reader data from an Excel file.
        | find_idx: Automatically find indices containing plate reader data
        | read_excel_kwargs: Arguments to pass to pd.read_excel()
        """
        if find_idx:
            self._findIdx()
        
        header = self._header
        usecols = self._usecols
        nrows = self._nrows
        index_col = self._index_col

        kwargs = dict(header=header, usecols=usecols, nrows=nrows, index_col=index_col)
        kwargs.update(read_excel_kwargs)
        data = pd.read_excel(self.filepath, self.sheet_name, **kwargs)
        
        # Label rows and columns
        if self.columns is not None:
            data.columns = self.columns
        
        if self.rows is not None:
            data.index = self.rows
        
        # Label rows and column types
        data = data.rename_axis(columns=self.column_type, index=self.row_type)

        # Save
        self._original_df = data
        self.df = data.copy()

    def removeWells(self, wells: List[Tuple]):
        """
        Remove wells - nullify specified cells in DataFrame
        | wells: List of well coordinates as tuples, ex: [('A', 1), ('B', 5)]
        """
        df = self.df
        for well_idx in wells:
            df.at[well_idx] = None

        self.removed_wells.extend(wells)

    def restoreWells(self):
        "Reload original data (following removeWells())"
        self.df = self._original_df.copy()
        self.removed_wells = []
        
    def wellApply(self, wells: List[Tuple], func: callable):
        "Apply a function to a selection of wells"
        return func([self.df.loc[well] for well in wells])

    @property
    def stacked(self):
        if self.df is None:
            self.loadPlateData()
        return self.df.stack().reset_index(name=self.measurement_name)

    def __repr__(self) -> str:
        return self.df.__repr__()
    
    def __eq__(self, other) -> bool:
        # Just check that dataframes are equal
        return (self.df.equals(other.df))

class CUEexperiment:
    def __init__(self, pre_od: Plate, post_od: Plate, pre_microresp: Plate, post_microresp: Plate, dilution: int, control_wells: List[Tuple]=None, culture_volume: float=_culture_volume, deepwell_volume: float=_deepwell_volume, bad_wells_od: List[Tuple]=None, bad_wells_microresp: List[Tuple]=None, name: str=None):
        f"""
        CUE Experiment containing data pertaining to the specified metric
        | {{pre, post}}_{{od, microresp}}: Plate objects for pre and post measurements
        | dilution: Dilution factor for biomass OD measurement
        | culture_volume: The volume (uL) of culture in each MicroResp well. Default: {culture_volume} uL
        | deepwell_volume: The volume (uL) of each well in the deepwell plate. Default: {deepwell_volume} uL
        | control_wells: List of wells that are cell-free controls for computing OD background, ex: [('A', 1), ('B', 2)]. Default: {_assumed_background_OD} as background OD.
        | bad_wells_{{od, microresp}}: Wells to remove from each experiment in the form for `Plate.removeWells()`
        | name: Experiment name
        
        Computes CUE according to Smith et al. 2021, Ecology Letters (doi/10.1111/ele.13840)

        Growth rate:
                                log(C1 / C0)
                            u = ------------
                                    t
        Where C0 and C1 are starting and final biomass, respectively, and t is the duration.

        Respiration rate:
                                    u R_tot
                            R = -----------------
                                C0 (exp(ut) - 1)
        
        Where R_tot is the total respiration measured (MicroResp) over the duration of the experiment.

        And CUE:
                                     u
                            CUE = -------
                                   u + R

        I've simplified this equation to (time-independent, assuming measurements are made over the same time):

                                        1
                            CUE = ------------
                                   1 +  R_tot
                                       -------
                                       C1 - C0
        """
        # Check Plate object
        for plate in pre_od, post_od, pre_microresp, post_microresp:
            if not isinstance(plate, Plate):
                raise TypeError('Must provide "Plate" object.')
        
        # Save
        self._dilution = dilution
        self._culture_volume = culture_volume
        self._deepwell_volume = deepwell_volume
        self._assumed_background_od = _assumed_background_OD
        if control_wells is None:
            # warnings.warn(f'No control wells listed. Background OD will be assumed to be {self._assumed_background_od}.')
            control_wells = []
        self._control_wells = control_wells
        self._pre_od = pre_od
        self._post_od = post_od
        self._pre_microresp = pre_microresp
        self._post_microresp = post_microresp
        self._bad_wells_od = bad_wells_od
        self._bad_wells_microresp = bad_wells_microresp
        self.name = name
        
        self._negative_delta_biomass_wells = []
        self._negative_delta_microresp_wells = []
        
        # Remove bad wells
        self._removeBadWells()

        # Adjust OD
        self._adjustOD()

        # Compute biomass C
        self._computeBiomassC()

        # Compute respiration C
        self._computeRespirationC()
        
        # Compute CUE
        self._computeCUE()
        
        # Stack data
        self._stack_data()

    def _removeBadWells(self):
        # Remove bad wells
        for bad_wells, plates in zip([self._bad_wells_od, self._bad_wells_microresp], [(self._pre_od, self._post_od), (self._pre_microresp, self._post_microresp)]):
            if bad_wells is not None:
                for plate in plates:
                    plate.removeWells(bad_wells)
    
    def _adjustOD(self):
        """Given control wells and dilution, adjust OD"""
        # Compute background OD
        self._pre_od_background = self._assumed_background_od
        self._post_od_background = self._assumed_background_od
        if self._control_wells:
            self._pre_od_background = self._pre_od.wellApply(self._control_wells, np.mean)
            self._post_od_background = self._post_od.wellApply(self._control_wells, np.mean)
        
        # Adjust OD based on empty controls
        od_pre, od_post = self._pre_od.df, self._post_od.df
        self._pre_od_adjusted = (od_pre - self._pre_od_background) * self._dilution
        self._post_od_adjusted = (od_post - self._post_od_background) * self._dilution

    def _computeBiomassC(self):
        "Convert OD to biomass C"
        self._pre_biomassC = od2biomassC(self._pre_od_adjusted, self._culture_volume)
        self._post_biomassC = od2biomassC(self._post_od_adjusted, self._culture_volume)
        
        # Compute delta biomass
        self._delta_biomassC = self._post_biomassC - self._pre_biomassC
        self.delta_biomassC = self._delta_biomassC.copy()
        
        # Check for negative delta biomass
        for r in self._delta_biomassC.index:
            for c in self._delta_biomassC.columns:
                well = (r, c)
                if well not in self._control_wells:
                    value = self._delta_biomassC.at[well]
                    if value < 0:
                        # warnings.warn('Negative biomass detected in at least one well. All cases set to null.')
                        # Document negative well
                        self._negative_delta_biomass_wells.append(well)
                        # Set value as null
                        self.delta_biomassC.at[well] = None

    def _computeRespirationC(self):
        "Convert MicroResp absorbance to respiration C according to manual pg. 15-16"
        # Normalize
        self._microresp_background = self._pre_microresp.df.mean().mean()
        microresp_normalized = self._post_microresp.df / self._pre_microresp.df * self._microresp_background
        
        # Reverse column names (because of how microresp is done, columns are flipped)
        self._microresp_normalized =  microresp_normalized.rename(columns=dict(zip(microresp_normalized.columns, reversed(microresp_normalized.columns))))
        
        # % CO2
        self._pct_co2 = -.2265 - 1.606/(1 - 6.771 * self._microresp_normalized)

        # Compute headspace
        self._microwell_volume = 300
        self._gel_volume = 150
        self._headspace = self._deepwell_volume + self._microwell_volume - self._gel_volume - self._culture_volume

        # Compute CO2-C mass
        self._temp = 25
        self._respirationC = self._pct_co2 / 100 * self._headspace * (44 / 22.4) * (12 / 44) * (273 / (273 + self._temp))
        self.respirationC = self._respirationC.copy()
        
        # Check instances where post > pre (absorbance values decrease with respiration)
        negative_microresp = self._post_microresp.df.gt(self._pre_microresp.df)
        for r in negative_microresp.index:
            for c in negative_microresp.columns:
                well = (r, c)
                if well not in self._control_wells:
                    value = negative_microresp.at[well]
                    if value:
                        # warnings.warn('Negative MicroResp detected in at least one well. All cases set to null.')
                        # Document negative well
                        self._negative_delta_microresp_wells.append(well)
                        # Set value as null
                        self.respirationC.at[well] = None

    @property
    def _negative_cue_wells(self):
        return self._negative_delta_biomass_wells + self._negative_delta_microresp_wells
    
    def _computeCUE(self):
        """
        Compute CUE according to my time-independent re-arrangement of Smith et al. 2021, Ecology Letters

        CUE = 1 / 1 + Rtot/(C1 - C0)
        """
        self._cue_no_null = 1 / (1 + self._respirationC / self._delta_biomassC)
        self.cue = 1 / (1 + self.respirationC / self.delta_biomassC)
        
    def _stack_data(self, attrs=['_delta_biomassC', '_respirationC']):
        "Create a stacked dataframe, `.stacked` with CUE and those listed in `attrs`"
        # Stack CUE
        cue_stacked = self.cue.stack().rename('cue')
        # Stack provided attributes
        attrs_stacked_data = []
        for attr in attrs:
            val = getattr(self, attr)
            if isinstance(val, Plate):
                val = val.df
            val_stacked = val.stack().rename(attr)
            attrs_stacked_data.append(val_stacked)
        attrs_stacked = pd.concat(attrs_stacked_data, axis=1)
        stacked = pd.concat([cue_stacked, attrs_stacked], axis=1)
        stacked = stacked.reindex(stacked.index.sort_values())
        self.stacked = stacked
        
    def __repr__(self) -> str:
        return self.cue.__repr__()

    def __eq__(self, other) -> bool:
        return self.cue.equals(other.cue)
    
class CUEexperiments:
    def __init__(self, od_filepaths: list, microresp_filepaths: list, dilutions, control_wells=None, culture_volumes: float=_culture_volume, deepwell_volumes: float=_deepwell_volume, bad_wells_od: dict={}, bad_wells_microresp: dict={}, names: List[str]=None):
        f"""A collection of `CUEexperiment`s

        Inputs:
        {{od, microresp}}_filepaths (list): A list of Excel filepaths for {{od, microresp}} where each workbook contains 
            a sheet for T0 (sheet 0) and T1 (sheet 1).
        dilutions: Either one number to use the same dilution for all OD measurements or a list for each
        control_wells: Wells for background 
            - None (default): Set background OD to {_assumed_background_OD}
            - List of tuples: Use list of tuples for each well for all experiments, 
            - List of lists of tuples: Use different control wells for each experiment
        culture_volumes (float, optional): Either one number or a list for each experiment. Default: {_culture_volume} uL.
        deepwell_volumes (float, optional): Either one number or a list for each experiment. Default: {_deepwell_volume} uL.
        bad_wells_{{od, microresp}}: Dictionary indicating wells to remove for OD and microresp (ex. {{0: [('A', 2), ('B', 1)]}}). Zero indexed!
        names: A list of names for each experiment or automatically name 0 through N (default).
        """
        ## SETUP
        self._od_filepaths = od_filepaths
        self._microresp_filepaths = microresp_filepaths
        if len(self._od_filepaths) != len(self._microresp_filepaths):
            raise ValueError('The same number of OD and MicroResp filepaths must be provided')
        self.n_experiments = len(self._od_filepaths)

        if isinstance(dilutions, (int, float)):
            dilutions = [dilutions] * self.n_experiments
        self._dilutions = tuple(dilutions)
        
        # Handle None control_wells
        if isinstance(control_wells, type(None)):
            control_wells = [control_wells] * self.n_experiments
        # Check for experiment-specific control wells
        controls_lists_of_lists = all(isinstance(el, list) for el in control_wells)
        if controls_lists_of_lists & (len(control_wells) != self.n_experiments):
            raise ValueError('If different `control_wells` are provided for each experiment there must be only one for each experiment')
        # Check for single list
        controls_single_list = all(isinstance(el, tuple) for el in control_wells)
        if controls_single_list:
            control_wells = [control_wells] * self.n_experiments
        self._control_wells = tuple(control_wells)
        
        if isinstance(culture_volumes, (int, float)):
            culture_volumes = [culture_volumes] * self.n_experiments
        self._culture_volumes = tuple(culture_volumes)
        
        if isinstance(deepwell_volumes, (int, float)):
            deepwell_volumes = [deepwell_volumes] * self.n_experiments
        self._deepwell_volumes = tuple(deepwell_volumes)
        
        if isinstance(names, type(None)):
            names = range(self.n_experiments)
        self.names = tuple(names)
        
        # Check lengths of dilutions, control_wells, culture_volumes, deepwell_volumes, and names
        for attr in ('_dilutions', '_control_wells', '_culture_volumes', '_deepwell_volumes', 'names'):
            check_len_n(self, attr, self.n_experiments)
        
        self._bad_wells_od = tuple(bad_wells_od.get(i) for i in range(self.n_experiments))
        self._bad_wells_microresp = tuple(bad_wells_microresp.get(i) for i in range(self.n_experiments))
        
        
        ## READ IN DATA
        self.pre_ods = None
        self.post_ods = None
        self.pre_microresps = None
        self.post_microresps = None
        self._read_excel_files()
        
        ## GENERATE CUE COLLECTION
        self.cues = None
        self._create_cues()
        
        # Stack data
        self.stacked = None
        self._stack_data()
        
    def _read_excel_files(self):
        for datatype in ['od', 'microresp']:
            filepath_attr = f'_{datatype}_filepaths'
            filepaths = getattr(self, filepath_attr)
            pre, post = [], []
            for filepath in filepaths:
                check_n_sheets(filepath, 2)
                for sheet, l in zip([0, 1], [pre, post]):
                    p = Plate(filepath, sheet, datatype)
                    l.append(p)

            pre = tuple(pre)
            pre_plates_attr = f'pre_{datatype}s'
            setattr(self, pre_plates_attr, pre)
            
            post = tuple(post)
            post_plates_attr = f'post_{datatype}s'
            setattr(self, post_plates_attr, post)
            
    def _create_cues(self):
        cues = []
        for pre_od, post_od, pre_microresp, post_microresp, dilution, control_wells, culture_volume, deepwell_volume, bad_wells_od, bad_wells_microresp, name in\
            zip(self.pre_ods, self.post_ods, self.pre_microresps, self.post_microresps, self._dilutions, self._control_wells, self._culture_volumes, self._deepwell_volumes, self._bad_wells_od, self._bad_wells_microresp, self.names):
            
            kwargs = dict(pre_od=pre_od, 
            post_od=post_od, 
            pre_microresp=pre_microresp, 
            post_microresp=post_microresp, 
            dilution=dilution, 
            control_wells=control_wells, 
            culture_volume=culture_volume, 
            deepwell_volume=deepwell_volume, 
            bad_wells_od=bad_wells_od, 
            bad_wells_microresp=bad_wells_microresp, 
            name=name
            )
            cue = CUEexperiment(**kwargs)
            cues.append(cue)
        self.cues = tuple(cues)
        
        
    def _stack_data(self, attrs=['_delta_biomassC', '_respirationC']):
        "Return a stacked dataframe with `attrs` and cue for each experiment"
        stacked_data = []
        for cue in self.cues:
            # If provided attributes aren't in stacked, add them
            if not pd.Index(attrs).isin(cue.stacked.columns).all():
                cue._stack_data(attrs)
            stacked = cue.stacked.reset_index()
            stacked['experiment'] = cue.name
            stacked_data.append(stacked)
        stacked = pd.concat(stacked_data)
        self.stacked = stacked
