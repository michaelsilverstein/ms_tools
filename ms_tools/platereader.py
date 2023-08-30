"""
Plate reader utility functions
"""

import pandas as pd
import numpy as np
from typing import List, Tuple
import warnings

class NegativeBiomassWarning(UserWarning):
    pass

class NegativeMicrorespWarning(UserWarning):
    pass

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

class Plate:
    """
    Excel plate reader data for single plate spectrophotometer reading
    """

    def __init__(self, filepath, sheet_name=0, measurement_name='measurement', columns=None, rows=None, column_type=None, row_type=None, plate_format=96, find_idx=True):
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
            df.at[well_idx] = np.nan

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

class CUEexperiment:
    def __init__(self, pre_od: Plate, post_od: Plate, pre_microresp: Plate, post_microresp: Plate, dilution: int, control_wells: List[Tuple]=None, culture_volume: float=500, deepwell_volume: float=2000, bad_wells_od: List[Tuple]=None, bad_wells_microresp: List[Tuple]=None):
        """
        CUE Experiment containing data pertaining to the specified metric
        | {pre, post}_{od, microresp}: Plate objects for pre and post measurements
        | dilution: Dilution factor for biomass OD measurement
        | culture_volume: The volume (uL) of culture in each MicroResp well. (Default: 500uL)
        | deepwell_volume: The volume (uL) of each well in the deepwell plate
        | control_wells: List of wells that are cell-free controls for computing OD background, ex: [('A', 1), ('B', 2)]. Default is all wells in row H.
        | bad_wells_{od, microresp}: Wells to remove from each experiment in the form for `Plate.removeWells()`

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
        self._assumed_background_od = 0.04
        if control_wells is None:
            warnings.warn(f'No control wells listed. Background OD will be assumed to be {self._assumed_background_od}.')
            control_wells = []
        self._control_wells = control_wells
        self._pre_od = pre_od
        self._post_od = post_od
        self._pre_microresp = pre_microresp
        self._post_microresp = post_microresp
        self._bad_wells_od = bad_wells_od
        self._bad_wells_microresp = bad_wells_microresp
        
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
        warnings.simplefilter('once', NegativeBiomassWarning)
        for r in self._delta_biomassC.index:
            for c in self._delta_biomassC.columns:
                well = (r, c)
                if well not in self._control_wells:
                    value = self._delta_biomassC.at[well]
                    if value < 0:
                        warnings.warn('Negative biomass detected in at least one well. All cases set to null.')
                        # Document negative well
                        self._negative_delta_biomass_wells.append(well)
                        # Set value as null
                        self.delta_biomassC.at[well] = np.nan

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
        warnings.simplefilter('once', NegativeMicrorespWarning)
        for r in negative_microresp.index:
            for c in negative_microresp.columns:
                well = (r, c)
                if well not in self._control_wells:
                    value = negative_microresp.at[well]
                    if value:
                        warnings.warn('Negative MicroResp detected in at least one well. All cases set to null.')
                        # Document negative well
                        self._negative_delta_microresp_wells.append(well)
                        # Set value as null
                        self.respirationC.at[well] = np.nan

    @property
    def _negative_cue_wells(self):
        return self._negative_delta_biomass_wells + self._negative_delta_microresp_wells
    
    def _computeCUE(self):
        """
        Compute CUE according to my time-independent re-arrangement of Smith et al. 2021, Ecology Letters

        CUE = 1 / 1 + Rtot/(C1 - C0)
        """
        self._cue_no_null = 1 / 1 + self._respirationC / self._delta_biomassC
        self.cue = 1 / 1 + self.respirationC / self.delta_biomassC
        
    def __repr__(self) -> str:
        return self.cue.__repr__()

    # def computeDeltaBiomass(self):
    #     "Compute the change in biomass from pre to post"
    #     # Change in OD
    #     self.delta_od = self._post_od_adjusted - self._pre_od_adjusted

    #     # Change in biomass
    #     # https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=3&id=108127
    #     vol = 150e-6
    #     # in ug
    #     self.delta_biomass = .56 * self.delta_od * vol * 1e6

    #     # Change in biomass -> C
    #     # https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100649&ver=8&trm=percent+carbon+e+coli&org=
    #     self.delta_biomass_c = self.delta_biomass * .47
    
    # def computeDeltaCO2(self):
    #     "Compute change in CO2 from pre to post"
        
    #     micro_pre, micro_post = self._pre_microresp.df, self._post_microresp.df
        
    #     # Normalize
    #     self._micro_background = micro_pre.mean().mean()
    #     self._micro_normalized = micro_post / micro_pre * self._micro_background

    #     # % CO2
    #     pct_co2 = -.2265 - 1.606/(1 - 6.771 * self._micro_normalized)
    #     # Reverse column names (because of how microresp is done, columns are flipped)
    #     self.pct_co2 = pct_co2.rename(columns=dict(zip(pct_co2.columns, reversed(pct_co2.columns))))

    #     # Mass C from Microresp manual
        
    #     # Estimate headspace
    #     vol = (1 + .15) - .3
    #     temp = 25
    #     # ug of C
    #     self.mass_co2 = self.pct_co2 * vol * (44/22.4) * (12/44) * (273/(273+temp))
        
    # def computeCUE(self):
    #     "Compute CUE from change in biomass and CO2"
        
    #     if not hasattr(self, 'delta_biomass_c'):
    #         self.computeDeltaBiomass()
        
    #     if not hasattr(self, 'mass_co2'):
    #         self.computeDeltaCO2()
        
    #     delta_biomass_c, mass_co2 = self.delta_biomass_c, self.mass_co2
        
    #     # Remove negative values
    #     for df in delta_biomass_c, mass_co2:
    #         df[df.lt(0)] = np.nan
        
    #     self.cue = delta_biomass_c / (delta_biomass_c + mass_co2)
