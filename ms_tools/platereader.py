"""
Plate reader utility functions
"""

import pandas as pd
import numpy as np
from typing import List, Tuple


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
            df.loc[well_idx] = np.nan

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

class CUEexperiment:
    def __init__(self, pre_biomass: Plate, post_biomass: Plate, pre_microresp: Plate, post_microresp: Plate, dilution: int, control_wells: List[Tuple]=None, bad_wells_biomass: List[Tuple]=None, bad_wells_microresp: List[Tuple]=None):
        """
        CUE Experiment containing data pertaining to the specified metric
        | {pre, post}_{biomass, microresp}: Plate objects for pre and post measurements
        | dilution: Dilution factor for biomass OD measurement
        | control_wells: List of wells that are cell-free controls, ex: [('A', 1), ('B', 2)]. Default is all wells in row H.
        | bad_wells_{biomass, microresp}: Wells to remove from each experiment in the form for `Plate.removeWells()`
        """
        # Check Plate object
        for plate in pre_biomass, post_biomass, pre_microresp, post_microresp:
            if not isinstance(plate, Plate):
                raise TypeError('Must provide "Plate" object.')

        # Remove bad wells
        for bad_wells, plates in zip([bad_wells_biomass, bad_wells_microresp], [(pre_biomass, post_biomass), (pre_microresp, post_microresp)]):
            if bad_wells is not None:
                for plate in plates:
                    plate.removeWells(bad_wells)
        
        # Save
        self.control_wells = control_wells
        self.pre_biomass = pre_biomass
        self.post_biomass = post_biomass
        self.pre_microresp = pre_microresp
        self.post_microresp = post_microresp
        self._bad_wells_biomass = bad_wells_biomass
        self._bad_wells_microresp = bad_wells_microresp

        self.dilution = dilution
        
        # Adjust OD
        self._adjustOD()
        
        # Compute CUE
        self.computeCUE()
    
    def _adjustOD(self):
        """Given control wells and dilution, adjust OD"""
        # Compute background OD
        self._pre_od_background = 0
        self._post_od_background = 0
        if self.control_wells is not None:
            self._pre_od_background = self.pre_biomass.wellApply(self.control_wells, np.mean)
            self._post_od_background = self.post_biomass.wellApply(self.control_wells, np.mean)
        
        # Adjust OD based on empty controls
        od_pre, od_post = self.pre_biomass.df, self.post_biomass.df
        self._od_pre_adjusted = (od_pre - self._pre_od_background) * self.dilution
        self._od_post_adjusted = (od_post - self._post_od_background) * self.dilution
        
    def computeGrowthRate(self):
        "Compute growth rate (Smith et al. 2021, Ecology Letters)"
        pass

            
    def computeDeltaBiomass(self):
        "Compute the change in biomass from pre to post"
        # Change in OD
        self.delta_od = self._od_post_adjusted - self._od_pre_adjusted

        # Change in biomass
        # https://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&v=3&id=108127
        vol = 150e-6
        # in ug
        self.delta_biomass = .56 * self.delta_od * vol * 1e6

        # Change in biomass -> C
        # https://bionumbers.hms.harvard.edu/bionumber.aspx?id=100649&ver=8&trm=percent+carbon+e+coli&org=
        self.delta_biomass_c = self.delta_biomass * .47
        
    def computeDeltaCO2(self):
        "Compute change in CO2 from pre to post"
        
        micro_pre, micro_post = self.pre_microresp.df, self.post_microresp.df
        
        # Normalize
        self._micro_background = micro_pre.mean().mean()
        self._micro_normalized = micro_post / micro_pre * self._micro_background

        # % CO2
        pct_co2 = -.2265 - 1.606/(1 - 6.771 * self._micro_normalized)
        # Reverse column names (because of how microresp is done, columns are flipped)
        self.pct_co2 = pct_co2.rename(columns=dict(zip(pct_co2.columns, reversed(pct_co2.columns))))

        # Mass C from Microresp manual
        
        # Estimate headspace
        vol = (1 + .15) - .3
        temp = 25
        # ug of C
        self.mass_co2 = self.pct_co2 * vol * (44/22.4) * (12/44) * (273/(273+temp))
        
    def computeCUE(self):
        "Compute CUE from change in biomass and CO2"
        
        if not hasattr(self, 'delta_biomass_c'):
            self.computeDeltaBiomass()
        
        if not hasattr(self, 'mass_co2'):
            self.computeDeltaCO2()
        
        delta_biomass_c, mass_co2 = self.delta_biomass_c, self.mass_co2
        
        # Remove negative values
        for df in delta_biomass_c, mass_co2:
            df[df.lt(0)] = np.nan
        
        self.cue = delta_biomass_c / (delta_biomass_c + mass_co2)
