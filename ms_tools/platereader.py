"""
Plate reader utility functions
"""

import pandas as pd


class Plate:
    """
    Excel plate reader data for single plate spectrophotometer reading
    """

    def __init__(self, filepath, sheet_name=0, measurement_name='measurement', columns=None, rows=None, column_type=None, row_type=None, plate_format=96):
        """
        filepath: Path to plate reader Excel file
        sheet_name: Name of sheet to reader from `filepath`. Must be specified.
        measurement_name: Name of measurement
        columns: Name of columns
        rows: Name of rows
        {column, row}_type: The type of data contained in the {columns, rows} (ex. Sample or condition)
        plate_format: Format of plate. Current models: 96
        """
        assert sheet_name is not None, '"sheet_name" must be specified'
        self.filepath = filepath
        self.df = None
        self.sheet_name = sheet_name
        self._original_df = pd.read_excel(self.filepath, self.sheet_name)
        self.measurement_name = measurement_name

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

        # Find indices
        self.starting_row_i, self.starting_column_i = None, None

        # Load plate data
        self._loadPlateData()

    def _findIdx(self):
        "Find plate reader table indices within Excel sheet"
        df = self._original_df

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

        self.starting_row_i, self.starting_column_i = starting_row_i, starting_column_i

    def _loadPlateData(self, find_idx=True, read_excel_kwargs={}):
        """
        Load plate reader data from an Excel file.
        | find_idx: Automatically find indices containing plate reader data
        | read_excel_kwargs: Arguments to pass to pd.read_excel()
        """
        header = None
        usecols = None
        nrows = None
        index_col = None
        if find_idx:
            self._findIdx()
            header = self.starting_row_i + 1
            usecols = range(self.starting_column_i,
                            self.starting_column_i + self._n_data_columns + 1)
            nrows = self._n_data_rows
            index_col = 0

        kwargs = dict(header=header, usecols=usecols, nrows=nrows, index_col=index_col)
        kwargs.update(read_excel_kwargs)
        data = pd.read_excel(self.filepath, self.sheet_name, **kwargs)

        # Label rows and columns
        if self.columns is not None:
            data.columns = self.columns
            data = data.rename_axis(columns=self.column_type)
        
        if self.rows is not None:
            data.index = self.rows
            data = data.rename_axis(index=self.row_type)

        # Save
        self.df = data

    #TODO: adjust for background and dilution

    @property
    def stacked(self):
        if self.df is None:
            self.loadPlateData()
        return self.df.stack().reset_index(name=self.measurement_name)
