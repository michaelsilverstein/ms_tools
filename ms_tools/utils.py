"Utility functions"

import pandas as pd

def check_len_n(obj, attr, n):
    "Check that `attr` from `obj` has length of n"
    val = getattr(obj, attr)
    if len(val) != n:
        raise ValueError(f'"{attr}" must be of length {n}')
    return True

def check_n_sheets(filepath, n):
    "Check that an Excel workbook contains `n` sheets"
    with pd.ExcelFile(filepath) as fh:
        n_sheets = len(fh.sheet_names)
    if n_sheets != n:
        raise ValueError(f'Excel file "{filepath}" must contain {n} sheets.')
    return True