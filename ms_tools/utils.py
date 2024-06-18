"Utility functions"

import pandas as pd
import scipy.stats as sps
import numpy as np

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

def zscore(x, pop_mean, pop_std, alternative='two-sided'):
    """
    Calculate the z-score and p-value for a point x being drawn from a normal distribution
    
    Inputs:
    | x: Sample to test
    | pop_mean: Population mean
    | pop_std: Population standard deviation
    | alternative: Hypothesis test
        * 'two-sided': Two-tailed test
        * 'less': Left-tailed test
        * 'greater': Right-tailed test
        
    Returns:
    | p: p-value given alternative hypothesis
    | z: z-score
    """
    if alternative not in ('two-sided', 'less', 'greater'):
        raise KeyError('"alternative" must be "two-sided", "less", or "greater"')
    
    # Compute cdf
    cdf = sps.norm(pop_mean, pop_std).cdf(x)
    # Compute z-score
    z = sps.norm.ppf(cdf)
    # Compute p-value
    if alternative == 'less':
        p = cdf
    elif alternative == 'two-sided':
        p = 2 * sps.norm.sf(abs(z))
    elif alternative == 'greater':
        p = 1 - cdf
    
    return p, z

def rarefy_present(df, n, axis=0):
    """Subsample `n` existing (True or > 0) entries from a DataFrame or array `df`.
    This will turn entries to "False" - it will not remove them from the object.

    Args:
        df (DataFrame or array): A DataFrame to subsample from
        n (int): Number of entries to subsample
        axis (int): DataFrame axis to subsample 
    """
    # Note of original input type
    initial_type = type(df)
    
    # Convert to DataFrame
    df = pd.DataFrame(df).copy()
    
    # Find present indices
    pres_idx = df.any(axis=axis).pipe(lambda x: x[x]).index
    
    # Number of indices to remove
    n_idx_pres = pres_idx.size
    n_to_remove = n_idx_pres - n
    
    # Select entries to remove
    to_remove = np.random.choice(range(n_idx_pres), n_to_remove, False)
    
    # Remove entries 
    df.loc[pres_idx[to_remove]] = False
    
    # Convert back to original
    df = initial_type(df)
    
    return df

def interpolate_fraction(x, y, t):
    """Interperolate a fraction `t` between two vectors

    Args:
        x (array): A vector of size n
        y (arrray): A vector of size n
        t (float): Fraction between `x` and `y`
            t = 0: x
            t = 1: y
    """
    
    interp = (1 - t) * np.array(x) + t * y
    
    return interp