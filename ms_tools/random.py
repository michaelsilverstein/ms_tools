"Functions for random sampling"
import numpy as np
import pandas as pd

def structured_binary_matrix(row_group_sizes, col_group_sizes, association_matrix, group_probability=.5, nongroup_probability=0, row_group_name=None, col_group_name=None, row_name=None, col_name=None):
    """Generate a structured binary matrix with specified associations between groups in rows and columns
    with tunable amounts of noise for group membership

    Args:
        {row, col}_group_sizes (array): Size of each group in rows. For example, for three groups of 
            sizes 100, 50, and 25, provide: [100, 50, 25]
        association_matrix (array): A n_rows X n_cols binary array specifying the associations between
            the rows and columns.
        group_probability (float): The probability of belonging to an associated group. A `group_sparsity`
        of 0 will result in no grouped entries and a `group_sparsity` of 1 will sample all group entries.
            Default = 0.5
        nongroup_probability (float): The probability of non-group entries.
            Default = 0.
        {row, col}{_group}_name (str): Name of group and individual row and column entries.
        
    Returns:
        x (array): A structured binary matrix
    """
    
    "Inputs"
    n_groups_col, n_groups_row = map(len, (col_group_sizes, row_group_sizes))
    n_cols, n_rows = map(sum, (col_group_sizes, row_group_sizes))
    association_matrix = np.array(association_matrix)

    if association_matrix.shape != (n_groups_row, n_groups_col):
        raise ValueError(f'Association matrix shape must match x and y group sizes {(n_groups_col, n_groups_row)}')
    
    if not row_group_name:
        row_group_name = 'row_group'
    if not col_group_name:
        col_group_name = 'col_group'
    if not row_name:
        row_name = 'row'
    if not col_name:
        col_name = 'col'

    "Initialize DataFrame"
    # Make an empty DataFrame with indices that allow for accessing groups
    idxs = {'col': [], 'row': []}
    for name, group_sizes in zip(['col', 'row'], [col_group_sizes, row_group_sizes]):
        x = 0
        abr = name[0].upper()
        for i,gs in enumerate(group_sizes):
            for _ in range(gs):
                idx = (f'{abr}G_{i}', f'{abr}_{x}')
                idxs[name].append(idx)
                x += 1

    col_idx = pd.MultiIndex.from_tuples(idxs['col'], names=[col_group_name, col_name])
    row_idx = pd.MultiIndex.from_tuples(idxs['row'], names=[row_group_name, row_name])

    M = pd.DataFrame(np.zeros((n_rows, n_cols)), row_idx, col_idx)
    
    "Apply grouping structure"
    col_groups = col_idx.get_level_values(col_group_name).unique()
    row_groups = row_idx.get_level_values(row_group_name).unique()
    row_groups_to_col_groups = dict(zip(row_groups, [col_groups[a] for a in association_matrix]))
    
    # Apply group and nongroup probability for each row-col group pairing
    for rg, cg in row_groups_to_col_groups.items():
        M.loc[rg, cg] = group_probability
        other_cg = list(set(col_groups).difference(set(cg)))
        M.loc[rg, other_cg] = nongroup_probability
    
    "Sample entries"
    x = np.random.rand(n_rows, n_cols) < M
    
    return x

def shuffle_rows(x, frac=None, n=None):
    """Shuffle a specified fraction or number of rows in `x` and return a new array. Unshuffled rows of x
    will remain in place.

    Args:
        x (array): Array to shuffle. Only the first axis of `x` will be shuffled.
        frac (float): Fraction of entries to shuffle. * 
        n (int): Number of entries to shuffle.  *
        
    * Only `frac` or `n` can be provided.
    """
    
    frac_provided = frac is not None
    n_provided = n is not None
    
    if frac_provided and n_provided:
        raise ValueError('Only `frac` or `n` can be supplied.')
    
    x = x.copy()
    
    # Number of rows in x
    n_entries = x.shape[0]
    
    # Number of entries to shuffle
    if frac_provided:
        n_to_shuffle = np.ceil(frac * n_entries).astype(int)
    else:
        n_to_shuffle = n
    
    # Entries to shuffle
    entries_to_shuffle = np.random.choice(n_entries, n_to_shuffle, replace=False)
    
    # Swap positions of selected entries
    new_positions = np.random.permutation(entries_to_shuffle)
    x[new_positions] = x[entries_to_shuffle]
    
    return x