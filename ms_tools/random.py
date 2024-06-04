"Functions for random sampling"
import numpy as np
import pandas as pd

def structured_binary_matrix(row_group_sizes, col_group_sizes, association_matrix, group_sparsity, noise):
    """Generate a structured binary matrix with specified associations between groups in rows and columns
    with tunable amounts of noise for group membership

    Args:
        {row, col}_group_sizes (array): Size of each group in rows. For example, for three groups of 
            sizes 100, 50, and 25, provide: [100, 50, 25]
        association_matrix (array): A n_rows X n_cols binary array specifying the associations between
            the rows and columns.
        group_sparsity (float): The probability of belonging to an associated group. A `group_sparsity`
        of 0 will result in no grouped entries and a `group_sparsity` of 1 will sample all group entries.
        noise (float): The probability of non-group entries.
        
    Returns:
        x (array): A structured binary matrix
    """
    
    group_sparsity = 1 - group_sparsity
    
    n_groups_col, n_groups_row = map(len, (col_group_sizes, row_group_sizes))
    n_cols, n_rows = map(sum, (col_group_sizes, row_group_sizes))
    association_matrix = np.array(association_matrix)

    if association_matrix.shape != (n_groups_row, n_groups_col):
        raise ValueError(f'Association matrix shape must match x and y group sizes {(n_groups_col, n_groups_row)}')

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

    col_idx = pd.MultiIndex.from_tuples(idxs['col'], names=['col_group', 'col'])
    row_idx = pd.MultiIndex.from_tuples(idxs['row'], names=['row_group', 'row'])

    M = pd.DataFrame(np.zeros((n_rows, n_cols)), row_idx, col_idx)
    
    "Apply grouping structure"
    col_groups = col_idx.get_level_values('col_group').unique()
    row_groups = row_idx.get_level_values('row_group').unique()
    row_groups_to_col_groups = dict(zip(row_groups, [col_groups[a] for a in association_matrix]))
    
    # Apply group sparsity and noise for each row-col group pairing
    for rg, cg in row_groups_to_col_groups.items():
        M.loc[rg, cg] = group_sparsity
        other_cg = list(set(col_groups).difference(set(cg)))
        M.loc[rg, other_cg] = noise
    
    "Sample entries"
    x = np.random.rand(n_rows, n_cols) < M
    
    return x