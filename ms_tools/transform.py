"Data transformation utilities"

import pandas as pd

def pair_md(df: pd.DataFrame, md: pd.DataFrame, sample: str = 'Sample', cols: list = None, wideform: bool = True):
    """
    Given a wideform pairwise distance table, `df` with index and columns containing sample names, add associated
    metadata for each sample from `md` in a longform table. Include all columns of metadata or those listed in `cols`
    If df is already a melted pairwise distance matrix, set `wideform=False`
    :param df: A pairwise distance matrix with index and columns containing sample names that appear in `md`
    :param md: A dataframe with columns || Sample | feature_1 | ... | feature_n ||, where 'Sample' matches those listed
        in `df`
    :param sample:
    :param cols: A list of feature columns to pair with distances ([feature_1, feature_2, ...])
    :param wideform: Assign to False if `df` is already melted
    :return combo: Longform pairwise distance data with metadata with columns:
        || to | from | distance | to_feature_1 | ... | from_feature_1 | ... | same_feature_1 | ... ||
    """
    """Checks"""
    # If a wideform, stack
    if wideform:
        df = df.rename_axis(index='sample1', columns='sample2').stack().reset_index(name='distance')
    # Check if 'to' and 'from' are in df
    if any([c not in df for c in ['sample1', 'sample2']]):
        raise KeyError('`df` must contain the columns "sample1" and "sample2"')
    # Check that sample is in md
    if sample not in md:
        raise KeyError('"%s" must be column in `md`' % sample)
    # Get feature columns
    if cols:
        missing_cols = [c for c in cols if c not in md]
        if missing_cols:
            raise KeyError('The following `cols` are not in `md`: "%s"' % missing_cols)

    # Otherwise, use all columns (either way, remove Sample)
    cols = [c for c in md if c!=sample]
    """Get metadata"""
    # Merge 'to' metadata
    to_md = df.merge(md.rename(columns={sample: 'sample1'}), how='left')[cols].rename(
        columns={c: 'sample1_%s' % c for c in cols})
    # Merge 'from' metadata
    from_md = df.merge(md.rename(columns={sample: 'sample2'}), how='left')[cols].rename(
        columns={c: 'sample2_%s' % c for c in cols})
    """Combine metadata"""
    # Combine
    combo = pd.concat([df, to_md, from_md], axis = 1)
    # Indicate what are the same
    combo['same_' + sample] = combo.sample1.eq(combo.sample2)
    for c in cols:
        same_col = 'same_' + c
        to_col, from_col = [x + c for x in ['sample1_', 'sample2_']]
        combo[same_col] = combo[to_col].eq(combo[from_col])
    return combo

def derep_pairs(data, sample1='sample1', sample2='sample2'):
    """
    Dereplicate paired data where each row describes one pair and both permutations are currently represented
    Input:
        data: Dataframe with columns | `sample1` | `sample2` | ... |
        sample{1, 2}: Name of columns
    Output:
        `data` with 1 row per pair
    """
    # Gather name of all samples
    samples = data[[sample1, sample2]].stack().unique()
    # Create one direction pairs
    pairs = pd.DataFrame([[samples[i], samples[j]] for i in range(len(samples)) for j in range(i, len(samples))], columns=[sample1, sample2])
    # Dereplicate
    derep = data.merge(pairs)
    return derep