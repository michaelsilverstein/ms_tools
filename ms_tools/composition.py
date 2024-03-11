"Compositional data utilities"

import numpy as np
from skbio.stats.composition import clr
from skbio.stats import subsample_counts
from scipy.spatial.distance import euclidean, jensenshannon

def rho(x, y, transform=False):
    """
    Correlation coefficient for compositional data by Erb and Notredame, 2016
    As described in the supplement of https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
    """
    if transform:
        # Perform CLR transform
        x_clr = clr(x + .5)
        y_clr = clr(y + .5)
    else:
        x_clr, y_clr = x, y
    # Compute rho
    rho = 1 - np.nanvar(x_clr - y_clr) / (np.nanvar(x_clr) + np.nanvar(y_clr))
    return rho

def dirichletCLR(x, n=100):
    """
    Dirichlet estimation of CLR transformed counts
    As described in the supplement of https://www.frontiersin.org/articles/10.3389/fmicb.2017.02224/full
    """
    # Monte Carlo sampling
    ds = np.random.dirichlet(x + .5, n)
    # Compute log
    log_ds = np.log(ds)
    # Compute expected value
    ev = (log_ds - log_ds.mean(1).reshape(-1, 1)).mean(0)
    return ev

def collapse(counts, classification, level):
    """
    Collapse `counts` to `level` as provided in `classification`.
    Requires `counts` and `classification` an indexed by OTU
    """
    # Get list of samples
    samples = counts.columns
    # Collapse
    collapsed = counts.join(classification[level]).groupby(level)[samples].sum()
    return collapsed

def rarefaction_curve(x, step_size=10):
    """Perform a rarefaction curve on `x` with the given step size

    Args:
        x (array): Array of counts
        step_size (int): Size of each step in depth
        
    Returns:
        steps, richness
    """
    x = np.array(x)
    depth = x.sum()
    steps = np.arange(step_size, depth, step_size)
    
    # Ensure up to depth
    if steps[-1] != depth:
        steps = np.append(steps, depth)
    
    # Calculate richness at each step size
    richnesses = np.array([(subsample_counts(x, step) > 0).sum() for step in steps])
    
    return steps, richnesses

def aitchison(x, y, pseudo=1, subset=True):
    """
    Compute the Aitchison distance on vectors x and y after adding a pseudocount to allow for log transformation.
    If subset = True, subset x and y to entries that are present in either or both
    Compatible with scipy.spatial.distance.pdist
    """
    x = np.array(x.copy())
    y = np.array(y.copy())
    
    if subset:
        union = (x > 0) | (y > 0)
        x = x[union]
        y = y[union]
    
    # Add pseudo count
    x += pseudo
    y += pseudo
    
    # Perform clr transform
    x_clr = clr(x)
    y_clr = clr(y)
    return euclidean(x_clr, y_clr)

def overlap_aitchison(x, y):
    """Compute Aitchison distance of overlap between `x` and `y` with metric of choice

    Args:
        x (array): 
        y (array): 
        metric (function, optional): Defaults to Euclidean
    """
    x = np.array(x)
    y = np.array(y)
    
    data = np.array([x, y])
    
    overlap = (data > 0).all(0)
    
    data_overlap_clr = clr(data[:, overlap])
    
    dist = euclidean(*data_overlap_clr)
    
    return dist

def overlap(x, y, relative=False):
    """Overlap distance as described by Bashin 2016 (10.1038/nature18301)
    
    For relative abundance profiles x and y
    
    Overlap(x, y) = sum_i (x_i, y_i) / 2, for shared taxa i.
    
    Overlap is like Jaccard, but robust to "OTU splitting".

    Args:
        {x, y}: Count (default) or relative profiles
        relative: If profiles are counts (default) or relative abundances
    """
    x = np.array(x)
    y = np.array(y)
    
    # If counts are provided, compute relative abundances
    if not relative:
        x = x / x.sum()
        y = y / y.sum()
    
    data = np.array([x, y])
    
    shared_idx = (data > 0).all(0)
    
    overlap_dist = data[:, shared_idx].mean(0).sum()
    
    return overlap_dist

def dissimilarity(x, y):
    """Dissimilarity distance as described by Bashin 2016 (10.1038/nature18301)
    
    Dissimilarity(x, y) = sqrt(JSD(x*, y*)) are count (or relative abundances)
    of shared taxa, renormalized to just those shared taxa

    Args:
        {x, y}: Count (or relative abundance) profiles
    """
    x = np.array(x)
    y = np.array(y)
    
    data = np.array([x, y])
    
    shared_idx = (data > 0).all(0)
    
    shared_data = data[:, shared_idx]
    shared_normalized = shared_data / shared_data.sum(1).reshape(-1, 1)
    
    dissimilarity_dist = np.sqrt(jensenshannon(*shared_normalized))
    
    return dissimilarity_dist