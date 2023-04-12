"Compositional data utilities"

import numpy as np
from skbio.stats.composition import clr
from skbio.stats import subsample_counts

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