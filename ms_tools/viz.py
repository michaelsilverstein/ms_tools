"Data visualization utilities"

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

def zebra(ax=None, orient='v', color='gray', alpha=.3, zorder=0, **kwargs):
    """
    Stripe figure - Color every other x position with `fill_between()`
    If no ax provided, use current ax
    Input:
    | ax: Axes handle
    | orient: 'v' for vertical stripes or 'h' for horizontal stripes
    | Any other argument accepted by `ax.fillbetween{x}()`
    
    Usage:
    1) On most recent plot
    > plt.plot(...)
    > zebra()
    
    2) On specific axis object
    > fig, ax = plt.subplots(...)
    > zebra(ax)
    """
    if not ax:
        ax = plt.gca()
    kwargs.update({'color': color, 'alpha': alpha, 'zorder': zorder})

    # Get lims to reset afterwards
    xlim, ylim = ax.get_xlim(), ax.get_ylim()

    if orient == 'v':
        # Choose x positions to color
        xs = ax.get_xticks()[::2]
        for x in xs:
            ax.fill_between((x - .5, x + .5), ylim[0], ylim[1], **kwargs)

    elif orient == 'h':
        # Choose y positions to color
        ys = ax.get_yticks()[::2]
        for y in ys:
            ax.fill_betweenx((y - .5, y + .5), xlim[0], xlim[1], **kwargs)
    else:
        raise ValueError("orient must be 'v' or 'h'")

    # Reset lims
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    return ax

def stackedbarplot(data, stack_order=None, palette=None, **barplot_kws):
    """
    Create a stacked barplot
    Inputs:
    | data <pd.DataFrame>: A wideform dataframe where the index is the variable to stack, the columns are different samples (x-axis), and the cells the counts (y-axis)
    | stack_order <array-like>: The order for bars to be stacked (Default: given order)
    | palette <array-like>: The colors to use for each value of `stack_order` (Default: husl)
    | barplot_kws: Arguments to pass to sns.barplot()
    
    Author: Michael Silverstein
    Usage: https://github.com/michaelsilverstein/Pandas-and-Plotting/blob/master/lessons/stacked_bar_chart.ipynb
    """
    # Order df
    if stack_order is None:
        stack_order = data.index
    # Create palette if none
    if palette is None:
        palette = dict(zip(stack_order, sns.husl_palette(len(stack_order))))
    # Compute cumsum
    cumsum = data.loc[stack_order].cumsum()
    # Melt for passing to seaborn
    cumsum_stacked = cumsum.stack().reset_index(name='count')
    # Get name of variable to stack and sample
    stack_name, sample_name = cumsum_stacked.columns[:2]
    
    # Plot bar plot
    for s in stack_order[::-1]:
        # Subset to this stack level
        d = cumsum_stacked[cumsum_stacked[stack_name].eq(s)]
        sns.barplot(x=sample_name, y='count', hue=stack_name, palette=palette, data=d, **barplot_kws)
    return plt.gca()

def draw_arrow(start, end, drawer='annotate', ax=None, **kwargs):
    """
    Draw arrow from arrays `start` to `end`
    
    `drawer` indicates which matplotlib function, 'annotate' or 'arrow', to use when drawing the arrows.
    |   'annotate' should maintain arrow properties even on axes with different scales (good for mutliple subplots)
    |   The only advantage that I can find with 'arrow' is that the axis limits will automatically adjust to include the arrow
    """
    
    assert drawer in ['annotate', 'arrow'], '`drawer` must be "annotate" or "arrow"'

    if not ax:
        ax = plt.gca()

    if drawer == 'annotate': 
        ax.annotate('', xytext=start, xy=end, arrowprops={**kwargs})

    if drawer == 'arrow':
        # Enforce np.array to allow for vector subtraction
        start, end = map(np.array, (start, end))
        # Calcualte change in x and y
        dx, dy = end - start

        ax.arrow(start[0], start[1], dx, dy, **kwargs)
    
def draw_arrow_series(coordinates, drawer='annotate', **kwargs):
    "Draw a series of arrows from a 2D array `coordinates`"
    for i in range(coordinates.shape[0] - 1):
        start = coordinates[i]
        end = coordinates[i + 1]
        draw_arrow(start, end, drawer, **kwargs)

def lowess_ci(x, y, n_iters=1000, ci=.95, ax=None, line_kwargs={}, fill_between_kwargs={}):
    """Visualize confidence interval around LOWESS curve

    Args:
        x (array): x data
        y (array): y data
        n_iters (int, optional): Number of bootstrap iterations. Defaults to 1000.
        ci (float, optional): Confidence interval size. Defaults to .95.
        ax: matplotlib axis
        {line, fill_between}_kwargs: Keyword arguments for the LOWESS line and CI interval
    """
    # Enforce numpy arrays
    x, y = map(np.array, (x, y))
    sorted_x = np.sort(x)
    
    # Perform Bootstrap
    smooth_ys = []
    for _ in range(n_iters):
        # Choose indices
        idx = np.random.choice(range(len(x)), len(x))
        
        # Sample data
        bootstrap_x, bootstrap_y = x[idx], y[idx]
        
        # Calculate lowess
        smooth_y = lowess(bootstrap_y, bootstrap_x, xvals=sorted_x)
        smooth_ys.append(smooth_y)
    
    # Compute confidence interval over bootstaps
    ci_low = (1 - ci) / 2
    ci_high = ci + ci_low
    
    lower_bound, upper_bound = np.quantile(smooth_ys, (ci_low, ci_high), axis=0)
    
    # Plot
    if not ax:
        ax = plt.gca()
        
    # CI
    ax.fill_between(sorted_x, lower_bound, upper_bound, **fill_between_kwargs)
    
    # LOWESS line
    smooth_x, smooth_y = lowess(y, x).T
    ax.plot(smooth_x, smooth_y, **line_kwargs)
