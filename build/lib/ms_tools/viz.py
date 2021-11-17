import matplotlib.pyplot as plt
import seaborn as sns

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