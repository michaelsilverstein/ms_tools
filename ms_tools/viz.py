import matplotlib.pyplot as plt

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