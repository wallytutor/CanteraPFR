# -*- coding: utf-8 -*-


def plotPFR(data, species, saveas, **kwargs):
    """ Provide formated plot of results.

    TODO
    ----
    Document keyword arguments.

    Parameters
    ----------
    data : pandas.DataFrame
        Table of integration solution.
    species : dict
        Dictionary of species to plot. Keys must match columns of `data`
        and values are used as labels (for allowing LaTeX names).
    saveas : path-like
        Path to graphics file.
    dpi : int, optional
        Output resolution (parsed to pyplot.savefig). Default is 300.
    """

    from matplotlib import pyplot

    cdef dict axes_default = {
        'x': 'Position ($\\mathrm{m}$)',
        'T': 'Temperature ($\\mathrm{K}$)',
        'rho': 'Density ($\\mathrm{kg\\,m^{3}}$)',
        'u': 'Velocity ($\\mathrm{m\\,s^{-1}}$)',
        'p': 'Pressure ($\\mathrm{Pa}$)',
        'Y': 'Mass fraction'
        }

    cdef int dpi = kwargs.get('dpi', 300)
    cdef str style = kwargs.get('style', 'bmh')
    cdef tuple figsize = kwargs.get('figsize', (10, 6))
    cdef dict axes = kwargs.get('axes', axes_default)

    cdef str xaxis = 'x'
    cdef str uaxis = 'u'
    cdef str paxis = 'p'
    cdef str Yaxis = 'Y'
    cdef str Taxis = 'T' if 'T' in data.columns else 'rho'

    pyplot.close('all')
    pyplot.style.use(style)
    pyplot.figure(figsize=figsize)

    pyplot.subplot(221)
    pyplot.plot(data[xaxis], data[Taxis])
    pyplot.xlabel(axes[xaxis])
    pyplot.ylabel(axes[Taxis])

    pyplot.subplot(222)
    pyplot.plot(data[xaxis], data[uaxis])
    pyplot.xlabel(axes[xaxis])
    pyplot.ylabel(axes[uaxis])

    pyplot.subplot(223)
    pyplot.plot(data[xaxis], data[paxis])
    pyplot.xlabel(axes[xaxis])
    pyplot.ylabel(axes[paxis])

    pyplot.subplot(224)
    for key, spec in species.items():
        pyplot.plot(data[xaxis], data[key], label=spec)
    pyplot.xlabel(axes[xaxis])
    pyplot.ylabel(axes[Yaxis])
    pyplot.legend()

    pyplot.tight_layout()
    pyplot.savefig(saveas, dpi=dpi)


def plot_adjmatrix(adjmat, saveas, overwrite=True, **kwargs):
    """ Plot adjacency matrix of graph.

    TODO
    ----
    Document keyword arguments.

    Parameters
    ----------
    adjmat : array-like
        Source graph's adjacency matrix.
    saveas : str
        Path to save file. Returns `None` if file exists.
    overwrite : bool, optional
        If `True` allows plot to be overwritten. Default is `False`.
    """

    from matplotlib import pyplot

    if not overwrite:
        from os.path import exists
        assert not exists(saveas), f' Graph cannot overwrite: {saveas}'

    cdef tuple figsize = kwargs.get('figsize', (6, 6))
    cdef str cmap = kwargs.get('cmap', 'Greys')
    cdef int dpi = kwargs.get('dpi', 300)

    pyplot.close('all')
    fig, ax = pyplot.subplots(figsize=figsize)
    ax.imshow(adjmat, cmap=cmap)
    ax.axes.get_xaxis().set_visible(False)
    ax.axes.get_yaxis().set_visible(False)
    fig.subplots_adjust(top=1.0, bottom=0.0, right=1.0, left=0.0)
    pyplot.savefig(saveas, dpi=dpi)


def plot_deghist(deghist, saveas, overwrite=False, **kwargs):
    """ Plot degrees of graph nodes.

    TODO
    ----
    Document keyword arguments.

    Parameters
    ----------
    deghist : array-like
        Source graph degree histogram intensities.
    saveas : str
        Path to save file. Returns `None` if file exists.
    overwrite : bool, optional
        If `True` allows plot to be overwritten. Default is `False`.
    """

    from numpy import arange
    from matplotlib import pyplot
    from matplotlib.ticker import MaxNLocator

    if not overwrite:
        from os.path import exists
        assert not exists(saveas), f' Graph cannot overwrite: {saveas}'

    cdef tuple figsize = kwargs.get('figsize', (6, 6))
    cdef str color = kwargs.get('cmap', 'k')
    cdef int dpi = kwargs.get('dpi', 300)
    cdef dict xlim = kwargs.get('xlim', None)

    pyplot.close('all')
    fig, ax = pyplot.subplots(figsize=figsize)
    ax.bar(arange(0, len(deghist), 1), deghist, color=color)
    ax.set_xlabel('Degree')
    ax.set_ylabel('Counts')
    if xlim is not None:
        ax.set_xlim(**xlim)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    pyplot.tight_layout()
    pyplot.savefig(saveas, dpi=dpi)
