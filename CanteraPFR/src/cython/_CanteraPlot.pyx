# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt


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
            Output resolution (parsed to plt.savefig). Default is 300.
    """

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

    plt.close('all')
    plt.style.use(style)
    plt.figure(figsize=figsize)

    plt.subplot(221)
    plt.plot(data[xaxis], data[Taxis])
    plt.xlabel(axes[xaxis])
    plt.ylabel(axes[Taxis])

    plt.subplot(222)
    plt.plot(data[xaxis], data[uaxis])
    plt.xlabel(axes[xaxis])
    plt.ylabel(axes[uaxis])

    plt.subplot(223)
    plt.plot(data[xaxis], data[paxis])
    plt.xlabel(axes[xaxis])
    plt.ylabel(axes[paxis])

    plt.subplot(224)
    for key, spec in species.items():
        plt.plot(data[xaxis], data[key], label=spec)
    plt.xlabel(axes[xaxis])
    plt.ylabel(axes[Yaxis])
    plt.legend()

    plt.tight_layout()
    plt.savefig(saveas, dpi=dpi)
