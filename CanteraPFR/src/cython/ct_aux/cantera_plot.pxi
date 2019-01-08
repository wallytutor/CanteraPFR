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
