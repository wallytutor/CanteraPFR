# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Mapping, Sequence
from matplotlib import pyplot
from matplotlib.ticker import MaxNLocator
from numpy import arange
import numpy as np
import pandas as pd


class Plotter:
    """ Plotting utility methods for reactors and graphs.

    Parameters
    ----------
    style : str = 'bmh'
        Matplotlib style sheet name.

    Returns
    -------
    Plotter
        Plotter instance.
    """

    __slots__ = ('_style',)

    def __init__(self, style: str = 'bmh') -> None:
        # XXX no documentation here, this is intentional
        self._style = style

    def plot_pfr(
        self,
        data: pd.DataFrame,
        species: Mapping[str, str],
        saveas: Path | str,
        dpi: int = 300,
        figsize: tuple[int, int] = (10, 6),
        axes: Mapping[str, str] | None = None
    ) -> None:
        """ Plot PFR integration profiles.

        Parameters
        ----------
        data : pd.DataFrame
            Table of reactor solution data.
        species : Mapping[str, str]
            Species column names and display labels.
        saveas : Path | str
            Destination image file path.
        dpi : int = 300
            Output image resolution.
        figsize : tuple[int, int] = (10, 6)
            Figure dimensions (width, height) in inches.
        axes : Mapping[str, str] | None = None
            Custom axis label mapping.

        Returns
        -------
        None
            Generates plot file on disk.
        """
        axes_default = {
            'x': 'Position ($\\mathrm{m}$)',
            'T': 'Temperature ($\\mathrm{K}$)',
            'rho': 'Density ($\\mathrm{kg\\,m^{3}}$)',
            'u': 'Velocity ($\\mathrm{m\\,s^{-1}}$)',
            'p': 'Pressure ($\\mathrm{Pa}$)',
            'Y': 'Mass fraction'
        }
        axis_map = dict(axes_default)

        if axes is not None:
            axis_map.update(axes)

        xaxis = 'x'
        uaxis = 'u'
        paxis = 'p'
        yaxis = 'Y'
        taxis = 'T' if 'T' in data.columns else 'rho'

        pyplot.close('all')
        pyplot.style.use(self._style)
        pyplot.figure(figsize=figsize)

        pyplot.subplot(221)
        pyplot.plot(data[xaxis], data[taxis])
        pyplot.xlabel(axis_map[xaxis])
        pyplot.ylabel(axis_map[taxis])

        pyplot.subplot(222)
        pyplot.plot(data[xaxis], data[uaxis])
        pyplot.xlabel(axis_map[xaxis])
        pyplot.ylabel(axis_map[uaxis])

        pyplot.subplot(223)
        pyplot.plot(data[xaxis], data[paxis])
        pyplot.xlabel(axis_map[xaxis])
        pyplot.ylabel(axis_map[paxis])

        pyplot.subplot(224)

        for key, spec_label in species.items():
            if key in data.columns:
                pyplot.plot(data[xaxis], data[key], label=spec_label)

        pyplot.xlabel(axis_map[xaxis])
        pyplot.ylabel(axis_map[yaxis])
        pyplot.legend()

        pyplot.tight_layout()
        out_path = Path(saveas)
        pyplot.savefig(out_path, dpi=dpi)

    def plot_adjmatrix(
        self,
        adjmat: np.ndarray,
        saveas: Path | str,
        overwrite: bool = True,
        dpi: int = 300,
        figsize: tuple[int, int] = (6, 6),
        cmap: str = 'Greys'
    ) -> None:
        """ Plot graph adjacency matrix.

        Parameters
        ----------
        adjmat : np.ndarray
            Graph adjacency matrix.
        saveas : Path | str
            Destination image file path.
        overwrite : bool = True
            Allow overwriting existing output file.
        dpi : int = 300
            Output image resolution.
        figsize : tuple[int, int] = (6, 6)
            Figure dimensions in inches.
        cmap : str = 'Greys'
            Colormap name.

        Returns
        -------
        None
            Generates plot file on disk.
        """
        out_path = Path(saveas)

        if not overwrite and out_path.exists():
            raise FileExistsError(f'Graph cannot overwrite: {out_path}')

        pyplot.close('all')
        fig, ax = pyplot.subplots(figsize=figsize)
        ax.imshow(adjmat, cmap=cmap)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        fig.subplots_adjust(top=1.0, bottom=0.0, right=1.0, left=0.0)
        pyplot.savefig(out_path, dpi=dpi)

    def plot_deghist(
        self,
        deghist: Sequence[int | float],
        saveas: Path | str,
        overwrite: bool = False,
        dpi: int = 300,
        figsize: tuple[int, int] = (6, 6),
        color: str = 'k',
        xlim: tuple[float, float] | None = None
    ) -> None:
        """ Plot degrees histogram of graph nodes.

        Parameters
        ----------
        deghist : Sequence[int | float]
            Histogram count sequence indexed by degree.
        saveas : Path | str
            Destination image file path.
        overwrite : bool = False
            Allow overwriting existing output file.
        dpi : int = 300
            Output image resolution.
        figsize : tuple[int, int] = (6, 6)
            Figure dimensions in inches.
        color : str = 'k'
            Bar plot color.
        xlim : tuple[float, float] | None = None
            X-axis limits (min, max).

        Returns
        -------
        None
            Generates plot file on disk.
        """
        out_path = Path(saveas)

        if not overwrite and out_path.exists():
            raise FileExistsError(f'Graph cannot overwrite: {out_path}')

        pyplot.close('all')
        fig, ax = pyplot.subplots(figsize=figsize)
        x_coords = arange(0, len(deghist), 1)
        ax.bar(x_coords, deghist, color=color)
        ax.set_xlabel('Degree')
        ax.set_ylabel('Counts')

        if xlim is not None:
            ax.set_xlim(xlim[0], xlim[1])

        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        pyplot.tight_layout()
        pyplot.savefig(out_path, dpi=dpi)


def plot_pfr(
        data: pd.DataFrame,
        species: Mapping[str, str],
        saveas: Path | str,
        **kwargs: object
    ) -> None:
    """ Plot PFR integration profiles.

    Parameters
    ----------
    data : pd.DataFrame
        Table of reactor solution data.
    species : Mapping[str, str]
        Species column names and display labels.
    saveas : Path | str
        Destination image file path.
    kwargs : object
        Additional plotting keyword arguments.

    Returns
    -------
    None
        Generates plot file on disk.
    """
    style = str(kwargs.get('style', 'bmh'))
    dpi = int(kwargs.get('dpi', 300))
    figsize = kwargs.get('figsize', (10, 6))
    axes = kwargs.get('axes', None)
    Plotter(style=style).plot_pfr(
        data,
        species,
        saveas,
        dpi=dpi,
        figsize=figsize,
        axes=axes
    )


def plotPFR(
        data: pd.DataFrame,
        species: Mapping[str, str],
        saveas: Path | str,
        **kwargs: object
    ) -> None:
    """ Plot PFR integration profiles.

    Parameters
    ----------
    data : pd.DataFrame
        Table of reactor solution data.
    species : Mapping[str, str]
        Species column names and display labels.
    saveas : Path | str
        Destination image file path.
    kwargs : object
        Additional plotting keyword arguments.

    Returns
    -------
    None
        Generates plot file on disk.
    """
    plot_pfr(data, species, saveas, **kwargs)


def plot_adjmatrix(
        adjmat: np.ndarray,
        saveas: Path | str,
        overwrite: bool = True,
        **kwargs: object
    ) -> None:
    """ Plot graph adjacency matrix.

    Parameters
    ----------
    adjmat : np.ndarray
        Graph adjacency matrix.
    saveas : Path | str
        Destination image file path.
    overwrite : bool = True
        Allow overwriting existing output file.
    kwargs : object
        Additional plotting keyword arguments.

    Returns
    -------
    None
        Generates plot file on disk.
    """
    dpi = int(kwargs.get('dpi', 300))
    figsize = kwargs.get('figsize', (6, 6))
    cmap = str(kwargs.get('cmap', 'Greys'))
    Plotter().plot_adjmatrix(
        adjmat,
        saveas,
        overwrite=overwrite,
        dpi=dpi,
        figsize=figsize,
        cmap=cmap
    )


def plot_deghist(
        deghist: Sequence[int | float],
        saveas: Path | str,
        overwrite: bool = False,
        **kwargs: object
    ) -> None:
    """ Plot degrees histogram of graph nodes.

    Parameters
    ----------
    deghist : Sequence[int | float]
        Histogram count sequence indexed by degree.
    saveas : Path | str
        Destination image file path.
    overwrite : bool = False
        Allow overwriting existing output file.
    kwargs : object
        Additional plotting keyword arguments.

    Returns
    -------
    None
        Generates plot file on disk.
    """
    dpi = int(kwargs.get('dpi', 300))
    figsize = kwargs.get('figsize', (6, 6))
    color = str(kwargs.get('color', kwargs.get('cmap', 'k')))
    xlim = kwargs.get('xlim', None)
    Plotter().plot_deghist(
        deghist,
        saveas,
        overwrite=overwrite,
        dpi=dpi,
        figsize=figsize,
        color=color,
        xlim=xlim
    )
