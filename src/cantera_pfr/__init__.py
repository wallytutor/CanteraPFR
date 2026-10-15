# -*- coding: utf-8 -*-
import os
from pathlib import Path
from .dimensionless import Dimensionless, Gr, Pe, Pr, Ra, Re, Sc
from .drg import SimplifyDRG, SimplifyDRGManager, SimplifySetup, StateSpace
from .graph import (
    BaseGraph,
    DirectedGraph,
    UndirectedGraph,
    analyse_graph
)
from .mechanism import MechanismFilter, filter_mechanism
from .plot import (
    Plotter,
    plot_adjmatrix,
    plot_deghist,
    plot_pfr,
    plotPFR
)

try:
    from .ct_pfr import PyPFR
    from .cpfr import CPFR, cPFR
except (ImportError, Exception):
    pass

__version__ = '0.2.0'


def set_env() -> Path:
    """ Configure environment variables and Cantera search paths.

    Returns
    -------
    Path
        Path to the package data directory.
    """
    data_dir = Path(__file__).resolve().parent / 'data'

    if data_dir.exists():
        data_str = str(data_dir)
        current = os.environ.get('CANTERA_DATA', '')
        os.environ['CANTERA_DATA'] = (
            f'{data_str}:{current}' if current else data_str
        )

        try:
            import cantera as ct
            ct.add_directory(data_str)
            ct.suppress_thermo_warnings()
        except (ImportError, Exception):
            pass

    return data_dir


set_env()

__all__ = [
    'Dimensionless',
    'Re',
    'Pr',
    'Sc',
    'Pe',
    'Gr',
    'Ra',
    'MechanismFilter',
    'filter_mechanism',
    'Plotter',
    'plot_pfr',
    'plotPFR',
    'plot_adjmatrix',
    'plot_deghist',
    'BaseGraph',
    'DirectedGraph',
    'UndirectedGraph',
    'analyse_graph',
    'StateSpace',
    'SimplifySetup',
    'SimplifyDRG',
    'SimplifyDRGManager',
    'PyPFR',
    'CPFR',
    'cPFR',
    'set_env',
    '__version__'
]
