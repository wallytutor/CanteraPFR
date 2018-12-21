# -*- coding: utf-8 -*-
"""
CanteraPFR base PFR model.

Implements a common PFR interface.
"""

import abc
import time
import numpy
import cantera
from scikits.odes import dae


class BaseSolution(abc.ABC):
    """ Common interface for solution retrieval.

    This interface specifies the minimum format for output of results. That
    means that all PFR classes should provide at least access to original IDA
    solution object, a positions array, a means to write a CSV file, and a
    method to retrieve the mole fractions of species. It assume that each of
    the solutions arrays start by the species. Inherited classes shall deal
    with what comes after the species. By default only mole fractions are
    recovered, but users can translate back to mass fractions by accessing the
    associated gas phase object.

    Parameters
    ----------
    gas : cantera.Solution
        The gas phase object used during integration.
    solution : namedtuple
        A namedtuple with attributes given in the table. For more details
        check `scikits.odes.dae`.

        =========== ==========================================
        Field       Meaning
        =========== ==========================================
        ``flag``    An integer flag (StatusEnumXXX)
        ``values``  Named tuple with entries array t and array y and
                    array ydot. y will correspond to y_retn value and
                    ydot to yp_retn!
        ``errors``  Named tuple with entries t and y and ydot of error
        ``roots``   Named tuple with entries array t and array y and
                    array ydot
        ``tstop``   Named tuple with entries array t and array y and
                    array ydot
        ``message`` String with message in case of an error
        =========== ==========================================
    """

    def __init__(self, gas, solution):
        self._gas, self._idasol = gas, solution
        self._position = self._idasol.values.t

        for i, Y in enumerate(self._idasol.values.y[:, :]):
            self._gas.TPY = 273.15, 101325, Y[:self._gas.n_species]
            self._idasol.values.y[i, :self._gas.n_species] = self._gas.X

    @abc.abstractmethod
    def to_csv(self, saveas):
        """ Standard method for CSV output. """
        pass

    def mole_fraction(self, name):
        """ Retrieve species mole fraction along reactor.

        Parameters
        ----------
        name : str
            Species name as stated in mechanism file.

        Returns
        -------
        numpy.array
            Array of mole fractions at each cell over domain.
        """

        idx = self._gas.species_index(name)
        return self._idasol.values.y[:, idx]

    @property
    def idasol(self):
        """ Original namedtuple returned by IDA solver. """

        return self._idasol

    @property
    def position(self):
        """ Positions of cells centers in meters. """

        return self._position

    @property
    def gas(self):
        """" Gas phase cantera.Solution object used for solution. """

        return self._gas


class BasePFR(abc.ABC):
    """ Common interface for plug-flow reactors.

    This base class provides a common interface to all PFR's in this package,
    allowing for a single interface for solution management.  It implies some
    basic requirements for derived classes to respect. First, all abstract
    methods must be supplied, otherwise interpreter will be unable to
    instantiate a class. Creation of any Cantera object shall be left to this
    base class, avoiding import of Cantera in any other module.

    Parameters
    ----------
    mech : str
        Path to gas phase mechanism in CTI or CTML format.
    phase : str
        Phase name in mechanism file.
    trans : str
        Transport model to use in mechanism file.
    """

    def __init__(self, mech, phase, trans, solmgr):

        self._solmgr = solmgr
        self._gas_constant = cantera.gas_constant

        try:
            self._gas = cantera.Solution(mech, phase, trans=trans)
        except:
            # TODO this is mostly for debug, Cantera deals internally!!
            print('Failed : phase {} transport {}'.format(phase, trans))
            self._gas = cantera.Solution(mech, phase, trans=None)

    def _initialize_viscosity(self, mufunc):
        """ Iniailize viscosity model. """

        print('\n Selecting viscosity model')

        if mufunc is not None:
            print('   Using user supplied viscosity')
            self._viscosity = mufunc
            self._viscosity(298.15)
        else:
            try:
                def mufunc(T):
                    return self._gas.viscosity

                mufunc(298.15)
                print('   Using viscosity from mechanism')
            except cantera.CanteraError:
                print('   Using viscosity fallback function')
                def mufunc(T):
                    return 3.957996309582866e-05

            self._viscosity = mufunc

    @abc.abstractmethod
    def report_derivatives(self, when):
        """ Common interface for reporting derivatives.

        This function shall be inherited by all classes and provides an output
        of space derivatives vector at a given time.
        """
        pass

    def integrate(self, L, d, **kwargs):
        """ Integrate model in space.

        Provides integration of the differential algebraic set of equations.
        Keyword arguments are parsed to `scikits.odes.dae` and must be
        compatible with solver `ida`.

        Note
        ----
        Do not provide keyword argument `old_api` to `ida`.

        Parameters
        ----------
        L : float
            Lenght of plug-flow reactor to integrate in meters.
        d : float
            Cell length for discretization in meters.

        Returns
        -------
        BaseSolution
            A solution object. See inherited class for details because the
            interface may change depending on the model.
        """

        t0 = time.time()
        coords = numpy.arange(0, L, d)
        self.report_derivatives('Initial')
        solver = dae('ida', self, old_api=False, **kwargs)
        solution = solver.solve(coords, self._vec0, self._vecp0)
        self.report_derivatives('Final')
        print('\n Solution took {:.6e} s'.format(time.time()-t0))
        return self._solmgr(self._gas, solution)
