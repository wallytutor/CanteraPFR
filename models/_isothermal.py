# -*- coding: utf-8 -*-
"""
CanteraPFR isothermal PFR model.

Implements an isothermal PFR.
"""

import numpy
import cantera
from scikits.odes import dae


class IsothermalPFR(object):
    """ Isothermal plug-flow reactor model.

    Provides an implementation of an isothermal plug-flow reactor built on top
    of Cantera. This class creates a gas phase, which is initialized to the
    given conditions and used to setup the DAE system to be solved.

    Parameters
    ----------
    mech : str
        Path to gas phase mechanism in CTI or CTML format. If phase name in
        file is not called 'gas', this should be specified by keyword argument
        `phase`. By default a `trans='multi'` will be used to retrieve the
        transport model from `cantera.Solution`. This can be overriden through
        keyword argument `trans`.
    TPX : tuple
        Tuple of temperature (in kelvin), pressure (in pascal), and gas molar
        composition in a format compatible with `cantera.Solution.TPX`. These
        values are used to establish the inlet boundary condition and compute
        mass flow rate.
    Q : float
        Volume flow rate at reference state (assumed to be 273.15 K and
        101325 Pa) given in cubic centimeters per minute. This unit was chosen
        because it is the most common used with laboratory scale mass flow
        controllers. Values of reference state variables can be overriden by
        keywords `T_ref` for temperature (in kelvin) and `P_ref` for pressure
        (in pascal).
    Ac : float
        Reactor cross-section area in square meters. The model assumes that
        this section is circular for momentum equation purposes.
    mufunc : function, optional
        If the model does not include transport data, the user may supply a
        function of temperature to model dynamic viscosity. Notice that the
        keyword argument `trans` shall be provided with value `None` in this
        case. Notice that currently when a mechanism fails to get the viscosity
        it will assume pure nitrogen. Default is `None`.
    """

    def __init__(self, mech, TPX, Q, Ac, mufunc=None, **kwargs):

        phase = kwargs.get('phase', 'gas')
        trans = kwargs.get('trans', 'multi')
        T_ref = kwargs.get('T_ref', 273.15)
        P_ref = kwargs.get('P_ref', cantera.one_atm)
        T0, P0, X0, self._Ac = *TPX, Ac

        self._gas = cantera.Solution(mech, phase, trans=trans)
        self._gas.TPX = T_ref, P_ref, X0
        self._rho_ref = self._gas.density
        self._Wk = self._gas.molecular_weights
        self._size = 3 + self._gas.n_species

        self._initialize_mu_func(mufunc)
        self._initialize_dae_problem(Q, T0, P0, X0)

    def _initialize_mu_func(self, mufunc):
        """ Iniailize viscosity model. """

        if mufunc is not None:
            self._mufunc = mufunc
        else:
            def mufunc(T):
                return 3.957996309582866e-05
            self._mufunc = mufunc

    def _initialize_dae_problem(self, Q, T0, P0, X0):
        """ Apply constraints to initialize derivatives. """

        self._gas.TPX = T0, P0, X0

        try:
            mu = self._gas.viscosity
        except cantera.CanteraError:
            mu = self._mufunc(self._gas.T)

        # Make code closer to mathematical formulation.
        rho0 = self._gas.density
        u0 = Q * self._rho_ref / (6.0e+07 * rho0 * self._Ac)
        Y0 = self._gas.Y
        RT = cantera.gas_constant * T0
        Wavg = self._gas.mean_molecular_weight
        wdot = self._gas.net_production_rates

        A = numpy.zeros((self._size, self._size))
        b = numpy.zeros(self._size)

        # Rows corresponding to species equations.
        for i in range(self._gas.n_species):
            A[i, 1+i] = rho0 * u0  # TODO why 1+i? I forgot!
        b[:-3] = wdot * self._Wk

        # Row of A corresponding to conservation of mass equation.
        A[-3, 0:2] = numpy.hstack((rho0, u0))
        b[-3] = 0

        # Row corresponding to momentum equation.
        A[-2, 0] = rho0 * u0
        A[-2, -1] = 1
        b[-2] = -8 * u0 * mu * numpy.pi / self._Ac

        # Row corresponding to state equation.
        coefYp = P0 * Wavg ** 2 / self._Wk
        A[-1, :] = numpy.hstack((0, RT, coefYp, -Wavg))
        b[-1] = 0

        self._vec0 = numpy.hstack((Y0, u0, rho0, P0))
        self._vecp0 = numpy.linalg.solve(A, b)

    def _unpack_ida(self, vec, vecp):
        """ Unpack IDA arrays into readable variables.

        Array of state `vec` and array of derivatives `vecp` are organized in
        the same order of physical parameters: velocity, density, species and
        pressure. The positions are provided by the defined slices.
        """

        Y = vec[:-3]
        u = vec[-3]
        rho = vec[-2]
        p = vec[-1]

        dYdz = vecp[:-3]
        dudz = vecp[-3]
        drhodz = vecp[-2]
        dpdz = vecp[-1]

        return (Y, u, rho, p), (dYdz, dudz, drhodz, dpdz)

    def __call__(self, z, vec, vecp, result):
        """ Residual equations for the problem.

        Parameters
        ----------
        z : float
            Independent variable (position).
        vec : numpy.array
            State of system expressed as [Y, u, rho, p].
        vecp : numpy.array
            Derivatives array [dYdz, dudz, drhodz, dpdz].
        result : numpy array
            Array of residuals (return of function).
        """

        (Y, u, rho, p), diffs = self._unpack_ida(vec, vecp)
        dYdz, dudz, drhodz, dpdz = diffs

        try:
            self._gas.set_unnormalized_mass_fractions(Y)
            self._gas.TP = self._gas.T, p
        except Exception as err:
            pass

        try:
            mu = self._gas.viscosity
        except cantera.CanteraError:
            mu = self._mufunc(self._gas.T)

        urho = u * rho
        result[:-3] = urho * dYdz - self._gas.net_production_rates * self._Wk
        result[-3] = rho * dudz + u * drhodz
        result[-2] = urho * dudz + dpdz + 8 * mu * u * numpy.pi / self._Ac
        result[-1] = self._gas.density - rho

    def integrate(self, L, d, **kwargs):
        """ Integrate model in space.

        Provides integration of the differential algebraic set of equations.
        Keyword arguments are parsed to `scikits.odes.dae` and must be
        compatible with solver `ida`.

        Note
        ----
        Do not provide keyword argument `old_api`.

        Parameters
        ----------
        L : float
            Lenght of plug-flow reactor to integrate in meters.
        d : float
            Cell length for discretization in meters.

        Returns
        -------
        Solution

            A wrapper class object with access to properties:

            - ``idasol``:  a namedtuple with attributes given in the table.
                For more details check `scikits.odes.dae`.

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

            - ``position``: Positions in meters corresponding to the center of
                the cells in discrete space.
            - ``velocity``: Axial velocity in meters per second (each cell).
            - ``density``: Density in kilograms per cubic meter (each cell).
            - ``pressure``: Pressu in pascal (each cell).

            Additionally, this object has a method `mole_fraction` which gets
            the name of a species as argument and returns an array of its mole
            fractions in the direction of reactor axis for each cell.
        """

        coords = numpy.arange(0, L, d)
        solver = dae('ida', self, old_api=False, **kwargs)
        solution = solver.solve(coords, self._vec0, self._vecp0)

        for i, Y in enumerate(solution.values.y[:, :]):
            self._gas.TPY = 273.15, 101325, Y[:self._gas.n_species]
            solution.values.y[i, :self._gas.n_species] = self._gas.X

        class Solution:
            idasol = solution
            position = solution.values.t
            velocity = solution.values.y[:, -3]
            density = solution.values.y[:, -2]
            pressure = solution.values.y[:, -1]

            def mole_fraction(this, name):
                # WARNING: using `this` instead of `self`!
                idx = self._gas.species_index(name)
                return this.idasol.values.y[:, idx]

        return Solution()
