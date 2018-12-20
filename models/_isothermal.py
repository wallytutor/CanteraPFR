# -*- coding: utf-8 -*-
"""
CanteraPFR isothermal PFR model.

Implements an isothermal PFR.
"""

import time
import numpy
import pandas
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
    ofreq : int, optional
        Interval of calls to output advance to screen. Default is `100`.
    """

    def __init__(self, mech, TPX, Q, Ac, mufunc=None, ofreq=100, **kwargs):

        phase = kwargs.get('phase', 'gas')
        trans = kwargs.get('trans', 'multi')
        T_ref = kwargs.get('T_ref', 273.15)
        P_ref = kwargs.get('P_ref', cantera.one_atm)
        T0, P0, X0, self._Ac, self._ofreq = *TPX, Ac, ofreq

        print(' Starting IsothermalPFR\n' +
              ' Mechanism ............. {}\n'.format(mech) +
              ' Phase ................. {}\n'.format(phase) +
              ' Transport ............. {}\n'.format(trans) +
              ' Reference temperature . {} K\n'.format(T_ref) +
              ' Reference pressure .... {} Pa\n'.format(P_ref) +
              ' Initial temperature ... {} K\n'.format(T0) +
              ' Initial pressure ...... {} Pa\n'.format(P0) +
              ' Initial composition ... {}\n'.format(X0) +
              ' Reactor cross-section . {} m\u00B2'.format(Ac)
              )

        try:
            self._gas = cantera.Solution(mech, phase, trans=trans)
        except:
            # TODO this is mostly for debug, Cantera deals internally!!
            print('Failed : phase {} transport {}'.format(phase, trans))
            self._gas = cantera.Solution(mech, phase, trans=None)

        self._gas.TPX = T_ref, P_ref, X0
        self._rho_ref = self._gas.density
        self._Wk = self._gas.molecular_weights
        self._size = 3 + self._gas.n_species
        self._ocount = 0

        self._initialize_mu_func(mufunc)
        self._initialize_dae_problem(Q, T0, P0, X0)

    def _initialize_mu_func(self, mufunc):
        """ Iniailize viscosity model. """

        if mufunc is not None:
            print(' Using user supplied viscosity')
            self._mufunc = mufunc
        else:
            try:
                def mufunc(T):
                    return self._gas.viscosity

                mufunc(298.15)
                print(' Using viscosity from mechanism')
            except cantera.CanteraError:
                print(' Using viscosity fallback function')
                def mufunc(T):
                    return 3.957996309582866e-05

            self._mufunc = mufunc

    def _initialize_dae_problem(self, Q, T0, P0, X0):
        """ Apply constraints to initialize derivatives.

        Compute initial derivatives and assembly initial state. Take care with
        indices: top-left corner start with species coefficients and last three
        rows are reserved for continuity, momentum and state equations. The
        following diagram tries to illustrate the shape of contraint matrix.
        For sake of clarity, we will fill even the zeros in last three rows so
        that comments can illustrate properly the derivatives they multiply.

           | m00  .   . . .     0      0    0 | Y0'  |   | rt0 * W0 |
           | 0    m11 . . .     0      0    0 | Y1'  |   | rt1 * W1 |
           | .    .   . . .     0      0    0 | .    |   | .        |
           | .    .   . . .     0      0    0 | .    | = | .        |
           | 0    .   . . mkk   0      0    0 | Yk'  |   | rtk * Wk |
           | 0    .   . . .     rho    u    0 | u'   |   | 0        |
           | 0    .   . . .     rho*u  0    1 | rho' |   | visc.    |
           | cY0  cY1 . . cYk   0      RT  -W | p'   |   | 0        |
        """

        self._gas.TPX = T0, P0, X0

        # Make code closer to mathematical formulation.
        rho0 = self._gas.density
        u0 = Q * self._rho_ref / (6.0e+07 * rho0 * self._Ac)
        Y0 = self._gas.Y
        mu = self._mufunc(self._gas.T)
        RT = cantera.gas_constant * T0
        Wavg = self._gas.mean_molecular_weight
        wdot = self._gas.net_production_rates

        # Create matrix/array for initial state.
        A = numpy.zeros((self._size, self._size))
        b = numpy.zeros(self._size)

        # Source terms.
        b[:-3] = wdot * self._Wk                     # species
        b[-3] = 0                                    # RHS continuity
        b[-2] = -8 * u0 * mu * numpy.pi / self._Ac   # RHS momentum
        b[-1] = 0                                    # RHS state

        # Species equations elements.
        k = self._gas.n_species
        A[range(k), range(k)] = rho0 * u0       # Yk'

        # Continuity equation elements.
        A[-3, :-3] = 0                          # Yk'
        A[-3, -3] = rho0                        # u'
        A[-3, -2] = u0                          # rho'
        A[-3, -1] = 0                           # p'

        # Momentum equation elements.
        A[-2, :-3] = 0                          # Yk'
        A[-2, -3] = rho0 * u0                   # u'
        A[-2, -3] = 0                           # rho'
        A[-2, -1] = 1                           # p'

        # State equation elements.
        A[-1, :-3] = P0 * Wavg ** 2 / self._Wk  # Yk'
        A[-1, -2] = 0                           # u'
        A[-1, -2] = RT                          # rho'
        A[-1, -1] = -Wavg                       # p'

        self._vec0 = numpy.hstack((Y0, u0, rho0, P0))
        self._vecp0 = numpy.linalg.solve(A, b)

        print('\n Initial derivatives\n' +
              '   u\'    {:+.6e} m.s\u207B\u00B2\n'.format(self._vecp0[-3]) +
              '   rho\'  {:+.6e} kg.m\u207B\u2074\n'.format(self._vecp0[-2]) +
              '   p\'    {:+.6e} Pa.m\u207B\u00B9\n'.format(self._vecp0[-1])
              )

    def __call__(self, z, vec, vecp, result):
        """ Residual equations for the problem.

        Array of state `vec` and array of derivatives `vecp` are organized in
        the same order of physical parameters: species, velocity, density, and
        pressure. The positions are provided by the defined slices.

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

        if not self._ocount % self._ofreq:
            print(' Solving at z = {:.4e}'.format(z))
        self._ocount += 1

        Y = vec[:-3]
        u = vec[-3]
        rho = vec[-2]
        p = vec[-1]

        dYdz = vecp[:-3]
        dudz = vecp[-3]
        drhodz = vecp[-2]
        dpdz = vecp[-1]

        try:
            self._gas.set_unnormalized_mass_fractions(Y)
            self._gas.TP = self._gas.T, p
        except Exception as err:
            pass

        urho = u * rho
        mu = self._mufunc(self._gas.T)
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

        TODO
        ----
        Make a standalone class for `Solution`.

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
            fractions in the direction of reactor axis for each cell. It also
            provide a method `to_csv` that accepts a path as argument to the
            output CSV file.
        """

        t0 = time.time()
        coords = numpy.arange(0, L, d)
        solver = dae('ida', self, old_api=False, **kwargs)
        solution = solver.solve(coords, self._vec0, self._vecp0)
        print('\n Solution took {:.6e} s'.format(time.time()-t0))

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

            def to_csv(this, saveas):
                queue = ['u', 'rho', 'p']
                columns = ['x'] + list(self._gas.species_names) + queue
                x = solution.values.t.reshape((solution.values.t.shape[0], 1))
                stack = numpy.hstack((x, solution.values.y))
                data = pandas.DataFrame(stack, columns=columns)
                data.to_csv(saveas, index=False)

        return Solution()
