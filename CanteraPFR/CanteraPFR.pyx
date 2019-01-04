# -*- coding: utf-8 -*-

from libcpp.string cimport string
from CanteraPFR cimport AdiabaticPFR
from CanteraPFR cimport HeatWallPFR
from CanteraPFR cimport IsothermalPFR
from CanteraPFR cimport SolvePFR

import os
import time
import ctypes
from numpy import arange
from numpy import array
from pandas import DataFrame

#

cdef class Closure:
    """ Just-in-time compiler for Python functions.

    Allows a Python function to be called from C/C++ library.
    Based on https://stackoverflow.com/questions/51044122.

    TODO
    ----
        This class is incoherent, once the interface of `_inner_fun` is fixed
        and does not follow what `ftype` requires. This must be wrapped.

    Parameters
    ----------
    python_fun : function
        Python function to call from C/C++.
    ftype : ctypes.CFUNCTYPE
        Interface to provide the function.
    """

    cdef object python_fun
    cdef object jit_wrap

    def __cinit__(self, python_fun, ftype):
        self.python_fun = python_fun
        self.jit_wrap = ftype(lambda *args: self.python_fun(*args))

    cdef func_t get_fun_ptr(self):
        return (<func_t *><size_t>ctypes.addressof(self.jit_wrap))[0]


cdef class PyPFR:
    """ Python wrapper to plug-flow reactor models.

    This class is built around several plug-flow reactor models.  The selection
    of these is made according to the provided model name:

        - `isothermal`: isothermal PFR (temperature held constant).
        - `adiabatic`: adiabatic PFR (only viscous losses).
        - `heatwall`: PFR with imposed wall temperature profile.

    TODO
    ----
    1. Deal with Python strings.
    2. Make standard graphical ouput.
    3. Add interface to user viscosity function.
    4. Make htc a function of gas speed.
    5. Allow vector of absolute tolerances.

    Parameters
    ----------
    rtype : str
        Key-word for model selection.
    mech : str
        Path to CTI/XML mechanism file.
    phase : str
        Phase name in mechanism file.
    Di : float
        Reactor diameter in meters.
    T0 : float
        Inlet temperature in kelvin.
    p0 : float
        Inlet pressure in pascal.
    X0 : str
        Inlet composition compatible with Cantera in string format.
    Q0 : float
        Volume flow rate at reference state (assumed to be 273.15 K and
        101325 Pa) given in cubic centimeters per minute. This unit was chosen
        because it is the most common used with laboratory scale mass flow
        controllers.
    htc : float, optional
        Wall convective heat transfer coefficient in watts per kelvin per square
        meter. If its value is zero, this is equivalent to an adiabatic reactor,
        regardless of wall temperature. This parameter is intended to be used
        with Tw. Must be supplied if `rtype='heatwall'`. Default is None.
    Tw : float or function `float(float)`, optional
        Wall temperature profile in kelvin. This can be supplied as a constant
        value or as a function of position along reactor axis (start coordinate
        at zero). Must be supplied if `rtype='heatwall'`. Default is None.

    Raises
    ------
    SystemExit
        If reactor type name is not recognized or if provided wall temperature
        function is not suitable for use with reactor model `heatwall`.
    """

    cdef CanteraPFR* obj
    cdef SolvePFR* sol
    cdef Closure cl

    def __cinit__(self, rtype, mech, phase, Di, T0, p0, X0, Q0, htc=None,
                  Tw=None):
        # cdef bytearray cpp_mech = bytearray(mech, 'utf8')
        # cdef string cpp_mech = <string> mech.encode('utf-8')

        print('*' * 79)
        if rtype.lower() == 'isothermal':
            self.obj = new IsothermalPFR(mech, phase, Di, T0, p0, X0, Q0)

        elif rtype.lower() == 'adiabatic':
            self.obj = new AdiabaticPFR(mech, phase, Di, T0, p0, X0, Q0)

        elif rtype.lower() == 'heatwall':
            assert htc is not None, 'Missing heat transfer coefficient'
            assert Tw is not None, 'Missing wall temperature'
            self.cl = self._get_closure(Tw)
            self.obj = new HeatWallPFR(mech, phase, Di, T0, p0, X0, Q0, htc,
                                       self.cl.get_fun_ptr())

        else:
            raise SystemExit(f'Unknown reactor type : {rtype}')

        self.sol = new SolvePFR(<CanteraPFR *> self.obj)
        self.sol.init(0.0)

    def __dealloc__(self):
        if type(self) is PyPFR:
            del self.obj
            del self.sol

    def _get_closure(self, Tw):
        ftype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)

        if isinstance(Tw, float):
            return Closure(lambda x: Tw, ftype)

        try:
            Tw(0.0)
        except Exception as err:
            raise SystemExit(f'Provide function is not suitable : {err}')

        return Closure(Tw, ftype)

    def set_tolerances(self, rtol, atol):
        """ Set solution tolerances.

        Parameters
        ----------
        rtol : float
            Relative tolerance applied to all variables.
        atol : float
            Absolute tolerance applied to all variables.
        """
        self.sol.setTolerances(rtol, atol)

    def set_max_num_steps(self, num):
        """ Set maximum number of integration steps.

        This represents the maximum number of internal steps taken by IDA while
        integrating the problem, not the number of outputs made by the user.

        Parameters
        ----------
        num : int
            Maximum number of integration steps.
        """
        self.sol.setMaxNumSteps(num)

    def set_initial_step_size(self, h0):
        """ Set length of initial integration step.

        Notice that a value too small here may lead to failure if the maximum
        number of steps as provided by :meth:`set_max_num_steps` is not
        consistent. Do not use this method if you are not facing failures at
        integration start or precision issues.

        Parameters
        ----------
        h0 : float
            Length of initial integration step.
        """
        self.sol.setInitialStepSize(h0)

    def set_stop_position(self, xstop):
        """ Set last coordinate to stop integration.

        Eventually IDA solver may trespass the last point to integrate.
        If using :meth:`manage_solution`, this is automatically done for you.

        Parameters
        ----------
        xstop : float
            Coordinate of last point in meters.
        """
        self.sol.setStopPosition(xstop)

    def solve(self, xout):
        """ Integrate for current time to the provided time.

        Notice that the solver will fail if `xout` is less than current step
        and no test is performed by this class.

        Parameters
        ----------
        xout : float
            Coordinate of point to return solution.

        Returns
        -------
        int
            State of solution as defined by IDA constants.
        """
        return self.sol.solve(xout)

    def solution(self, num):
        """ Returns solution component by its index.

        Solution array depends on the used model, so no general description can
        be provided here. For advanced usage only.

        Parameters
        ----------
        num : int
            Index of solution component to return.

        Returns
        -------
        float
            Value of solution component.
        """
        return self.sol.solution(num)

    def manage_solution(self, Lr, dx, saveas, outfreq=100):
        """ Manage integration of reactorself.

        Parameters
        ----------
        Lr : float
            Total length of reactor in meters.
        dx : float
            Step to save results.
        saveas : str
            Path to results file.
        outfreq : int, optional
            Count of steps to provide screen output. Default is `100`.

        Raises
        ------
        AssertionError
            If parameters are unphysical (negative lengths), if step is below
            total length, if output file exists, or negative output frequency.

        Returns
        -------
        pandas.DataFrame
            Table of results for further manipulation.
        """

        assert Lr > 0, 'Reactor length must be positive'
        assert dx > 0, 'Saving step must be positive'
        assert Lr > dx, 'Reactor length must be above saving step'
        assert not os.path.exists(saveas), 'Output file already exists'
        assert outfreq > 0, 'Output frequency must be positive'

        columns = self.sol.variablesNames()
        columns = ['x'] + [c.decode('utf-8') for c in columns]
        sol0 = array([0] + list(self.sol.solutionVector()))

        # TODO convert to mole fractions!
        data = DataFrame(columns=columns)
        data.loc[0, columns] = sol0

        saveat = arange(0.0, Lr+dx/2, dx)
        self.sol.setStopPosition(Lr)

        t0 = time.time()
        for i, x in enumerate(saveat[1:]):
            if not i % outfreq:
                print(f'Currently solving to reach {x:.4f} m')

            if self.sol.solve(x) == -99:
                raise SystemExit('Failure during integration')

            sol0 = array([x] + list(self.sol.solutionVector()))
            data.loc[i+1, columns] = sol0

        print(f'Integration took {time.time()-t0:.5f} s')
        data.to_csv(saveas, index=False)

        return data
