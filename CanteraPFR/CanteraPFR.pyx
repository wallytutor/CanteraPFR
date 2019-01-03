# -*- coding: utf-8 -*-
from libcpp.string cimport string
from CanteraPFR cimport AdiabaticPFR
from CanteraPFR cimport HeatWallPFR
from CanteraPFR cimport IsothermalPFR
from CanteraPFR cimport SolvePFR

import time
import ctypes
from numpy import arange
from numpy import array
from pandas import DataFrame

# See https://stackoverflow.com/questions/51044122
# https://github.com/JohannesBuchner/PyMultiNest/issues/5

cdef class Closure:
    cdef object python_fun
    cdef object jit_wrap

    def __cinit__(self, python_fun):
        self.python_fun = python_fun
        ftype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
        self.jit_wrap = ftype(self.inner_fun)

    def inner_fun(self, double arg):
        return self.python_fun(arg)

    cdef func_t get_fun_ptr(self):
        return (<func_t *><size_t>ctypes.addressof(self.jit_wrap))[0]


cdef class PyPFR:
    """ Python wrapper to plug-flow reactor models.

    This class is built around several plug-flow reactor models.  The selection
    of these is made according to the provided model name.

    TODO
    ----
    Provide docstrings to all methods.
    Make standard graphical ouput.
    Include theory in rst documentation.

    Parameters
    ----------
    """
    cdef CanteraPFR* obj
    cdef SolvePFR* sol
    cdef Closure cl

    def __cinit__(self, rtype, mech, phase, Di, T0, p0, X0, Q0, htc=None,
                  Tw=None):
        # TODO deal with Python strings!
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

            if isinstance(Tw, float):
                self.cl = Closure(lambda x: Tw)
            else:
                # Some test here!
                self.cl = Closure(Tw)

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

    def set_tolerances(self, rtol, atol):
        self.sol.setTolerances(rtol, atol)

    def set_max_num_steps(self, num):
        self.sol.setMaxNumSteps(num)

    def set_initial_step_size(self, h0):
        self.sol.setInitialStepSize(h0)

    def set_stop_position(self, tstop):
        self.sol.setStopPosition(tstop)

    def solve(self, tout):
        return self.sol.solve(tout)

    def solution(self, num):
        return self.sol.solution(num)

    def manage_solution(self, Lr, dx, saveas, outfreq=100):
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
