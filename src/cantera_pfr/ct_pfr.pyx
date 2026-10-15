# -*- coding: utf-8 -*-
# cython: language_level=3

__all__ = ['PyPFR']

import ctypes
import time
from pathlib import Path
import numpy as np
from pandas import DataFrame

from libcpp.string cimport string
from libcpp.vector cimport vector
from .ct_pfr cimport (
    CanteraPFR,
    ConstAreaPFR,
    AdiabaticPFR,
    HeatWallPFR,
    IsothermalPFR,
    SolvePFR,
    func_t,
)


cdef class Closure:
    """ JIT callback wrapper for Python functions called from C++. """

    cdef object python_fun
    cdef object jit_wrap

    def __cinit__(self, python_fun, ftype):
        self.python_fun = python_fun
        self.jit_wrap = ftype(lambda *args: self.python_fun(*args))

    cdef func_t get_fun_ptr(self):
        cdef size_t ptr = <size_t>ctypes.cast(
            self.jit_wrap, ctypes.c_void_p
        ).value
        return (<func_t>ptr)


cdef class PyPFR:
    """ Python wrapper to plug-flow reactor models.

    Parameters
    ----------
    rtype : str
        Reactor model type ('isothermal', 'adiabatic', 'heatwall').
    mech : str
        Path to mechanism YAML file.
    phase : str
        Phase name in mechanism file.
    Di : float
        Reactor diameter in meters.
    T0 : float
        Inlet temperature in kelvin.
    p0 : float
        Inlet pressure in pascal.
    X0 : str
        Inlet composition in Cantera format.
    Q0 : float
        Inlet volume flow rate in sccm.
    htc : float, optional
        Wall convective heat transfer coefficient in W/(m^2 K).
    Tw : float or callable, optional
        Wall temperature in kelvin (constant or function of position).
    """

    cdef CanteraPFR* obj
    cdef SolvePFR* sol
    cdef Closure cl

    def __cinit__(
        self,
        str rtype,
        str mech,
        str phase,
        double Di,
        double T0,
        double p0,
        str X0,
        double Q0,
        object htc=None,
        object Tw=None
    ):
        cdef string cpp_mech = mech.encode('utf-8')
        cdef string cpp_phase = phase.encode('utf-8')
        cdef string cpp_X0 = X0.encode('utf-8')
        cdef double c_Di = Di
        cdef double c_T0 = T0
        cdef double c_p0 = p0
        cdef double c_Q0 = Q0
        cdef double c_htc = 0.0

        r_type = rtype.lower()

        if r_type == 'isothermal':
            self.obj = <CanteraPFR*> new IsothermalPFR(
                cpp_mech, cpp_phase, c_Di, c_T0, c_p0, cpp_X0, c_Q0
            )
        elif r_type == 'adiabatic':
            self.obj = <CanteraPFR*> new AdiabaticPFR(
                cpp_mech, cpp_phase, c_Di, c_T0, c_p0, cpp_X0, c_Q0
            )
        elif r_type == 'heatwall':
            if htc is None:
                raise ValueError('Missing heat transfer coefficient')

            if Tw is None:
                raise ValueError('Missing wall temperature')

            c_htc = float(htc)

            if isinstance(Tw, (int, float)):
                self.obj = <CanteraPFR*> new HeatWallPFR(
                    cpp_mech, cpp_phase, c_Di, c_T0, c_p0, cpp_X0, c_Q0,
                    c_htc, <double>Tw
                )
            else:
                self.cl = self._get_closure(Tw)
                self.obj = <CanteraPFR*> new HeatWallPFR(
                    cpp_mech, cpp_phase, c_Di, c_T0, c_p0, cpp_X0, c_Q0,
                    c_htc, self.cl.get_fun_ptr()
                )
        else:
            raise ValueError(f'Unknown reactor type : {rtype}')

        self.sol = new SolvePFR(self.obj)

    def __dealloc__(self):
        if self.sol != NULL:
            del self.sol
            self.sol = NULL

        if self.obj != NULL:
            del self.obj
            self.obj = NULL

    def _get_closure(self, Tw):
        ftype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)

        if isinstance(Tw, (int, float)):
            return Closure(lambda x: float(Tw), ftype)

        try:
            Tw(0.0)
        except Exception as err:
            raise ValueError(f'Provided function is not suitable : {err}')

        return Closure(Tw, ftype)

    def set_tolerances(self, double rtol, double atol):
        """ Set integration relative and absolute tolerances. """
        self.sol.setTolerances(rtol, atol)

    def set_max_num_steps(self, unsigned int num):
        """ Set maximum number of internal integration steps. """
        self.sol.setMaxNumSteps(num)

    def set_initial_step_size(self, double h0):
        """ Set length of initial integration step. """
        self.sol.setInitialStepSize(h0)

    def set_stop_position(self, double xstop):
        """ Set position to stop integration. """
        self.sol.setStopPosition(xstop)

    def solve(self, double xout) -> int:
        """ Integrate to target position. """
        return self.sol.solve(xout)

    def solution(self, size_t num) -> double:
        """ Return solution component by index. """
        return self.sol.solution(num)

    def solution_vector(self):
        """ Return complete solution vector at current position. """
        cdef vector[double] vec = self.sol.solutionVector()
        return np.array(vec, dtype=np.float64)

    def solutionVector(self):
        """ Alias for solution_vector. """
        return self.solution_vector()

    def variables_names(self):
        """ Return list of model variable names. """
        cdef vector[string] names = self.sol.variablesNames()
        return [
            n.decode('utf-8') if isinstance(n, bytes) else str(n)
            for n in names
        ]

    def variablesNames(self):
        """ Alias for variables_names. """
        return self.variables_names()

    def write_results(self, str saveas):
        """ Write buffered solution results to file. """
        cdef string cpp_saveas = saveas.encode('utf-8')
        self.sol.writeResults(cpp_saveas)

    def manage_solution(
        self,
        double Lr,
        double dx,
        object saveas=None,
        int outfreq=100,
        bint overwrite=False
    ):
        """ Manage full integration along reactor length.

        Parameters
        ----------
        Lr : float
            Total length of reactor in meters.
        dx : float
            Step size to record results in meters.
        saveas : str or Path, optional
            Path to results CSV file.
        outfreq : int, optional
            Screen output frequency.
        overwrite : bool, optional
            Whether to overwrite existing output file.

        Returns
        -------
        DataFrame
            Table of results along reactor axis.
        """
        assert Lr > 0, 'Reactor length must be positive'
        assert dx > 0, 'Saving step must be positive'
        assert Lr > dx, 'Reactor length must be greater than step'
        assert outfreq > 0, 'Output frequency must be positive'

        if saveas is not None and not overwrite:
            if Path(saveas).exists():
                raise FileExistsError(f'Output file already exists: {saveas}')

        columns = ['x'] + self.variables_names()
        sol0 = np.concatenate(([0.0], self.solution_vector()))

        saveat = np.arange(0.0, Lr + dx / 2.0, dx)
        self.sol.setStopPosition(Lr)

        rows = [sol0]
        t0 = time.time()

        for i, x in enumerate(saveat[1:]):
            if not (i % outfreq):
                print(f'Currently solving to reach {x:.4f} m')

            ret = self.sol.solve(x)

            if ret != 0:
                raise RuntimeError(
                    f'Failure during integration at x = {x:.4f} m (code {ret})'
                )

            row = np.concatenate(([x], self.solution_vector()))
            rows.append(row)

        print(f'Integration took {time.time() - t0:.5f} s')
        data = DataFrame(rows, columns=columns)

        if saveas is not None:
            data.to_csv(saveas, index=False)

        return data
