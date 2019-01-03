# -*- coding: utf-8 -*-
from libcpp.string cimport string
from CanteraPFR cimport AdiabaticPFR
from CanteraPFR cimport HeatWallPFR
from CanteraPFR cimport IsothermalPFR
from CanteraPFR cimport SolvePFR

import time
from numpy import arange
from numpy import array
from pandas import DataFrame


cdef class PyPFR:
    cdef CanteraPFR* obj
    cdef SolvePFR* sol

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
            self.obj = new HeatWallPFR(mech, phase, Di, T0, p0, X0, Q0, htc, Tw)
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
