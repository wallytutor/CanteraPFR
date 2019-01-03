# -*- coding: utf-8 -*-
from libcpp.string cimport string
from CanteraPFR cimport ResidJacEval
from CanteraPFR cimport AdiabaticPFR
from CanteraPFR cimport HeatWallPFR
from CanteraPFR cimport IsothermalPFR
from CanteraPFR cimport SolvePFR


cdef class PyIsothermalPFR:
    cdef CanteraPFR* obj
    cdef SolvePFR* sol

    def __cinit__(self, mech, phase, Di, T0, p0, X0, Q0):
        # TODO deal with Python strings!
        # cdef bytearray cpp_mech = bytearray(mech, 'utf8')
        # cdef string cpp_mech = <string> mech.encode('utf-8')
        self.obj = new IsothermalPFR(mech, phase, Di, T0, p0, X0, Q0)
        self.sol = new SolvePFR(<CanteraPFR *> self.obj)
        print('Solver ready')

    def __dealloc__(self):
        if type(self) is PyIsothermalPFR:
            del self.obj
            del self.sol

    def setTolerances(self, rtol, atol):
        self.sol.setTolerances(rtol, atol)

    def setMaxNumSteps(self, num):
        self.sol.setMaxNumSteps(num)

    def setInitialStepSize(self, h0):
        self.sol.setInitialStepSize(h0)

    def setStopTime(self, tstop):
        self.sol.setStopTime(tstop)

    def solve(self, tout):
        return self.sol.solve(tout)

    def solution(self, num):
        return self.sol.solution(num)
