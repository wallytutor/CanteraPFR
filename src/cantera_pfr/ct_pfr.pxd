# -*- coding: utf-8 -*-
from libcpp.string cimport string
from libcpp.vector cimport vector


ctypedef double (*func_t_none)()
ctypedef double (*func_t)(double)


cdef extern from "cantera/numerics/FuncEval.h" namespace "Cantera":
    cdef cppclass FuncEval:
        pass


cdef extern from "CanteraPFR.hpp" namespace "Cantera":
    cdef cppclass CanteraPFR(FuncEval):
        CanteraPFR(const string& mech, string phase, double T0, double p0,
                   string X0, unsigned neqs_extra_) except+
        unsigned getSpeciesIndex(string name) const
        double getIntEnergyMass() const
        void setViscosityFunc(func_t_none mu)
        size_t nEquations() const
        vector[string] variablesNames() const


cdef extern from "ConstAreaPFR.hpp" namespace "Cantera":
    cdef cppclass ConstAreaPFR(CanteraPFR):
        ConstAreaPFR(const string& mech, string phase, double Di, double T0,
                     double p0, string X0, double Q0,
                     unsigned neqs_extra_) except+


cdef extern from "AdiabaticPFR.hpp" namespace "Cantera":
    cdef cppclass AdiabaticPFR(ConstAreaPFR):
        AdiabaticPFR(const string& mech, string phase, double Di, double T0,
                     double p0, string X0, double Q0) except+


cdef extern from "HeatWallPFR.hpp" namespace "Cantera":
    cdef cppclass HeatWallPFR(ConstAreaPFR):
        HeatWallPFR(const string& mech, string phase, double Di, double T0,
                    double p0, string X0, double Q0, double htc,
                    func_t Tw) except+
        HeatWallPFR(const string& mech, string phase, double Di, double T0,
                    double p0, string X0, double Q0, double htc,
                    double Tw_const) except+


cdef extern from "IsothermalPFR.hpp" namespace "Cantera":
    cdef cppclass IsothermalPFR(ConstAreaPFR):
        IsothermalPFR(const string& mech, string phase, double Di,
                      double T0, double p0, string X0, double Q0) except+


cdef extern from "SolvePFR.hpp" namespace "Cantera":
    cdef cppclass SolvePFR:
        SolvePFR(CanteraPFR* pfr) except+
        void setTolerances(double reltol, double abstol)
        void setMaxNumSteps(unsigned n)
        void setInitialStepSize(double h0)
        void setStopPosition(double tstop)
        int solve(double tout)
        double solution(size_t k) const
        vector[double] solutionVector()
        vector[string] variablesNames() const
        void writeResults(const string& saveas)
