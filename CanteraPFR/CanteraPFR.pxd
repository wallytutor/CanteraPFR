# -*- coding: utf-8 -*-
from libcpp.string cimport string


cdef extern from "cantera/numerics/ResidJacEval.h" namespace "Cantera":
    cdef cppclass ResidJacEval:
        ResidJacEval(double atol = 1.0e-13) except+

    cdef enum ResidEval_Type_Enum:
        Base_ResidEval = 0,
        JacBase_ResidEval,
        JacDelta_ResidEval,
        Base_ShowSolution,
        Base_LaggedSolutionComponents


cdef extern from "CanteraPFR/CanteraPFR.hpp" namespace "Cantera":
    cdef cppclass CanteraPFR(ResidJacEval):
        CanteraPFR(const string& mech, string phase, double T0, double p0,
                   string X0, unsigned neqs_extra_) except+
        unsigned getSpeciesIndex(string name) const
        double getIntEnergyMass() const
        # void setViscosityFunc(std::function<double()> const& mu)
        int getInitialConditions(const double t0, double *const y,
                                 double *const ydot)
        int evalResidNJ(const double t, const double delta_t,
                        const double* const y, const double* const ydot,
                        double* const resid,
                        const ResidEval_Type_Enum evalType = Base_ResidEval,
                        const int id_x = -1, const double delta_x = 0.0)


cdef extern from "CanteraPFR/ConstAreaPFR.hpp" namespace "Cantera":
    cdef cppclass ConstAreaPFR(CanteraPFR):
        CanteraPFR(const string& mech, string phase, double Di, double T0,
                   double p0, string X0, double Q0,
                   unsigned neqs_extra_) except+


cdef extern from "CanteraPFR/AdiabaticPFR.hpp" namespace "Cantera":
    cdef cppclass AdiabaticPFR(ConstAreaPFR):
        AdiabaticPFR(const string& mech, string phase, double Di, double T0,
                     double p0, string X0, double Q0) except+


cdef extern from "CanteraPFR/HeatWallPFR.hpp" namespace "Cantera":
    cdef cppclass HeatWallPFR(ConstAreaPFR):
        HeatWallPFR(const string& mech, string phase, double Di, double T0,
                    double p0, string X0, double Q0, double htc,
                    double Tw) except+


cdef extern from "CanteraPFR/IsothermalPFR.hpp" namespace "Cantera":
    cdef cppclass IsothermalPFR(ConstAreaPFR):
        IsothermalPFR(const string& mech, string phase, double Di,
                      double T0, double p0, string X0, double Q0) except+


cdef extern from "CanteraPFR/SolvePFR.hpp" namespace "Cantera":
    cdef cppclass SolvePFR:
        SolvePFR(CanteraPFR* pfr) except+
        void setTolerances(double reltol, double abstol)
        void setMaxNumSteps(int n)
        void setInitialStepSize(double h0)
        void setStopTime(double tstop)
        int solve(double tout)
        double solution(int k) const
