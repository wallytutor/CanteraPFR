// ***************************************************************************
// Provides an interface to IDA_Solver
//
// Author : Walter Dal'Maz Silva
// Date   : December 30th 2018
// ***************************************************************************

#ifndef __SOLVEPFR_HPP__
#define __SOLVEPFR_HPP__

#include "CanteraPFR/CanteraPFR.hpp"

namespace Cantera {

class SolvePFR
{
public:
    SolvePFR(CanteraPFR* pfr)
    {
        try
        {
            m_solver = new IDA_Solver {*pfr};
            m_solver->init(0.0);
            m_solver->setJacobianType(0);
            m_solver->setDenseLinearSolver();
        }
        catch (CanteraError& err)
        {
            std::cerr << err.what() << std::endl;
        }

    }

    ~SolvePFR()
    {
        if (m_solver != nullptr) delete m_solver;
    }

    // void init(doublereal t0)
    // {
    //     m_solver->init(t0);
    // }

    void setTolerances(doublereal rtol, doublereal atol)
    {
        m_solver->setTolerances(rtol, atol);
    }

    void setMaxNumSteps(unsigned maxsteps)
    {
        m_solver->setMaxNumSteps(maxsteps);
    }

    void setInitialStepSize(doublereal h0)
    {
        m_solver->setInitialStepSize(h0);
    }

    void setStopTime(doublereal tstop)
    {
        m_solver->setStopTime(tstop);
    }

    int solve(doublereal tout)
    {
        return m_solver->solve(tout);
    }

    doublereal solution(unsigned num)
    {
        return m_solver->solution(num);
    }

protected:
    //! Pointer to IDA solver.
    IDA_Solver* m_solver;

}; // (class SolvePFR)

} // (namespace Cantera)

#endif // (__SOLVEPFR_HPP__)

// ***************************************************************************
//                                  EOF
// ***************************************************************************
