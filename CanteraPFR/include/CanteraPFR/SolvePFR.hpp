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
        m_vec.resize(pfr->nEquations());
        m_var = pfr->variablesNames();

        try
        {
            m_solver = new IDA_Solver {*pfr};
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

    void init(doublereal t0)
    {
        m_solver->init(t0);
    }

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

    void setStopPosition(doublereal tstop)
    {
        m_solver->setStopTime(tstop);
    }

    int solve(doublereal tout)
    {
        // TODO Manage return codes from IDA_Solver.solve.
        try
        {
            return m_solver->solve(tout);
        }
        catch (CanteraError& err)
        {
            std::cerr << err.what() << std::endl;
            return -99;
        }
    }

    doublereal solution(unsigned num) const
    {
        return m_solver->solution(num);
    }

    std::vector<doublereal> solutionVector()
    {
        // TODO make this with STL algorithm.
        const doublereal* sol = m_solver->solutionVector();
        for (unsigned i = 0; i != m_vec.size(); ++i) { m_vec[i] = sol[i]; }
        return m_vec;
    }

    std::vector<std::string> variablesNames() const
    {
        return m_var;
    }

protected:
    //! To provide access to results.
    std::vector<doublereal> m_vec;

    //! Provides access to variables names.
    std::vector<std::string> m_var;

    //! Pointer to IDA solver.
    IDA_Solver* m_solver;

}; // (class SolvePFR)

} // (namespace Cantera)

#endif // (__SOLVEPFR_HPP__)

// ***************************************************************************
//                                  EOF
// ***************************************************************************
