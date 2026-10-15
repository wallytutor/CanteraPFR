// ***************************************************************************
// Provides an interface to IdasIntegrator
//
// Author : Walter Dal'Maz Silva
// Date   : December 30th 2018
// ***************************************************************************

#ifndef __SOLVEPFR_HPP__
#define __SOLVEPFR_HPP__

#include <iostream>
#include <fstream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "CanteraPFR.hpp"
#include "cantera/numerics/IdasIntegrator.h"

namespace Cantera {

class SolvePFR
{
public:
    SolvePFR(CanteraPFR* pfr)
    {
        m_neq = pfr->nEquations();
        m_vec.resize(m_neq);
        m_var = pfr->variablesNames();

        try
        {
            m_solver = std::make_unique<IdasIntegrator>();
            m_solver->setLinearSolverType("DENSE");
            m_solver->initialize(0.0, *pfr);
        }
        catch (CanteraError& err)
        {
            std::cerr << err.what() << std::endl;
        }
    }

    ~SolvePFR() = default;

    void setTolerances(double rtol, double atol)
    {
        if (m_solver) {
            m_solver->setTolerances(rtol, atol);
        }
    }

    void setMaxNumSteps(unsigned maxsteps)
    {
        if (m_solver) {
            m_solver->setMaxSteps(static_cast<int>(maxsteps));
        }
    }

    void setInitialStepSize(double h0)
    {
        // Initial step size is managed by IdasIntegrator
    }

    void setStopPosition(double tstop)
    {
        // Stop position is handled at integration target
    }

    int solve(double xout)
    {
        int retcode = 0;

        if (!m_ss_started && m_solver)
        {
            for (const auto& var : m_var) { m_ss << var << ","; }
            m_ss << "x\n";

            for (size_t i = 0; i != m_neq; ++i) {
                m_ss << m_solver->solution(i) << ",";
            }
            m_ss << 0.0 << "\n";

            m_ss_started = true;
        }

        try
        {
            if (m_solver) {
                m_solver->integrate(xout);
            } else {
                return -1;
            }
        }
        catch (CanteraError& err)
        {
            std::cerr << err.what() << std::endl;
            retcode = -99;
        }

        for (size_t i = 0; i != m_neq; ++i) {
            m_ss << m_solver->solution(i) << ",";
        }
        m_ss << xout << "\n";

        return retcode;
    }

    double solution(size_t num) const
    {
        return m_solver ? m_solver->solution(num) : 0.0;
    }

    std::vector<double> solutionVector()
    {
        if (m_solver) {
            const double* sol = m_solver->solution();
            for (size_t i = 0; i != m_vec.size(); ++i) {
                m_vec[i] = sol[i];
            }
        }
        return m_vec;
    }

    std::vector<std::string> variablesNames() const
    {
        return m_var;
    }

    void writeResults(const std::string & saveas)
    {
        std::ofstream ofs(saveas, std::ios::out | std::ios::binary);
        if (!ofs)
        {
            throw std::runtime_error("Cannot output to file");
        }
        ofs << m_ss.str();
        ofs.close();
    }

protected:
    //! To provide access to results.
    std::vector<double> m_vec;

    //! Provides access to variables names.
    std::vector<std::string> m_var;

    //! Pointer to IDA solver.
    std::unique_ptr<IdasIntegrator> m_solver;

    //! Buffer for results file.
    std::stringstream m_ss;

    //! Check if writer was started.
    bool m_ss_started = false;

    //! Number of equations in reactor model.
    size_t m_neq = 0;

}; // (class SolvePFR)

} // (namespace Cantera)

#endif // (__SOLVEPFR_HPP__)
