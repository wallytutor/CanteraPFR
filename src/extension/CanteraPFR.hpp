// ***************************************************************************
// Provides a generic plug-flow reactor model.
//
// Author : Walter Dal'Maz Silva
// Date   : December 30th 2018
// ***************************************************************************

#ifndef __CANTERAPFR_HPP__
#define __CANTERAPFR_HPP__

#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include "cantera/base/config.h"
#include "cantera/base/global.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/ThermoPhase.h"
#include "cantera/kinetics/Kinetics.h"
#include "cantera/transport/Transport.h"
#include "cantera/numerics/FuncEval.h"
#include "cantera/numerics/IdasIntegrator.h"
#include "cantera/numerics/eigen_dense.h"

#ifndef VISCOSITY_DEFAULT
#define VISCOSITY_DEFAULT 3.957996309582866e-05
#endif

#ifndef SUPPRESS_WARNINGS
#define SUPPRESS_WARNINGS true
#endif

namespace Cantera
{

static inline double circleArea(double Di)
{
    return Cantera::Pi * Di * Di / 4;
}

static inline double sccmTocmps(double sccm)
{
    return sccm / 60000000.0;
}

class CanteraPFR : public FuncEval
{
public:
    CanteraPFR(std::string const& mech,
               std::string phase,
               double T0,
               double p0,
               std::string X0,
               unsigned neqs_extra_)
        : FuncEval{}
    {
        suppress_thermo_warnings(SUPPRESS_WARNINGS);

        std::cout << std::boolalpha
                  << "\nIntegrating PFR"
                  << "\nUsing Cantera  : " << CANTERA_VERSION
                  << std::endl;

        m_sol = newSolution(mech, phase);
        m_gas = m_sol->thermo();
        m_kin = m_sol->kinetics();
        m_trn = m_sol->transport();

        m_gas->setState_TPX(273.15, Cantera::OneAtm, X0);
        m_rho_ref = m_gas->density();

        m_gas->setState_TPX(T0, p0, X0);

        nspec_gas_ = static_cast<unsigned>(m_gas->nSpecies());
        neq_ = nspec_gas_ + neqs_extra_;

        m_W.resize(nspec_gas_);
        m_wdot.resize(nspec_gas_);
        m_gas->getMolecularWeights(m_W.data());

        m_var = m_gas->speciesNames();

        try
        {
            std::cout << "Getting transport manager..." << std::endl;
            if (m_trn) {
                m_mu = [this]() { return m_trn->viscosity(); };
                m_mu();
                std::cout << "Using mechanism viscosity..." << std::endl;
            } else {
                throw CanteraError("CanteraPFR", "No transport manager available.");
            }
        }
        catch (CanteraError &err)
        {
            std::cerr << err.what() << std::endl;
            std::cout << "Fallback to default viscosity..." << std::endl;
            m_mu = [](){ return VISCOSITY_DEFAULT; };
        }
    }

    ~CanteraPFR() override = default;

    virtual int getInitialConditions(const double t0,
                                     double *const y,
                                     double *const ydot) = 0;

    virtual int evalResidNJ(const double t,
                            const double delta_t,
                            const double* const y,
                            const double* const ydot,
                            double* const resid,
                            const int id_x = -1,
                            const double delta_x = 0.0) = 0;

    size_t neq() const override
    {
        return neq_;
    }

    size_t nEquations() const
    {
        return neq_;
    }

    void getStateDae(double* y, double* ydot) override
    {
        getInitialConditions(0.0, y, ydot);
    }

    void evalDae(double t, double* y, double* ydot, double* p, double* residual) override
    {
        evalResidNJ(t, 0.0, y, ydot, residual);
    }

    void getConstraints(double* constraints) override
    {
        std::fill(constraints, constraints + neq_, 1.0);
    }

    unsigned getSpeciesIndex(std::string name) const
    {
        return static_cast<unsigned>(m_kin->kineticsSpeciesIndex(name));
    }

    double getIntEnergyMass() const
    {
        return m_gas->intEnergy_mass();
    }

    void setViscosityFunc(std::function<double()> const& mu)
    {
        m_mu = mu;
    }

    std::vector<std::string> variablesNames() const
    {
        return m_var;
    }

protected:
    //! Managed pointer to Solution object.
    std::shared_ptr<Solution> m_sol;

    //! Pointer to the gas phase thermo object.
    std::shared_ptr<ThermoPhase> m_gas;

    //! Pointer to the kinetics manager.
    std::shared_ptr<Kinetics> m_kin;

    //! Pointer to the transport manager.
    std::shared_ptr<Transport> m_trn;

    //! Viscosity function interface.
    std::function<double()> m_mu;

    //! Species molar weights.
    std::vector<double> m_W;

    //! Species net production rates.
    std::vector<double> m_wdot;

    //! Species names and variables.
    std::vector<std::string> m_var;

    //! Number of gas phase species.
    unsigned nspec_gas_;

    //! Total number of equations.
    unsigned neq_;

    //! Reference state inlet density.
    double m_rho_ref;

}; // (class CanteraPFR)

} // (namespace Cantera)

#endif // (__CANTERAPFR_HPP__)
