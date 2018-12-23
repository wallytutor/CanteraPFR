// ***************************************************************************
// Provides an isothermal plug-flow reactor model.
//
// Author : Walter Dal'Maz Silva
// Date   : December 21st 2018
//
// TODO provide analytical Jacobian.
// ***************************************************************************

#ifndef __IsothermalPFR_HPP__
#define __IsothermalPFR_HPP__

#include <functional>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "cantera/IdealGasMix.h"
#include "cantera/numerics/IDA_Solver.h"
#include "cantera/numerics/eigen_dense.h"
#include "cantera/transport.h"

#ifndef VISCOSITY_DEFAULT
#define VISCOSITY_DEFAULT 3.957996309582866e-05
#endif

namespace Cantera
{

class IsothermalPFR : public ResidJacEval
{
public:
    IsothermalPFR(std::string const& mech,
                  std::string phase,
                  doublereal Di,
                  doublereal T0,
                  doublereal p0,
                  std::string X0,
                  doublereal Q0)
        : ResidJacEval{}
    {
        m_gas = new IdealGasMix {mech, phase};
        getTransportManager();

        m_gas->setState_TPX(273.15, Cantera::OneAtm, X0);
        doublereal rhor = m_gas->density();

        m_gas->setState_TPX(T0, p0, X0);
        doublereal rho0 = m_gas->density();

        m_T0 = T0;
        m_Ac = Cantera::Pi * Di * Di / 4;
        m_u0 = Q0 * rhor / (60000000 * rho0 * m_Ac);

        neq_ = m_gas->nSpecies() + 3;
        m_W.resize(neq_-3);
        m_wdot.resize(neq_-3);
        m_gas->getMolecularWeights(m_W);

        std::cout << "\nInitial temperature (K) . " << m_gas->temperature()
                  << "\nInitial pressure (Pa) ... " << m_gas->pressure()
                  << "\nInitial velocity (m/s) .. " << m_u0
                  << "\nNumber of equations ..... " << neq_
                  << std::endl;
    }

    ~IsothermalPFR()
    {
        if (m_gas != nullptr) delete m_gas;
        if (m_trn != nullptr) delete m_trn;
    }

    int getInitialConditions(const doublereal t0,
                             doublereal *const y,
                             doublereal *const ydot);

    int evalResidNJ(const doublereal t,
                    const doublereal delta_t,
                    const doublereal* const y,
                    const doublereal* const ydot,
                    doublereal* const resid,
                    const ResidEval_Type_Enum evalType = Base_ResidEval,
                    const int id_x = -1,
                    const doublereal delta_x = 0.0);

    size_t getSpeciesIndex(std::string name)
    {
        return m_gas->kineticsSpeciesIndex(name);
    }

    void setViscosityFunc(std::function<doublereal()>& mu)
    {
        m_mu = mu;
    }

private:
    //! Pointer to the gas phase object.
    IdealGasMix *m_gas;

    //! Pointer to the transport manager.
    Transport* m_trn = nullptr;

    //! Viscosity function interface.
    std::function<doublereal()> m_mu;

    //! Reactor cross-sectional area.
    doublereal m_Ac;

    //! Reactor temperature.
    doublereal m_T0;

    //! Inlet velocity.
    doublereal m_u0;

    //! Species molar weights.
    std::vector<doublereal> m_W;

    //! Species net production rates.
    std::vector<doublereal> m_wdot;

    //! Try to retrieve viscosity from mechanism.
    void getTransportManager()
    {
        try {
            std::cout << "Getting transport manager..." << std::endl;
            m_trn = newDefaultTransportMgr(m_gas);
            m_mu = [&](){ return m_trn->viscosity(); };
            m_mu();
        } catch (CanteraError &err) {
            std::cerr << err.what() << std::endl;
            std::cout << "Fallback to default viscosity..." << std::endl;
            m_mu = viscosity_default;
        }
    }

    //! Default viscosity function.
    static doublereal viscosity_default()
    {
        return VISCOSITY_DEFAULT;
    }
}; // (class IsothermalPFR)

} // (namespace Cantera)

#endif // (__IsothermalPFR_HPP__)

// ***************************************************************************
//                                  EOF
// ***************************************************************************
