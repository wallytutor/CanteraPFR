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

static inline doublereal circleArea(doublereal Di)
{
    return Cantera::Pi * Di * Di / 4;
}

static inline doublereal sccmTocmps(doublereal sccm)
{
    return sccm / 60000000;
}

class CanteraPFR : public ResidJacEval
{
public:
    CanteraPFR(std::string const& mech,
               std::string phase,
               doublereal T0,
               doublereal p0,
               std::string X0,
               unsigned neqs_extra_)
        : ResidJacEval{}
    {
        m_gas = new IdealGasMix {mech, phase};
        m_gas->setState_TPX(273.15, Cantera::OneAtm, X0);
        m_rho_ref = m_gas->density();

        m_gas->setState_TPX(T0, p0, X0);

        nspec_gas_ = m_gas->nSpecies();
        neq_ = nspec_gas_ + neqs_extra_;

        m_W.resize(nspec_gas_);
        m_wdot.resize(nspec_gas_);
        m_gas->getMolecularWeights(m_W);

        try
        {
            std::cout << "Getting transport manager..." << std::endl;
            m_trn = newDefaultTransportMgr(m_gas);
            m_mu = [&](){ return m_trn->viscosity(); };
            m_mu();
            std::cout << "Using mechanism viscosity..." << std::endl;
        }
        catch (CanteraError &err)
        {
            std::cerr << err.what() << std::endl;
            std::cout << "Fallback to default viscosity..." << std::endl;
            m_mu = [&](){ return VISCOSITY_DEFAULT; };
        }
    }

    ~CanteraPFR()
    {
        if (m_gas != nullptr) delete m_gas;
        if (m_trn != nullptr) delete m_trn;
    }

    virtual int getInitialConditions(const doublereal t0,
                                     doublereal *const y,
                                     doublereal *const ydot) = 0;

    virtual int evalResidNJ(const doublereal t,
                    const doublereal delta_t,
                    const doublereal* const y,
                    const doublereal* const ydot,
                    doublereal* const resid,
                    const ResidEval_Type_Enum evalType = Base_ResidEval,
                    const int id_x = -1,
                    const doublereal delta_x = 0.0) = 0;

    unsigned getSpeciesIndex(std::string name) const
    {
        return m_gas->kineticsSpeciesIndex(name);
    }

    doublereal getIntEnergyMass() const
    {
        return m_gas->intEnergy_mass();
    }

    void setViscosityFunc(std::function<doublereal()> const& mu)
    {
        m_mu = mu;
    }

protected:
  //! Pointer to the gas phase object.
  IdealGasMix *m_gas = nullptr;

  //! Pointer to the transport manager.
  Transport *m_trn = nullptr;

  //! Viscosity function interface.
  std::function<doublereal()> m_mu;

  //! Species molar weights.
  std::vector<doublereal> m_W;

  //! Species net production rates.
  std::vector<doublereal> m_wdot;

  //! Number of gas phase species.
  unsigned nspec_gas_;

  //! Reference state inlet density.
  doublereal m_rho_ref;

}; // (class CanteraPFR)

} // (namespace Cantera)

#endif // (__CANTERAPFR_HPP__)

// ***************************************************************************
//                                  EOF
// ***************************************************************************
