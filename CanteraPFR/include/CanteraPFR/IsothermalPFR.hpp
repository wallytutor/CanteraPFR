// ***************************************************************************
// Provides an isothermal plug-flow reactor model.
//
// Author : Walter Dal'Maz Silva
// Date   : December 21st 2018
// ***************************************************************************

#ifndef __ISOTHERMALPFR_HPP__
#define __ISOTHERMALPFR_HPP__

#include "CanteraPFR/ConstAreaPFR.hpp"

namespace Cantera
{

class IsothermalPFR : public ConstAreaPFR
{
public:
    IsothermalPFR(std::string const& mech,
                  std::string phase,
                  doublereal Di,
                  doublereal T0,
                  doublereal p0,
                  std::string X0,
                  doublereal Q0)
        : ConstAreaPFR{mech, phase, Di, T0, p0, X0, Q0, neqs_extra_},
          idx0{nspec_gas_+0},
          idx1{nspec_gas_+1},
          idx2{nspec_gas_+2},
          m_T0{T0}
    {
        m_var.push_back("u");
        m_var.push_back("rho");
        m_var.push_back("p");

        std::cout << "\nStarting solver : " << "IsothermalPFR"
                  << "\nInitial temperature (K) . " << m_gas->temperature()
                  << "\nInitial pressure (Pa) ... " << m_gas->pressure()
                  << "\nInitial velocity (m/s) .. " << m_u0
                  << "\nNumber of equations ..... " << neq_
                  << std::endl;
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

private:
    //! Number of extra equations.
    static constexpr unsigned neqs_extra_ = 3;

    //! Index of extra equations.
    const unsigned idx0 = 0, idx1 = 0, idx2 = 0;

    //! Reactor temperature.
    const doublereal m_T0 = 0.0;

}; // (class IsothermalPFR)

} // (namespace Cantera)

#endif // (__ISOTHERMALPFR_HPP__)

// ***************************************************************************
//                                  EOF
// ***************************************************************************
