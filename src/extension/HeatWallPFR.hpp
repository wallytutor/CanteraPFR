// ***************************************************************************
// Provides a plug-flow reactor model with wall exchanges.
//
// Author : Walter Dal'Maz Silva
// Date   : December 31st 2018
// ***************************************************************************

#ifndef __HEATWALLPFR_HPP__
#define __HEATWALLPFR_HPP__

#include "CanteraPFR/ConstAreaPFR.hpp"

namespace Cantera
{

class HeatWallPFR : public ConstAreaPFR
{
public:
    HeatWallPFR(std::string const& mech,
                 std::string phase,
                 doublereal Di,
                 doublereal T0,
                 doublereal p0,
                 std::string X0,
                 doublereal Q0,
                 doublereal htc,
                 std::function<doublereal(doublereal)> const& Tw)
        : ConstAreaPFR{mech, phase, Di, T0, p0, X0, Q0, neqs_extra_},
          idx0{nspec_gas_+0},
          idx1{nspec_gas_+1},
          idx2{nspec_gas_+2},
          idx3{nspec_gas_+3},
          m_PoA{4.0 / Di},
          m_htc{htc},
          Tw{Tw}
    {
        m_hbar.resize(nspec_gas_);
        m_var.push_back("u");
        m_var.push_back("rho");
        m_var.push_back("p");
        m_var.push_back("T");

        std::cout << "\nStarting solver : " << "HeatWallPFR"
                  << "\nInitial temperature (K) . " << m_gas->temperature()
                  << "\nInitial pressure (Pa) ... " << m_gas->pressure()
                  << "\nInitial velocity (m/s) .. " << m_u0
                  << "\nNumber of equations ..... " << neq_
                  << std::endl;
    }

    HeatWallPFR(std::string const& mech,
                 std::string phase,
                 doublereal Di,
                 doublereal T0,
                 doublereal p0,
                 std::string X0,
                 doublereal Q0,
                 doublereal htc,
                 doublereal Tw)
        : HeatWallPFR{mech, phase, Di, T0, p0, X0, Q0, htc,
                      [&](doublereal x){ return Tw; }} {}

    // TODO provide struct interface.
    // HeatWallPFR(const ReactorConfig & r) {}

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
    static constexpr unsigned neqs_extra_ = 4;

    //! Index of extra equations.
    const unsigned idx0 = 0, idx1 = 0, idx2 = 0, idx3 = 0;

    //! Species molar enthalpies.
    std::vector<doublereal> m_hbar;

    //! Ratio perimeter per cross section.
    const doublereal m_PoA = 0.0;

    //! Global heat transfer coefficient.
    const doublereal m_htc = 0.0;

    //! Wall temperature in terms of position.
    std::function<doublereal(doublereal)> Tw;

    //! Wall heat exchange term.
    const doublereal wallHeatExchange(doublereal x, doublereal T) const
    {
        return m_htc * m_PoA * (Tw(x) - T);
    }

}; // (class HeatWallPFR)

} // (namespace Cantera)

#endif // (__HEATWALLPFR_HPP__)

// ***************************************************************************
//                                  EOF
// ***************************************************************************
