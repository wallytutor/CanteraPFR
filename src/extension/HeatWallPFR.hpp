// ***************************************************************************
// Provides a plug-flow reactor model with wall exchanges.
//
// Author : Walter Dal'Maz Silva
// Date   : December 31st 2018
// ***************************************************************************

#ifndef __HEATWALLPFR_HPP__
#define __HEATWALLPFR_HPP__

#include "ConstAreaPFR.hpp"

namespace Cantera
{

class HeatWallPFR : public ConstAreaPFR
{
public:
    HeatWallPFR(std::string const& mech,
                std::string phase,
                double Di,
                double T0,
                double p0,
                std::string X0,
                double Q0,
                double htc,
                std::function<double(double)> const& Tw)
        : ConstAreaPFR{mech, phase, Di, T0, p0, X0, Q0, neqs_extra_},
          idx0{nspec_gas_ + 0},
          idx1{nspec_gas_ + 1},
          idx2{nspec_gas_ + 2},
          idx3{nspec_gas_ + 3},
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
                double Di,
                double T0,
                double p0,
                std::string X0,
                double Q0,
                double htc,
                double Tw_const)
        : HeatWallPFR{mech, phase, Di, T0, p0, X0, Q0, htc,
                      [Tw_const](double){ return Tw_const; }} {}

    int getInitialConditions(const double t0,
                             double *const y,
                             double *const ydot) override;

    int evalResidNJ(const double t,
                    const double delta_t,
                    const double* const y,
                    const double* const ydot,
                    double* const resid,
                    const int id_x = -1,
                    const double delta_x = 0.0) override;

private:
    //! Number of extra equations.
    static constexpr unsigned neqs_extra_ = 4;

    //! Index of extra equations.
    const unsigned idx0 = 0, idx1 = 0, idx2 = 0, idx3 = 0;

    //! Species molar enthalpies.
    std::vector<double> m_hbar;

    //! Ratio perimeter per cross section.
    const double m_PoA = 0.0;

    //! Global heat transfer coefficient.
    const double m_htc = 0.0;

    //! Wall temperature in terms of position.
    std::function<double(double)> Tw;

    //! Wall heat exchange term.
    double wallHeatExchange(double x, double T) const
    {
        return m_htc * m_PoA * (Tw(x) - T);
    }

}; // (class HeatWallPFR)

} // (namespace Cantera)

#endif // (__HEATWALLPFR_HPP__)
