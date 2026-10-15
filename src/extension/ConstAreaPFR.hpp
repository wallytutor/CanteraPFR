// ***************************************************************************
// Provides a generic constant area plug-flow reactor model.
//
// Author : Walter Dal'Maz Silva
// Date   : December 30th 2018
// ***************************************************************************

#ifndef __CONSTAREAPFR_HPP__
#define __CONSTAREAPFR_HPP__

#include "CanteraPFR.hpp"

namespace Cantera
{

class ConstAreaPFR : public CanteraPFR
{
public:
    ConstAreaPFR(std::string const& mech,
                 std::string phase,
                 double Di,
                 double T0,
                 double p0,
                 std::string X0,
                 double Q0,
                 unsigned neqs_extra_)
        : CanteraPFR{mech, phase, T0, p0, X0, neqs_extra_},
          m_Ac{circleArea(Di)},
          m_u0{setVelocity(Q0)}
        {}

    virtual int getInitialConditions(const double t0,
                                     double *const y,
                                     double *const ydot) override = 0;

    virtual int evalResidNJ(const double t,
                            const double delta_t,
                            const double* const y,
                            const double* const ydot,
                            double* const resid,
                            const int id_x = -1,
                            const double delta_x = 0.0) override = 0;

    //! Compute inlet velocity.
    virtual double setVelocity(double Q0) const
    {
        return (m_rho_ref / m_gas->density()) * sccmTocmps(Q0) / m_Ac;
    }

    //! Pressure drop model of viscous loss.
    virtual double viscousLoss(double u) const
    {
        return 8.0 * m_mu() * u * Cantera::Pi / m_Ac;
    }

protected:
    //! Reactor cross-sectional area.
    const double m_Ac = 0.0;

    //! Inlet velocity.
    const double m_u0 = 0.0;

}; // (class ConstAreaPFR)

} // (namespace Cantera)

#endif // (__CONSTAREAPFR_HPP__)
