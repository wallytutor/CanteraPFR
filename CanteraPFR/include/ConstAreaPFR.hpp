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
                 doublereal Di,
                 doublereal T0,
                 doublereal p0,
                 std::string X0,
                 doublereal Q0,
                 unsigned neqs_extra)
        : CanteraPFR{mech, phase, T0, p0, X0, neqs_extra},
          m_Ac{circleArea(Di)},
          m_u0{setVelocity(Q0)}
        {}

    int getInitialConditions(const doublereal t0,
                             doublereal *const y,
                             doublereal *const ydot) = 0;

    int evalResidNJ(const doublereal t,
                    const doublereal delta_t,
                    const doublereal* const y,
                    const doublereal* const ydot,
                    doublereal* const resid,
                    const ResidEval_Type_Enum evalType = Base_ResidEval,
                    const int id_x = -1,
                    const doublereal delta_x = 0.0) = 0;

    //! Compute inlet velocity.
    virtual const doublereal setVelocity(doublereal Q0) const
    {
        return (m_rho_ref / m_gas->density()) * sccmTocmps(Q0) / m_Ac;
    }

    //! Pressure drop model of viscous loss.
    virtual const doublereal viscousLoss(doublereal u) const
    {
        // TODO Re number test.
        return 8 * m_mu() * u * Cantera::Pi / m_Ac;
    }

protected:
    //! Reactor cross-sectional area.
    const doublereal m_Ac = 0.0;

    //! Inlet velocity.
    const doublereal m_u0 = 0.0;

}; // (class ConstAreaPFR)

} // (namespace Cantera)

#endif // (__CONSTAREAPFR_HPP__)

// ***************************************************************************
//                                  EOF
// ***************************************************************************
