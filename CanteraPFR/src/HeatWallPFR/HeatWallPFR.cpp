// ***************************************************************************
// Provides a plug-flow reactor model with wall exchanges.
//
// Author : Walter Dal'Maz Silva
// Date   : December 31st 2018
// ***************************************************************************

#include "HeatWallPFR.hpp"

int Cantera::HeatWallPFR::getInitialConditions(const doublereal t0,
    doublereal *const y, doublereal *const ydot)
{
    const doublereal T0 = m_gas->temperature();
    const doublereal P0 = m_gas->pressure();
    const doublereal rho0 = m_gas->density();
    const doublereal Wavg = m_gas->meanMolecularWeight();
    const doublereal RT = T0 * Cantera::GasConstant;
    const doublereal rho0R = rho0 * Cantera::GasConstant;
    const doublereal rhoUCp = rho0 * m_u0 * m_gas->cp_mass();
    doublereal hdot = 0;
    doublereal wall = 0 ; //wallHeatExchange(t0, T0);

    m_gas->getMassFractions(y);
    m_gas->getNetProductionRates(&m_wdot[0]);
    m_gas->getPartialMolarEnthalpies(&m_hbar[0]);

    Eigen::MatrixXd A(neq_, neq_);
    Eigen::VectorXd b(neq_);

    y[idx0] = m_u0;
    y[idx1] = rho0;
    y[idx2] = P0;
    y[idx3] = T0;

    for (unsigned k = 0; k != idx0; ++k)
    {
        // Energy contribution.
        hdot += m_wdot[k] * m_hbar[k];

        // For species equations.
        A(k, k) = rho0 * m_u0;
        b(k) = m_wdot[k] * m_W[k];

        // Yk' for other equations, exceptionally here!
        A(idx2, k) = P0 * Wavg * Wavg / m_W[k];
    }

    // Continuity equation elements.
    A(idx0, idx0) = rho0;           // u'
    A(idx0, idx1) = m_u0;           // rho'
    A(idx0, idx2) = 0;              // p'
    A(idx0, idx3) = 0;              // T'

    // Momentum equation elements.
    A(idx1, idx0) = rho0 * m_u0;    // u'
    A(idx1, idx1) = 0;              // rho'
    A(idx1, idx2) = 1;              // p'
    A(idx1, idx3) = 0;              // T'

    // State equation elements.
    A(idx2, idx0) = 0;              // u'
    A(idx2, idx1) = RT;             // rho'
    A(idx2, idx2) = -Wavg;          // p'
    A(idx2, idx3) = rho0R;          // T'

    // Energy equation elements.
    A(idx3, idx0) = 0;              // u'
    A(idx3, idx1) = 0;              // rho'
    A(idx3, idx2) = 0;              // p'
    A(idx3, idx3) = rhoUCp;         // T'

    b(idx0) = 0;                    // RHS continuity
    b(idx1) = -viscousLoss(m_u0);   // RHS momentum
    b(idx2) = 0;                    // RHS state
    b(idx3) = -hdot + wall;         // RHS energy

    Eigen::VectorXd x = A.fullPivLu().solve(b);
    Eigen::VectorXd::Map(ydot, x.rows()) = x;

    return 0;
}

int Cantera::HeatWallPFR::evalResidNJ(const doublereal t,
    const doublereal delta_t, const doublereal* const y,
    const doublereal* const ydot, doublereal* const resid,
    const Cantera::ResidEval_Type_Enum evalType, const int id_x,
    const doublereal delta_x)
{
    const doublereal u = y[idx0];
    const doublereal r = y[idx1];
    const doublereal p = y[idx2];
    const doublereal T = y[idx3];

    const doublereal dudz = ydot[idx0];
    const doublereal drdz = ydot[idx1];
    const doublereal dpdz = ydot[idx2];
    const doublereal dTdz = ydot[idx3];

    const doublereal Cp = m_gas->cp_mass();
    doublereal hdot = 0.0;

    m_gas->setMassFractions_NoNorm(y);
    m_gas->setState_TP(T, p);
    m_gas->getNetProductionRates(&m_wdot[0]);
    m_gas->getPartialMolarEnthalpies(&m_hbar[0]);

    for (unsigned k = 0; k != idx0; ++k)
    {
        resid[k] = u * r * ydot[k] - m_wdot[k] * m_W[k];
        hdot += m_wdot[k] * m_hbar[k];
    }

    resid[idx0] = r * dudz + u * drdz;
    resid[idx1] = u * r * dudz + dpdz + viscousLoss(u);
    resid[idx2] = m_gas->density() - r;
    resid[idx3] = r * u * Cp * dTdz + hdot - wallHeatExchange(t, T);

    return 0;
}

// ***************************************************************************
//                                  EOF
// ***************************************************************************
