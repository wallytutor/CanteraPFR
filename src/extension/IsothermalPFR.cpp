// ***************************************************************************
// Provides an isothermal plug-flow reactor model.
//
// Author : Walter Dal'Maz Silva
// Date   : December 21st 2018
// ***************************************************************************

#include "CanteraPFR/IsothermalPFR.hpp"

int Cantera::IsothermalPFR::getInitialConditions(const doublereal t0,
    doublereal *const y, doublereal *const ydot)
{
    const doublereal P0 = m_gas->pressure();
    const doublereal rho0 = m_gas->density();
    const doublereal Wavg = m_gas->meanMolecularWeight();
    const doublereal RT = m_gas->temperature() * Cantera::GasConstant;

    m_gas->getMassFractions(y);
    m_gas->getNetProductionRates(&m_wdot[0]);

    Eigen::MatrixXd A(neq_, neq_);
    Eigen::VectorXd b(neq_);

    y[idx0] = m_u0;
    y[idx1] = rho0;
    y[idx2] = P0;

    for (unsigned i = 0; i != idx0; ++i)
    {
        // For species equations.
        A(i, i) = rho0 * m_u0;
        b(i) = m_wdot[i] * m_W[i];

        // Yk' for other equations, exceptionally here!
        A(idx2, i) = P0 * Wavg * Wavg / m_W[i];
    }

    // Continuity equation elements.
    A(idx0, idx0) = rho0;           // u'
    A(idx0, idx1) = m_u0;           // rho'
    A(idx0, idx2) = 0;              // p'

    // Momentum equation elements.
    A(idx1, idx0) = rho0 * m_u0;    // u'
    A(idx1, idx1) = 0;              // rho'
    A(idx1, idx2) = 1;              // p'

    // State equation elements.
    A(idx2, idx0) = 0;              // u'
    A(idx2, idx1) = RT;             // rho'
    A(idx2, idx2) = -Wavg;          // p'

    b(idx0) = 0;                    // RHS continuity
    b(idx1) = -viscousLoss(m_u0);   // RHS momentum
    b(idx2) = 0;                    // RHS state

    Eigen::VectorXd x = A.fullPivLu().solve(b);
    Eigen::VectorXd::Map(ydot, x.rows()) = x;

    return 0;
}

int Cantera::IsothermalPFR::evalResidNJ(const doublereal t,
    const doublereal delta_t, const doublereal* const y,
    const doublereal* const ydot, doublereal* const resid,
    const Cantera::ResidEval_Type_Enum evalType, const int id_x,
    const doublereal delta_x)
{
    const doublereal u = y[idx0];
    const doublereal r = y[idx1];
    const doublereal p = y[idx2];

    const doublereal dudz = ydot[idx0];
    const doublereal drdz = ydot[idx1];
    const doublereal dpdz = ydot[idx2];

    m_gas->setMassFractions_NoNorm(y);
    m_gas->setState_TP(m_T0, p);
    m_gas->getNetProductionRates(&m_wdot[0]);

    for (unsigned k = 0; k != idx0; ++k)
    {
        resid[k] = u * r * ydot[k] - m_wdot[k] * m_W[k];
    }

    resid[idx0] = r * dudz + u * drdz;
    resid[idx1] = u * r * dudz + dpdz + viscousLoss(u);
    resid[idx2] = m_gas->density() - r;

    return 0;
}

// ***************************************************************************
//                                  EOF
// ***************************************************************************
