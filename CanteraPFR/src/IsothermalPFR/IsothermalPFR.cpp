// ***************************************************************************
// Provides an isothermal plug-flow reactor model.
//
// Author : Walter Dal'Maz Silva
// Date   : December 21st 2018
// ***************************************************************************

#include "IsothermalPFR.hpp"

int Cantera::IsothermalPFR::getInitialConditions(const doublereal t0,
    doublereal *const y, doublereal *const ydot)
{
    doublereal P0 = m_gas->pressure();
    doublereal rho0 = m_gas->density();
    doublereal Wavg = m_gas->meanMolecularWeight();
    doublereal RT = m_gas->temperature() * Cantera::GasConstant;
    doublereal loss = 8 * m_mu() * m_u0 * Cantera::Pi / m_Ac;

    m_gas->getNetProductionRates(&m_wdot[0]);
    m_gas->getMassFractions(y);

    Eigen::MatrixXd A(neq_, neq_);
    Eigen::VectorXd b(neq_);

    y[neq_-3] = m_u0;
    y[neq_-2] = rho0;
    y[neq_-1] = P0;

    for (int i = 0; i != neq_-3; ++i)
    {
        // For species equations.
        A(i, i) = rho0 * m_u0;
        b(i) = m_wdot[i] * m_W[i];

        // Yk' for other equations, exceptionally here!
        A(neq_-1, i) = P0 * Wavg * Wavg / m_W[i];
    }

    // Continuity equation elements.
    A(neq_-3, neq_-3) = rho0;           // u'
    A(neq_-3, neq_-2) = m_u0;           // rho'
    A(neq_-3, neq_-1) = 0;              // p'

    // Momentum equation elements.
    A(neq_-2, neq_-3) = rho0 * m_u0;    // u'
    A(neq_-2, neq_-2) = 0;              // rho'
    A(neq_-2, neq_-1) = 1;              // p'

    // State equation elements.
    A(neq_-1, neq_-3) = 0;              // u'
    A(neq_-1, neq_-2) = RT;             // rho'
    A(neq_-1, neq_-1) = -Wavg;          // p'

    b(neq_-3) = 0;                      // RHS continuity
    b(neq_-2) = -loss;                  // RHS momentum
    b(neq_-1) = 0;                      // RHS state

    Eigen::VectorXd x = A.fullPivLu().solve(b);

    for (int i = 0; i != neq_; ++i)
    {
        ydot[i] = x(i);
    }

    return 0;
}

int Cantera::IsothermalPFR::evalResidNJ(const doublereal t,
    const doublereal delta_t, const doublereal* const y,
    const doublereal* const ydot, doublereal* const resid,
    const Cantera::ResidEval_Type_Enum evalType, const int id_x,
    const doublereal delta_x)
{
    doublereal u = y[neq_-3];
    doublereal r = y[neq_-2];
    doublereal p = y[neq_-1];

    doublereal dudz = ydot[neq_-3];
    doublereal drdz = ydot[neq_-2];
    doublereal dpdz = ydot[neq_-1];

    m_gas->setMassFractions_NoNorm(y);
    m_gas->setState_TP(m_T0, p);
    m_gas->getNetProductionRates(&m_wdot[0]);

    doublereal loss = 8 * m_mu() * u * Cantera::Pi / m_Ac;

    for (int k = 0; k != neq_-3; ++k)
    {
        resid[k] = u * r * ydot[k] - m_wdot[k] * m_W[k];
    }

    resid[neq_-3] = r * dudz + u * drdz;
    resid[neq_-2] = u * r * dudz + dpdz + loss;
    resid[neq_-1] = m_gas->density() - r;

    return 0;
}


// ***************************************************************************
//                                  EOF
// ***************************************************************************
