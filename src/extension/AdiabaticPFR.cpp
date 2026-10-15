// ***************************************************************************
// Provides an adiabatic plug-flow reactor model.
//
// Author : Walter Dal'Maz Silva
// Date   : December 21st 2018
// ***************************************************************************

#include "AdiabaticPFR.hpp"

int Cantera::AdiabaticPFR::getInitialConditions(const double t0,
    double *const y, double *const ydot)
{
    const double T0 = m_gas->temperature();
    const double P0 = m_gas->pressure();
    const double rho0 = m_gas->density();
    const double Wavg = m_gas->meanMolecularWeight();
    const double RT = T0 * Cantera::GasConstant;
    const double rho0R = rho0 * Cantera::GasConstant;
    const double rhoUCp = rho0 * m_u0 * m_gas->cp_mass();
    double hdot = 0.0;

    m_gas->getMassFractions(y);
    m_kin->getNetProductionRates(m_wdot.data());
    m_gas->getPartialMolarEnthalpies(m_hbar.data());

    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(neq_, neq_);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(neq_);

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
    A(idx0, idx2) = 0.0;            // p'
    A(idx0, idx3) = 0.0;            // T'

    // Momentum equation elements.
    A(idx1, idx0) = rho0 * m_u0;    // u'
    A(idx1, idx1) = 0.0;            // rho'
    A(idx1, idx2) = 1.0;            // p'
    A(idx1, idx3) = 0.0;            // T'

    // State equation elements.
    A(idx2, idx0) = 0.0;            // u'
    A(idx2, idx1) = RT;             // rho'
    A(idx2, idx2) = -Wavg;          // p'
    A(idx2, idx3) = rho0R;          // T'

    // Energy equation elements.
    A(idx3, idx0) = 0.0;            // u'
    A(idx3, idx1) = 0.0;            // rho'
    A(idx3, idx2) = 0.0;            // p'
    A(idx3, idx3) = rhoUCp;         // T'

    b(idx0) = 0.0;                  // RHS continuity
    b(idx1) = -viscousLoss(m_u0);   // RHS momentum
    b(idx2) = 0.0;                  // RHS state
    b(idx3) = -hdot;                // RHS energy

    Eigen::VectorXd x = A.fullPivLu().solve(b);
    Eigen::VectorXd::Map(ydot, x.rows()) = x;

    return 0;
}

int Cantera::AdiabaticPFR::evalResidNJ(const double t,
    const double delta_t, const double* const y,
    const double* const ydot, double* const resid,
    const int id_x, const double delta_x)
{
    const double u = y[idx0];
    const double r = y[idx1];
    const double p = y[idx2];
    const double T = y[idx3];

    const double dudz = ydot[idx0];
    const double drdz = ydot[idx1];
    const double dpdz = ydot[idx2];
    const double dTdz = ydot[idx3];

    m_gas->setMassFractions_NoNorm(y);
    m_gas->setState_TP(T, p);
    m_kin->getNetProductionRates(m_wdot.data());
    m_gas->getPartialMolarEnthalpies(m_hbar.data());

    const double Cp = m_gas->cp_mass();
    double hdot = 0.0;

    for (unsigned k = 0; k != idx0; ++k)
    {
        resid[k] = u * r * ydot[k] - m_wdot[k] * m_W[k];
        hdot += m_wdot[k] * m_hbar[k];
    }

    resid[idx0] = r * dudz + u * drdz;
    resid[idx1] = u * r * dudz + dpdz + viscousLoss(u);
    resid[idx2] = m_gas->density() - r;
    resid[idx3] = r * u * Cp * dTdz + hdot;

    return 0;
}
