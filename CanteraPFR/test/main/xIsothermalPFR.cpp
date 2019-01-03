// ***************************************************************************
// Provides an isothermal plug-flow reactor model.
//
// Author : Walter Dal'Maz Silva
// Date   : December 21st 2018
// ***************************************************************************

#include "CanteraPFR/IsothermalPFR.hpp"

int main()
try
{
    std::cout << std::boolalpha
              << "\nStarting solver : " << "IsothermalPFR"
              << "\n Using Sundials : " << CT_SUNDIALS_VERSION
              << "\n Usign LAPACK   : " << bool(CT_SUNDIALS_USE_LAPACK)
              << std::endl;

    // TODO all these parameters should be read from a file.
    double x = 0.00;
    double L = 0.40;
    double dx = 0.01;
    double Di = 0.028;
    double T0 = 1173.0;
    double p0 = 5000.0;
    double Q0 = 222.0;
    std::string X0 = "N2:0.64, C2H2:0.3528, CH3COCH3:6.48e-03, CH4:7.2e-04";
    std::string mech = "test/CT-hydrocarbon-norinaga-2009-mech.xml";
    // std::string mech = "test/CT-hydrocarbon-dalmazsi-2017-mech.xml";
    std::string phase = "gas";
    double rtol = 1.0e-12;
    double atol = 1.0e-20;
    unsigned maxsteps = 50000;
    double dx0 = 1.0e-05;

    // TODO specify type!
    auto t0 = std::chrono::system_clock::now();

    Cantera::IsothermalPFR pfr {mech, phase, Di, T0, p0, X0, Q0};
    Cantera::IDA_Solver solver {pfr};
    solver.init(x);
    solver.setTolerances(rtol, atol);
    solver.setMaxNumSteps(maxsteps);
    solver.setJacobianType(0);
    solver.setDenseLinearSolver();
    solver.setInitialStepSize(dx0);
    solver.setStopTime(L);

    size_t id0 = pfr.getSpeciesIndex("C2H2");
    size_t id1 = pfr.getSpeciesIndex("H2");
    int neq = pfr.nEquations();

    while (x < L)
    {
        x = x + std::min(L-x, dx);
        solver.solve(x);

        std::cout << std::scientific << x << " "
                  << solver.solution(id0) << " "
                  << solver.solution(id1) << " "
                  << solver.solution(neq-3) << " "
                  << solver.solution(neq-2) << " "
                  << solver.solution(neq-1) << " "
                  << std::endl;
    }

    // TODO specify types
    auto t1 = std::chrono::system_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0);
    std::cout << "\nCalculation took " << std::fixed << std::setprecision(3)
              << dt.count() / 1000.0 << " s" << std::endl;

    return 0;
}
catch (Cantera::CanteraError &err)
{
    std::cerr << err.what() << std::endl;
}

// ***************************************************************************
//                                  EOF
// ***************************************************************************
