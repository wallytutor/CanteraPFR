// ***************************************************************************
// Provides a plug-flow reactor model with wall exchanges.
//
// Author : Walter Dal'Maz Silva
// Date   : December 31st 2018
// ***************************************************************************

#include "CanteraPFR/HeatWallPFR.hpp"

void example01()
{
    // TODO all these parameters should be read from a file.
    double x = 0.00;
    double L = 0.40;
    double dx = 0.010;
    double Di = 0.028;
    double T0 = 800;
    double p0 = 5000.0;
    double Q0 = 222.0;
    double htc = 10.0;
    double Tw = 1173.0;
    // std::function<double(double)> Tw = [](double) { return 1173.0; };
    std::string X0 = "N2:0.64, C2H2:0.3528, CH3COCH3:6.48e-03, CH4:7.2e-04";
    // std::string mech = "test/CT-hydrocarbon-norinaga-2009-mech.xml";
    std::string mech = "test/CT-hydrocarbon-dalmazsi-2017-mech.xml";
    std::string phase = "gas";
    double rtol = 1.0e-12;
    double atol = 1.0e-20;
    unsigned maxsteps = 50000;
    double dx0 = 1.0e-05;

    // TODO specify type!
    auto t0 = std::chrono::system_clock::now();

    Cantera::HeatWallPFR pfr {mech, phase, Di, T0, p0, X0, Q0, htc, Tw};
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
                  << solver.solution(neq-4) << " "
                  << solver.solution(neq-3) << " "
                  << solver.solution(neq-2) << " "
                  << solver.solution(neq-1) << " "
                  << pfr.getIntEnergyMass() << " "
                  << std::endl;
    }

    // TODO specify types
    auto t1 = std::chrono::system_clock::now();
    auto dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0);
    std::cout << "\nCalculation took " << std::fixed << std::setprecision(3)
              << dt.count() / 1000.0 << " s" << std::endl;
}

void example_thesis()
{

}

int main()
try
{
    std::cout << std::boolalpha
              << "\nStarting solver : " << "HeatWallPFR"
              << "\n Using Sundials : " << CT_SUNDIALS_VERSION
              << "\n Usign LAPACK   : " << bool(CT_SUNDIALS_USE_LAPACK)
              << std::endl;

    // example01();

    // pars = {
    //     '773'  : (0.04132785, 0.36586941, 1.92089872, 12.41516606),
    //     '873'  : (0.03457862, 0.39032227, 1.41582889, 9.79102679),
    //     '973'  : (0.02537489, 0.39703098, 0.99659743, 9.77523826),
    //     '1023' : (0.02528152, 0.40339555, 0.88494798, 10.55513796),
    //     '1073' : (0.02507178, 0.40847247, 0.81631547, 11.98899245),
    //     '1123' : (0.02497517, 0.40832661, 0.80065655, 11.97005813),
    //     '1173' : (0.02492942, 0.40810172, 0.78913918, 11.91548263),
    //     '1223' : (0.02596356, 0.40572591, 0.85168097, 11.01722351),
    //     '1273' : (0.02682903, 0.40342913, 0.91051192, 10.36909121)
    // }

    const double Ta = 300.0;
    const double Tc = 1173.0;
    const double Ts = 400.0;
    const double x1 = 0.02492942, x2 = 0.40810172;
    const double m1 = 0.78913918, m2 = 11.91548263;

    std::function<double(double)> Tw = [&](double x) {
        double term1 = 1 - std::exp(-std::pow(x / x1, m1));
        double term2 = 1 - std::exp(-std::pow(x / x2, m2));
        double wallT = Ta + (Tc - Ta) * term1 - (Tc - Ts) * term2;
        return 0.97 * wallT;
    };

    double x = 0.00;
    double L = 0.45;
    double dx = 0.005;
    double Di = 0.028;
    double T0 = 300.00;
    double p0 = 10000.0;
    double Q0 = 222.0;
    double htc = 10.0;

    std::string X0 = "N2:0.64, C2H2:0.3528, CH3COCH3:6.48e-03, CH4:7.2e-04";

    std::string mech;
    mech = "test/CT-hydrocarbon-norinaga-2009-mech.xml";
    // mech = "test/CT-hydrocarbon-dalmazsi-2017-mech.xml";
    std::string phase = "gas";

    double rtol = 1.0e-06;
    double atol = 1.0e-15;
    unsigned maxsteps = 50000;
    double dx0 = 1.0e-05;

    // TODO specify type!
    auto t0 = std::chrono::system_clock::now();

    Cantera::HeatWallPFR pfr {mech, phase, Di, T0, p0, X0, Q0, htc, Tw};
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
                  << solver.solution(neq-4) << " "
                  << solver.solution(neq-3) << " "
                  << solver.solution(neq-2) << " "
                  << solver.solution(neq-1) << " "
                  << pfr.getIntEnergyMass() << " "
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
