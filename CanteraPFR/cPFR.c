// ***************************************************************************
// Interface for plug-flow reactor models.
//
// Author : Walter Dal'Maz Silva
// Date   : January 14th 2019
// ***************************************************************************

#include "CanteraPFR/HeatWallPFR.hpp"


int cppHeatWallPFR(const std::string & mech, const std::string & phase,
    const std::string & X0, const double Di, const double T0, const double p0,
    const double Q0, const double htc, const std::function<double(double)>& Tw,
    const std::string & saveas, const double length, const double step,
    const double rtol = 1.0e-09, const double atol = 1.0e-15,
    const unsigned maxsteps = 10000, const double initstep = 1.0e-05)
{
    using namespace Cantera;

    try
    {
        HeatWallPFR pfr {mech, phase, Di, T0, p0, X0, Q0, htc, Tw};
        IDA_Solver solver {pfr};

        solver.init(0.0);
        solver.setTolerances(rtol, atol);
        solver.setMaxNumSteps(maxsteps);
        solver.setJacobianType(0);
        solver.setDenseLinearSolver();
        solver.setInitialStepSize(initstep);
        solver.setStopTime(length);

        double x = 0.0;
        int neq = pfr.nEquations();

        while (x < length)
        {
            x = x + std::min(length - x, step);
            solver.solve(x);

            // TODO write to file.
            std::cout << std::scientific << x << " "
                      << solver.solution(neq-4) << " "
                      << solver.solution(neq-3) << " "
                      << solver.solution(neq-2) << " "
                      << solver.solution(neq-1) << " "
                      << pfr.getIntEnergyMass() << " "
                      << std::endl;
        }
    }
    catch (const CanteraError & err)
    {
        std::cerr << err.what() << std::endl;
        return -1;
    }

    return 0;
}


extern "C"
{
    typedef double (*HTCFunc_t)(double);


    int cHeatWallPFR(const char mech[], const char phase[], const char X0[],
        const double Di, const double T0, const double p0, const double Q0,
        const double htc, const HTCFunc_t Tw, const char saveas[],
        const double length, const double step, const double rtol,
        const double atol, const unsigned maxsteps, const double initstep)
    {
        std::string cpp_mech {mech};
        std::string cpp_phase {phase};
        std::string cpp_X0 {X0};
        std::string cpp_saveas {saveas};

        std::cout << "\ncHeatWallPFR interface"
                  << "\nMechanism .............. " << cpp_mech
                  << "\nPhase name ............. " << cpp_phase
                  << "\nInlet composition ...... " << cpp_X0
                  << "\nInlet temperature ...... " << T0
                  << "\nInlet pressure ......... " << p0
                  << "\nInlet flow rate ........ " << Q0
                  << "\nWall HTC ............... " << htc
                  << "\nInlet wall temperature . " << Tw(0.0)
                  << "\nReactor diameter ....... " << Di
                  << "\nReactor length ......... " << length
                  << "\nSave step .............. " << step
                  << "\nRelative tolerance ..... " << rtol
                  << "\nAbsolute tolerance ..... " << atol
                  << "\nInitial step ........... " << initstep
                  << "\nMaximum no. of steps ... " << maxsteps
                  << std::endl;

        return cppHeatWallPFR(cpp_mech, cpp_phase, cpp_X0, Di, T0, p0, Q0, htc,
            Tw, cpp_saveas, length, step, rtol, atol, maxsteps, initstep);
    }


    double Tw(double x)
    {
        const double Ta = 300.0;
        const double Tc = 1173.0;
        const double Ts = 400.0;
        const double x1 = 0.02492942, x2 = 0.40810172;
        const double m1 = 0.78913918, m2 = 11.91548263;
        double term1 = 1 - std::exp(-std::pow(x / x1, m1));
        double term2 = 1 - std::exp(-std::pow(x / x2, m2));
        double wallT = Ta + (Tc - Ta) * term1 - (Tc - Ts) * term2;
        return 0.97 * wallT;
    }


    int test_PFR()
    {
        char saveas[] = "test_PFR.csv";
        char mech[] = "CT-hydrocarbon-dalmazsi-2017-mech.xml";
        char phase[] = "gas";
        char X0[] = "N2:0.64, C2H2:0.36";
        double Di = 0.028;
        double T0 = 300.0;
        double p0 = 5000.0;
        double Q0 = 222.0;
        double htc = 10.0;
        double length = 0.45;
        double step = 0.001;
        double rtol = 1.0e-06;
        double atol = 1.0e-20;
        unsigned maxsteps = 10000;
        double initstep = 1.0e-05;

        return cHeatWallPFR(mech, phase, X0, Di, T0, p0, Q0, htc, Tw,
            saveas, length, step, rtol, atol, maxsteps, initstep);
    }

} // extern "C"

// ***************************************************************************
//                                  EOF
// ***************************************************************************
