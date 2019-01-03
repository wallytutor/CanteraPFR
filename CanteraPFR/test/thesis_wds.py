# -*- coding: utf-8 -*-

import os
import time
from CanteraPFR import PyIsothermalPFR as iPFR


mech = b"CT-hydrocarbon-dalmazsi-2017-mech.xml"
phase = b"gas"

Di = 0.028
T0 = 1173.0
p0 = 5000.0
Q0 = 222.0
X0 = b"N2:0.64, C2H2:0.3528, CH3COCH3:6.48e-03, CH4:7.2e-04"

Lr = 0.40
dx = 0.01

rtol = 1.0e-12
atol = 1.0e-20
maxsteps = 50000
dx0 = 1.0e-05

prob = iPFR(mech, phase, Di, T0, p0, X0, Q0)

# Cantera::IDA_Solver solver {pfr};
# solver.init(x);
# solver.setTolerances(rtol, atol);
# solver.setMaxNumSteps(maxsteps);
# solver.setJacobianType(0);
# solver.setDenseLinearSolver();
# solver.setInitialStepSize(dx0);
# solver.setStopTime(L);
#
# size_t id0 = pfr.getSpeciesIndex("C2H2");
# size_t id1 = pfr.getSpeciesIndex("H2");
# int neq = pfr.nEquations();
#
# while (x < L)
# {
#     x = x + std::min(L-x, dx);
#     solver.solve(x);
#
#     std::cout << std::scientific << x << " "
#               << solver.solution(id0) << " "
#               << solver.solution(id1) << " "
#               << solver.solution(neq-3) << " "
#               << solver.solution(neq-2) << " "
#               << solver.solution(neq-1) << " "
#               << std::endl;
# }
