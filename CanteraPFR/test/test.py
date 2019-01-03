# -*- coding: utf-8 -*-
from CanteraPFR import PyPFR


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

rtype = 'isothermal'
prob = PyPFR(rtype, mech, phase, Di, T0, p0, X0, Q0)
prob.set_tolerances(rtol, atol)
prob.set_max_num_steps(maxsteps)
prob.set_initial_step_size(dx0)
prob.manage_solution(Lr, dx, f'solution_{rtype}.csv', outfreq=10)

rtype = 'adiabatic'
prob = PyPFR(rtype, mech, phase, Di, T0, p0, X0, Q0)
prob.set_tolerances(rtol, atol)
prob.set_max_num_steps(maxsteps)
prob.set_initial_step_size(dx0)
prob.manage_solution(Lr, dx, f'solution_{rtype}.csv', outfreq=10)

htc = 10
Tw1 = 973.0

rtype = 'heatwall'
prob = PyPFR(rtype, mech, phase, Di, T0, p0, X0, Q0, htc=htc, Tw=Tw1)
prob.set_tolerances(rtol, atol)
prob.set_max_num_steps(maxsteps)
prob.set_initial_step_size(dx0)
prob.manage_solution(Lr, dx, f'solution_{rtype}_1.csv', outfreq=10)

# def Tw2(x):
#     return 973.0
#
# rtype = 'heatwall'
# prob = PyPFR(rtype, mech, phase, Di, T0, p0, X0, Q0, htc=htc, Tw=Tw2)
# prob.set_tolerances(rtol, atol)
# prob.set_max_num_steps(maxsteps)
# prob.set_initial_step_size(dx0)
# prob.manage_solution(Lr, dx, f'solution_{rtype}_2.csv', outfreq=10)
