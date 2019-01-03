# -*- coding: utf-8 -*-
import os
import numpy
from CanteraPFR import PyPFR


mech = os.path.join(b"data", b"CT-hydrocarbon-dalmazsi-2017-mech.xml")
# mech = os.path.join(b"data", b"CT-hydrocarbon-norinaga-2009-mech.cti")
phase = b"gas"

Di = 0.028
T0 = 1173.0
p0 = 10000.0
Q0 = 222.0
X0 = b"N2:0.64, C2H2:0.3528, CH3COCH3:6.48e-03, CH4:7.2e-04"

Lr = 0.45
dx = 0.005

rtol = 1.0e-08
atol = 1.0e-08
maxsteps = 50000
dx0 = 1.0e-05


def Tw2(x):
    """ Wall temperature in terms of position. """

    Ta, Tc, Ts = 300.0, 1173.0, 400.0
    x1, x2, m1, m2 = 0.02492942, 0.40810172, 0.78913918, 11.91548263
    term1 = 1 - numpy.exp(-(x / x1) ** m1)
    term2 = 1 - numpy.exp(-(x / x2) ** m2)
    wallT = Ta + (Tc - Ta) * term1 - (Tc - Ts) * term2
    return 0.97 * wallT


def case01():
    """ Usage of isothermal PFR model. """

    rtype = 'isothermal'
    prob = PyPFR(rtype, mech, phase, Di, T0, p0, X0, Q0)
    prob.set_tolerances(rtol, atol)
    prob.set_max_num_steps(maxsteps)
    prob.set_initial_step_size(dx0)
    prob.manage_solution(Lr, dx, f'solution_{rtype}.csv', outfreq=10)


def case02():
    """ Usage of adiabatic PFR model. """

    rtype = 'adiabatic'
    prob = PyPFR(rtype, mech, phase, Di, T0, p0, X0, Q0)
    prob.set_tolerances(rtol, atol)
    prob.set_max_num_steps(maxsteps)
    prob.set_initial_step_size(dx0)
    prob.manage_solution(Lr, dx, f'solution_{rtype}.csv', outfreq=10)


def case03():
    """ Usage of PFR model with functional wall temperature. """

    htc = 10
    T0 = 300.0

    rtype = 'heatwall'
    prob = PyPFR(rtype, mech, phase, Di, T0, p0, X0, Q0, htc=htc, Tw=Tw2)
    prob.set_tolerances(rtol, atol)
    prob.set_max_num_steps(maxsteps)
    prob.set_initial_step_size(dx0)
    prob.manage_solution(Lr, dx, f'solution_{rtype}_2.csv', outfreq=10)


def case04():
    """ Usage of PFR model with constant wall temperature. """

    htc = 10
    Tw1 = 1273.0

    rtype = 'heatwall'
    prob = PyPFR(rtype, mech, phase, Di, T0, p0, X0, Q0, htc=htc, Tw=Tw1)
    prob.set_tolerances(rtol, atol)
    prob.set_max_num_steps(maxsteps)
    prob.set_initial_step_size(dx0)
    prob.manage_solution(Lr, dx, f'solution_{rtype}_1.csv', outfreq=10)


case01()
case02()
case03()
case04()
