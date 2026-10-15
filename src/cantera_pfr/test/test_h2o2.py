# -*- coding: utf-8 -*-
import numpy
import cantera
from matplotlib import pyplot
from scipy.integrate import ode
from CanteraPFR.cPFR import CPFR


def test_CPFR():
    mech = 'h2o2.cti'
    phase = 'ohmech'
    T0 = 1500.0
    p0 = 101325.0
    X0 = "H2:2, O2:1, AR:0.1"

    gas = cantera.Solution(mech, phase)

    gas.TPX = 273.15, 101325.0, X0
    rho_ref = gas.density

    gas.TPX = T0, p0, X0
    rho0 = gas.density

    Ac = 1.0e-04
    u0 = 0.006

    # Use flow rate in SCCM!
    Q0 = 6.0e+07 * rho0 * u0 * Ac / rho_ref
    Di = (4 * Ac / numpy.pi) ** 0.5

    length = 1.5e-07
    htc = 0.0
    Tw = T0

    r = CPFR(mech, phase, Di, T0, p0, X0, Q0)
    r.set_wall_conditions(htc, tw=Tw)
    r.solve(length)


def test_VODE():
    """ PFR represented as a marching ODE. """
    mech = 'h2o2.xml'
    gas = cantera.Solution(mech)

    T0 = 1500              # [K]
    p0 = cantera.one_atm   # [Pa]
    gas.TPX = T0, p0, "H2:2, O2:1, AR:0.1"

    length = 1.5e-7   # *approximate* PFR length [m]
    u0 = 0.006        # inflow velocity [m/s]
    Ac = 1.0e-4       # cross-sectional area [m**2]
    Nx = 2000         # number of grid
    dx = length / Nx

    mass_flow_rate = u0 * gas.density * Ac
    N = gas.n_species  # number of gas species
    vecp = numpy.zeros((N+1, 1))

    def rhseqn(t, vec):
        Y = vec[0:N]
        T = vec[-1]

        rho = gas.density
        u = mass_flow_rate / (rho * Ac)
        h_g = gas.partial_molar_enthalpies

        gas.set_unnormalized_mass_fractions(Y)
        gas.TP = T, p0

        wdot_g = gas.net_production_rates
        W_g = gas.molecular_weights

        for k in range(N):
            vecp[k] = wdot_g[k] * W_g[k] / (rho * u)

        vecp[-1] = -numpy.sum(wdot_g * h_g) / (rho * u * gas.cp)
        return vecp

    x = 0.0
    count = 0
    Y0 = numpy.hstack((gas.Y, T0))

    solver = ode(rhseqn)
    solver.set_integrator('vode', method='bdf')
    solver.set_initial_value(Y0, x)

    while x <= length:
        solver.integrate(solver.t + dx)

        if not solver.successful():
            print(f'Failed at x = {solver.t}')

        if not count % 100:
            print(f'x = {solver.t:.3e}\tT = {solver.y[-1]:.3f}')

        x = solver.t
        count += 1

    print(f'x = {solver.t:.3e}\tT = {solver.y[-1]:.3f}')


if __name__ == '__main__':
    test_CPFR()
    test_VODE()
