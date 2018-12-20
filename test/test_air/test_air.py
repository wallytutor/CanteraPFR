# -*- coding: utf-8 -*-
"""
Test CanteraPFR isothermal PFR model.
"""

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.pardir, os.pardir, os.pardir)))

import numpy
from pandas import DataFrame
from matplotlib import pyplot
from CanteraPFR import IsothermalPFR

Q0 = 500
T0 = 1723.0
p0 = 0.1 * 101325
X0 = {'O2': 0.70, 'N2':0.28, 'AR': 0.01, 'NO': 0.01}

Di = 0.028
Ac = numpy.pi * Di ** 2 / 4

kwargs = dict(first_step_size=1.0e-08, atol=1.0e-20,
              rtol=1.0e-13, max_steps=5000)

solver = IsothermalPFR('air.cti', (T0, p0, X0), Q0, Ac,
                       phase='air', trans='Mix')
solution = solver.integrate(0.40, 0.002, **kwargs)
solution.to_csv('test_air.csv')

if solution.idasol.flag == 0:
    pyplot.close('all')
    pyplot.style.use('bmh')

    x = solution.position
    fig, ax = pyplot.subplots(2, 2, figsize=(15, 9), dpi=150)

    ax[0, 0].plot(x, solution.velocity, color='C0')
    ax[0, 0].set_xlabel('Distance ($m$)')
    ax[0, 0].set_ylabel('Velocity ($m\\,s^{-1}$)')

    ax[0, 1].plot(x, solution.density, color='C1')
    ax[0, 1].set_xlabel('Distance ($m$)')
    ax[0, 1].set_ylabel('Density ($\\mathregular{kg\\,m^{-3}}$)')
    ax[0, 1].ticklabel_format(axis='y', style='sci', scilimits=(-2, 2))

    ax[1, 0].plot(x, solution.pressure, color='C2')
    ax[1, 0].set_xlabel('Distance (m)')
    ax[1, 0].set_ylabel('Pressure (Pa)')

    ax[1, 1].plot(x, solution.mole_fraction('O'), label='$O$')
    ax[1, 1].plot(x, solution.mole_fraction('NO2'), label='$NO_2$')
    ax[1, 1].legend(loc='best')
    ax[1, 1].set_xlabel('Distance (m)')
    ax[1, 1].set_ylabel('Mole Fraction')

    fig.tight_layout()
    pyplot.savefig('test_air.png', dpi=300)
