# -*- coding: utf-8 -*-
"""
Test CanteraPFR isothermal PFR model.
"""

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.pardir, os.pardir, os.pardir)))

import numpy
from matplotlib import pyplot
from CanteraPFR import IsothermalPFR

mechanisms = ['CT-hydrocarbon-dalmazsi-2017-mech.cti',
              'CT-hydrocarbon-norinaga-2009-mech.cti']
mechanisms = [os.path.join(os.pardir, m) for m in mechanisms]

for mech in mechanisms:
    basename = os.path.basename(mech)
    saveas = 'cantera_pfr-{}.png'.format(basename.split('.')[0])

    kwargs = dict(first_step_size=1.0e-06, atol=1.0e-10, rtol=1.0e-08,
                  max_steps=1000)

    x0 = 0.360
    x1 = 0.018 * x0
    x2 = 0.002 * x0

    T0 = 1173.0
    p0 = 5000.0
    X0 = {'N2': 0.64, 'C2H2': x0, 'CH3COCH3': x1, 'CH4': x2}

    Q = 222

    Di = 0.028
    Ac = numpy.pi * Di ** 2 / 4

    solver = IsothermalPFR(mech, (T0, p0, X0), Q, Ac)
    solution = solver.integrate(0.40, 0.005, **kwargs)
    x = solution.position

    if solution.idasol.flag == 0:
        pyplot.close('all')
        pyplot.style.use('bmh')

        fig, ax = pyplot.subplots(3, 2, figsize=(15, 9), dpi=150)

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

        ax[2, 0].plot(x, solution.mole_fraction('C4H2'), label='$C_4H_2$')
        ax[2, 0].plot(x, solution.mole_fraction('C4H4'), label='$C_4H_4$')
        ax[2, 0].plot(x, solution.mole_fraction('CO'), label='$CO$')
        ax[2, 0].legend(fontsize=7, loc='upper right')
        ax[2, 0].set_xlabel('Distance ($m$)')
        ax[2, 0].set_ylabel('Mole Fraction')

        ax[1, 1].plot(x, solution.mole_fraction('C2H2'), label='$C_2H_2$')
        ax[1, 1].legend(loc='best')
        ax[1, 1].set_xlabel('Distance (m)')
        ax[1, 1].set_ylabel('Mole Fraction')

        ax[2, 1].plot(x, solution.mole_fraction('H2'), label='$H_2$')
        ax[2, 1].plot(x, solution.mole_fraction('CH4'), label='$CH_4$')
        ax[2, 1].plot(x, solution.mole_fraction('C6H6'), label='$C_6H_6$')
        ax[2, 1].legend(fontsize=7, loc='upper right')
        ax[2, 1].set_xlabel('Distance ($m$)')
        ax[2, 1].set_ylabel('Mole Fraction')

        fig.tight_layout()
        pyplot.savefig(saveas, dpi=300)
