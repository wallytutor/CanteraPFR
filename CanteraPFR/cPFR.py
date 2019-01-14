# -*- coding: utf-8 -*-
"""Call directly from C."""

import os
import time
import ctypes
import numpy


def _set_paths():
    mdir = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
    os.environ['CANTERA_DATA'] = os.path.join(mdir, 'data')
    return mdir


class CPFR(object):
    """ Interface to PFR solver.

    Parameters
    ----------

    """

    def __init__(self, mech, phase, Di, T0, p0, X0, Q0):
        lib = 'libCanteraPFR_shared.so'
        self._clib = ctypes.CDLL(os.path.join(_set_paths(), 'lib', lib))
        self._func = self._clib['cHeatWallPFR']

        self._mech = ctypes.create_string_buffer(mech.encode('utf-8'))
        self._phase = ctypes.create_string_buffer(phase.encode('utf-8'))
        self._X0 = ctypes.create_string_buffer(X0.encode('utf-8'))
        self._Di = ctypes.c_double(Di)
        self._T0 = ctypes.c_double(T0)
        self._p0 = ctypes.c_double(p0)
        self._Q0 = ctypes.c_double(Q0)

    def set_tolerances(self, rtol=1.0e-06, atol=1.0e-15):
        """Set relative and absolute tolerances."""
        self._rtol = ctypes.c_double(rtol)
        self._atol = ctypes.c_double(atol)

    def set_steps(self, initstep=1.0e-05, maxsteps=10000):
        """Set initial step and maximum number of steps."""
        self._initstep = ctypes.c_double(initstep)
        self._maxsteps = ctypes.c_uint(maxsteps)

    def set_wall_conditions(self, htc, tw=None):
        """Set heat transfer coefficient and wall temperature."""
        self._htc = ctypes.c_double(htc)

        # TODO parse sympy compatible strings!
        # TODO accept pure C-functions!

        if tw is None:
            tw = float(self._T0.value)

        if isinstance(tw, (ctypes.c_double, float)):
            self._twptr = tw

            def Tw(x):
                return self._twptr
        else:
            Tw = tw

        self._pytwall = Tw
        ftype = ctypes.CFUNCTYPE(ctypes.c_double, ctypes.c_double)
        self._cctwall = ftype(lambda x: self._pytwall(x))

    def solve(self, length, step=None, saveas=None):
        """Integrate problem."""
        try:
            assert hasattr(self, '_htc'),\
                'Please set HTC before proceeding'

            assert hasattr(self, '_cctwall'),\
                'Please set wall temperature before proceeding'
        except (AssertionError) as err:
            print(err)
            return 1

        if not hasattr(self, '_rtol') or not hasattr(self, '_atol'):
            self.set_tolerances()

        if not hasattr(self, '_initstep') or not hasattr(self, '_maxsteps'):
            self.set_steps()

        if step is None:
            step = ctypes.c_double(length / 100)

        if saveas is None:
            # TODO add datetime.
            saveas = 'test.csv'

        saveas = ctypes.c_wchar_p(saveas)
        length = ctypes.c_double(length)

        t0 = time.time()
        cresult = self._func(self._mech, self._phase, self._X0, self._Di,
                             self._T0, self._p0, self._Q0, self._htc,
                             self._cctwall, saveas, length, step, self._rtol,
                             self._atol, self._maxsteps, self._initstep)
        print(f'Calculation took {time.time()-t0} s')
        return cresult


def test_PFR():
    os.environ['CANTERA_DATA'] = 'data:$CANTERA_DATA'
    lib = os.path.join('lib', 'libCanteraPFR_shared.so')
    module = ctypes.CDLL(lib)
    module['test_PFR']()


def test_full():
    mech = "CT-hydrocarbon-dalmazsi-2017-mech.xml"
    phase = "gas"
    T0 = 300.0
    p0 = 5000.0
    X0 = "N2:0.64, C2H2:0.36"
    Q0 = 222.0
    Di = 0.028
    length = 0.45
    htc = 10.0

    def Tw(x):
        """ Wall temperature in terms of position. """
        Ta, Tc, Ts = 300.0, 1173.0, 400.0
        x1, x2, m1, m2 = 0.02492942, 0.40810172, 0.78913918, 11.91548263
        term1 = 1 - numpy.exp(-(x / x1) ** m1)
        term2 = 1 - numpy.exp(-(x / x2) ** m2)
        wallT = Ta + (Tc - Ta) * term1 - (Tc - Ts) * term2
        return 0.97 * wallT

    r = CPFR(mech, phase, Di, T0, p0, X0, Q0)

    r.set_wall_conditions(htc, tw=Tw)
    r.solve(length)

    r.set_wall_conditions(htc, tw=None)
    r.solve(length)

    r.set_wall_conditions(htc, tw=900.0)
    r.solve(length)


def main():
    """Call to main test."""
    # test_PFR()
    test_full()


if __name__ == '__main__':
    main()
