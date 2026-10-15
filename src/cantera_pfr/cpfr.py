# -*- coding: utf-8 -*-
""" High-level interface to plug-flow reactor models. """

import time
from pathlib import Path
from typing import Callable, Union

from .ct_pfr import PyPFR

TwType = Union[float, Callable[[float], float]]


class CPFR:
    """ Plug-flow reactor solver interface compatible with legacy CPFR API.

    Parameters
    ----------
    mech : str
        Mechanism file path or name.
    phase : str
        Phase name in the mechanism file.
    Di : float
        Reactor inner diameter in meters.
    T0 : float
        Inlet temperature in kelvin.
    p0 : float
        Inlet pressure in pascal.
    X0 : str
        Inlet composition in Cantera format.
    Q0 : float
        Inlet volume flow rate in sccm.
    """

    __slots__ = (
        'mech',
        'phase',
        'Di',
        'T0',
        'p0',
        'X0',
        'Q0',
        'rtol',
        'atol',
        'initstep',
        'maxsteps',
        'htc',
        'tw',
    )

    def __init__(
        self,
        mech: str,
        phase: str,
        Di: float,
        T0: float,
        p0: float,
        X0: str,
        Q0: float
    ) -> None:
        self.mech = str(mech)
        self.phase = str(phase)
        self.Di = float(Di)
        self.T0 = float(T0)
        self.p0 = float(p0)
        self.X0 = str(X0)
        self.Q0 = float(Q0)

        self.rtol = 1.0e-06
        self.atol = 1.0e-15
        self.initstep = 1.0e-05
        self.maxsteps = 10000
        self.htc = 0.0
        self.tw: TwType = float(T0)

    def set_tolerances(
        self,
        rtol: float = 1.0e-06,
        atol: float = 1.0e-15
    ) -> None:
        """ Set relative and absolute tolerances.

        Parameters
        ----------
        rtol : float = 1.0e-06
            Relative tolerance.
        atol : float = 1.0e-15
            Absolute tolerance.
        """
        self.rtol = float(rtol)
        self.atol = float(atol)

    def set_steps(
        self,
        initstep: float = 1.0e-05,
        maxsteps: int = 10000
    ) -> None:
        """ Set initial step size and maximum number of steps.

        Parameters
        ----------
        initstep : float = 1.0e-05
            Initial step length in meters.
        maxsteps : int = 10000
            Maximum number of internal solver steps.
        """
        self.initstep = float(initstep)
        self.maxsteps = int(maxsteps)

    def set_wall_conditions(
        self,
        htc: float,
        tw: Union[float, Callable[[float], float], None] = None
    ) -> None:
        """ Set wall heat transfer coefficient and wall temperature profile.

        Parameters
        ----------
        htc : float
            Wall heat transfer coefficient in W/(m^2 K).
        tw : float or callable, optional
            Wall temperature in kelvin (constant or function of position).
        """
        self.htc = float(htc)

        if tw is None:
            self.tw = float(self.T0)
        else:
            self.tw = tw

    def solve(
        self,
        length: float,
        step: Union[float, None] = None,
        saveas: Union[str, Path, None] = None
    ) -> int:
        """ Integrate problem along reactor length.

        Parameters
        ----------
        length : float
            Total reactor length in meters.
        step : float, optional
            Saving step in meters. Defaults to length / 100.
        saveas : str or Path, optional
            Path to output CSV results file.

        Returns
        -------
        int
            Exit status (0 on success).
        """
        dx = float(length / 100.0) if step is None else float(step)

        out_path = Path(saveas) if saveas is not None else Path('test.csv')

        print(
            f'\nCPFR solver interface'
            f'\nMechanism .............. {self.mech}'
            f'\nPhase name ............. {self.phase}'
            f'\nInlet composition ...... {self.X0}'
            f'\nInlet temperature ...... {self.T0}'
            f'\nInlet pressure ......... {self.p0}'
            f'\nInlet flow rate ........ {self.Q0}'
            f'\nWall HTC ............... {self.htc}'
            f'\nReactor diameter ....... {self.Di}'
            f'\nReactor length ......... {length}'
            f'\nSave step .............. {dx}'
            f'\nRelative tolerance ..... {self.rtol}'
            f'\nAbsolute tolerance ..... {self.atol}'
            f'\nInitial step ........... {self.initstep}'
            f'\nMaximum no. of steps ... {self.maxsteps}'
        )

        pfr = PyPFR(
            rtype='heatwall',
            mech=self.mech,
            phase=self.phase,
            Di=self.Di,
            T0=self.T0,
            p0=self.p0,
            X0=self.X0,
            Q0=self.Q0,
            htc=self.htc,
            Tw=self.tw,
        )

        pfr.set_tolerances(self.rtol, self.atol)
        pfr.set_max_num_steps(self.maxsteps)
        pfr.set_initial_step_size(self.initstep)
        pfr.set_stop_position(length)

        t0 = time.time()
        x = 0.0

        while x < length:
            x = x + min(length - x, dx)
            retcode = pfr.solve(x)

            if retcode != 0:
                print(f'Solver error {retcode} at x = {x}')
                return retcode

        if saveas is not None:
            pfr.write_results(str(out_path))

        print(f'Calculation took {time.time() - t0:.5f} s')
        return 0


cPFR = CPFR

__all__ = ['CPFR', 'cPFR']
