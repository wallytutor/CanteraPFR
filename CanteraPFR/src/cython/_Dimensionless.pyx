# -*- coding: utf-8 -*-

import statistics


cdef class Dimensionless(object):
    """ Compute dimensionless numbers for chemical reactors.

    TODO
    ----
    1. Differentiate boundary layer from bulk numbers.
    2. Provide vectorized interfaces for plotting.
    3. Inherit from module so initialization has more sense.

    Parameters
    ----------
    gas : cantera.Solution
        Reactor gas phase.
    L : float
        Characteristic length [m].
    U : float
        Characteristic speed [m/s].
    dT : float
        Temperature gradient [K].
    Q : float
        Volume flow rate [m**3/s].
    """
    def __init__(self, gas, L, U, dT, Q=0):
        alpha_k = gas.thermal_conductivity / (gas.density * gas.cp_mass)
        alpha_D = statistics.mean(gas.mix_diff_coeffs)
        beta = gas.thermal_expansion_coeff
        nu = gas.viscosity / gas.density

        self._Re = self.Re(U, L, nu)
        self._Gr = self.Gr(beta, L, dT, nu)

        self.Pr_a = self.Pr(nu, alpha_k)
        self.Pe_a = self.Pe(U, L, alpha_k)
        self.Ra_a = self._Gr * self.Pr_a  # TODO make function
        self.No_a = [self.Pr_a, self.Pe_a, self.Ra_a]

        self.Pr_D = self.Pr(nu, alpha_D)
        self.Pe_D = self.Pe(U, L, alpha_D)
        self.No_D = [self.Pr_D, self.Pe_D]

        self.props = (gas.density, gas.viscosity, L, U)
        self.No = [self._Re, self._Gr] + self.No_a + self.No_D

    @staticmethod
    def Re(U, L, nu):
        """ Returns Reynolds number.

        Parameters
        ----------
        U : float or array(flow)
            Fluid velocity [m/s].
        L : float
            Problem characteristic length [m].
        nu : float
            Fluid kinematic velocity [m**2/s].
        """
        return U * L / nu

    @staticmethod
    def Gr(beta, L, dT, nu, g=9.81):
        """ Returns Grashof number.

        Parameters
        ----------
        beta : float
            Fluid thermal expansion coefficient [1/K].
        L : float
            Problem characteristic length [m].
        dT : float
            Difference between wall and fluid temperature [K].
        nu : float
            Fluid kinematic velocity [m**2/s].
        """
        return beta * g * L**3 * dT / nu**2

    @staticmethod
    def Pr(nu, kappa):
        """ Returns Prandtl number.

        Parameters
        ----------
        nu : float
            Fluid kinematic viscosity [m**2/s].
        kappa : float.
            Relevant transport coefficient (heat or species) [m**2/s].
        """
        return nu / kappa

    @staticmethod
    def Pe(U, L, kappa):
        """ Returns Peclet number.

        Parameters
        ----------
        U : float or array(flow)
            Fluid velocity [m/s].
        L : float
            Problem characteristic length [m].
        kappa : float.
            Relevant transport coefficient (heat or species) [m**2/s].
        """
        return U * L / kappa

    def __str__(self):
        """ Print formatter. """
        gaspro = (" Gas density....... {:.6e} kg/m\u00B3\n"
                  " Gas viscosity..... {:.6e} Pa.s\n"
                  " Strip velocity.... {:.6e} m/s\n"
                  " Strip length...... {:.6e} m\n\n").format(*self.props)
        output = (" Re   = {:6e}\n"
                  " Gr   = {:6e}\n"
                  " Pr_a = {:6e}\n"
                  " Pe_a = {:6e}\n"
                  " Ra_a = {:6e}\n"
                  " Pr_D = {:6e}\n"
                  " Pe_D = {:6e}\n\n").format(*self.No)
        return gaspro + output

    def report(self):
        """ Print formatted object data. """
        print(self)
