# -*- coding: utf-8 -*-
from typing import Sequence
import cantera as ct


class Dimensionless:
    """ Evaluate dimensionless numbers for reacting flows.

    Parameters
    ----------
    gas : ct.Solution
        Cantera solution object representing gas phase.

    Returns
    -------
    Dimensionless
        Configured evaluator instance.
    """

    __slots__ = ('_gas',)

    def __init__(self, gas: ct.Solution) -> None:
        # XXX no documentation here, this is intentional
        self._gas = gas

    def _get_mu(self, mu):
        try:
            return self._gas.viscosity
        except Exception:
            if mu is None:
                raise ValueError('Missing dynamic viscosity')

            return mu

    def re(
        self,
        uz: float,
        d: float,
        mu: float | None = None
    ) -> float:
        """ Compute Reynolds number.

        Parameters
        ----------
        uz : float
            Fluid velocity in meters per second.
        d : float
            Characteristic diameter in meters.
        mu : float | None = None
            Dynamic viscosity in pascal second.

        Returns
        -------
        float
            Reynolds number.
        """
        return self._gas.density * uz * d / self._get_mu(mu)

    def pr(
        self,
        mu: float | None = None,
        k: float | None = None
    ) -> float:
        """ Compute Prandtl number.

        Parameters
        ----------
        mu : float | None = None
            Dynamic viscosity in pascal second.
        k : float | None = None
            Thermal conductivity in W/(m*K).

        Returns
        -------
        float
            Prandtl number.
        """
        cp = self._gas.cp_mass

        try:
            return cp * self._get_mu(mu) / self._gas.thermal_conductivity
        except Exception:
            if k is None:
                raise ValueError('Missing thermal conductivity')

            return cp * self._get_mu(mu) / k

    def sc(
        self,
        mu: float | None = None,
        dab: float | None = None
    ) -> tuple[float, float]:
        """ Compute Schmidt number bounds.

        Parameters
        ----------
        mu : float | None = None
            Dynamic viscosity in pascal second.
        dab : float | None = None
            Species diffusivity in square meters per second.

        Returns
        -------
        tuple[float, float]
            Minimum and maximum Schmidt numbers.
        """
        try:
            diff_coeffs: Sequence[float] = self._gas.mix_diff_coeffs
            kappa_min = float(min(diff_coeffs))
            kappa_max = float(max(diff_coeffs))
        except Exception:
            if dab is None:
                raise ValueError('Missing diffusivity')

            kappa_min = float(dab)
            kappa_max = float(dab)

        nu = self._get_mu(mu) / self._gas.density
        return nu / kappa_max, nu / kappa_min

    def pe(
        self,
        uz: float,
        d: float,
        mu: float | None = None,
        dab: float | None = None,
        k: float | None = None
    ) -> tuple[float, tuple[float, float]]:
        """ Compute thermal and chemical Peclet numbers.

        Parameters
        ----------
        uz : float
            Fluid velocity in meters per second.
        d : float
            Characteristic diameter in meters.
        mu : float | None = None
            Dynamic viscosity in pascal second.
        dab : float | None = None
            Species diffusivity in square meters per second.
        k : float | None = None
            Thermal conductivity in W/(m*K).

        Returns
        -------
        tuple[float, tuple[float, float]]
            Thermal and chemical Peclet numbers.
        """
        re_val = self.re(uz, d, mu=mu)
        sc_val = self.sc(mu=mu, dab=dab)
        pr_val = self.pr(mu=mu, k=k)
        pe_th = re_val * pr_val
        pe_ch = (re_val * sc_val[0], re_val * sc_val[1])
        return pe_th, pe_ch

    def gr(
        self,
        d: float,
        tw: float,
        mu: float | None = None,
        g: float = 9.81
    ) -> float:
        """ Compute Grashof number.

        Parameters
        ----------
        d : float
            Characteristic diameter in meters.
        tw : float
            Wall temperature in kelvin.
        mu : float | None = None
            Dynamic viscosity in pascal second.
        g : float = 9.81
            Gravity acceleration in m/s^2.

        Returns
        -------
        float
            Grashof number.
        """
        beta = self._gas.thermal_expansion_coeff
        nu = self._get_mu(mu) / self._gas.density
        return beta * g * (tw - self._gas.T) * (d ** 3) / (nu ** 2)

    def ra(
        self,
        d: float,
        tw: float,
        mu: float | None = None,
        k: float | None = None,
        g: float = 9.81
    ) -> float:
        """ Compute Rayleigh number.

        Parameters
        ----------
        d : float
            Characteristic diameter in meters.
        tw : float
            Wall temperature in kelvin.
        mu : float | None = None
            Dynamic viscosity in pascal second.
        k : float | None = None
            Thermal conductivity in W/(m*K).
        g : float = 9.81
            Gravity acceleration in m/s^2.

        Returns
        -------
        float
            Rayleigh number.
        """
        return self.gr(d, tw, mu=mu, g=g) * self.pr(mu=mu, k=k)


def Re(
        gas: ct.Solution,
        Uz: float,
        D: float,
        mu: float | None = None
    ) -> float:
    """ Compute Reynolds number.

    Parameters
    ----------
    gas : ct.Solution
        Cantera gas solution.
    Uz : float
        Fluid velocity in meters per second.
    D : float
        Characteristic diameter in meters.
    mu : float | None = None
        Dynamic viscosity in pascal second.

    Returns
    -------
    float
        Reynolds number.
    """
    return Dimensionless(gas).re(Uz, D, mu=mu)


def Pr(
        gas: ct.Solution,
        mu: float | None = None,
        k: float | None = None
    ) -> float:
    """ Compute Prandtl number.

    Parameters
    ----------
    gas : ct.Solution
        Cantera gas solution.
    mu : float | None = None
        Dynamic viscosity in pascal second.
    k : float | None = None
        Thermal conductivity in W/(m*K).

    Returns
    -------
    float
        Prandtl number.
    """
    return Dimensionless(gas).pr(mu=mu, k=k)


def Sc(
        gas: ct.Solution,
        mu: float | None = None,
        Dab: float | None = None
    ) -> tuple[float, float]:
    """ Compute Schmidt number bounds.

    Parameters
    ----------
    gas : ct.Solution
        Cantera gas solution.
    mu : float | None = None
        Dynamic viscosity in pascal second.
    Dab : float | None = None
        Diffusivity in square meters per second.

    Returns
    -------
    tuple[float, float]
        Minimum and maximum Schmidt numbers.
    """
    return Dimensionless(gas).sc(mu=mu, dab=Dab)


def Pe(
        gas: ct.Solution,
        Uz: float,
        D: float,
        mu: float | None = None,
        Dab: float | None = None,
        k: float | None = None
    ) -> tuple[float, tuple[float, float]]:
    """ Compute thermal and chemical Peclet numbers.

    Parameters
    ----------
    gas : ct.Solution
        Cantera gas solution.
    Uz : float
        Fluid velocity in meters per second.
    D : float
        Characteristic diameter in meters.
    mu : float | None = None
        Dynamic viscosity in pascal second.
    Dab : float | None = None
        Diffusivity in square meters per second.
    k : float | None = None
        Thermal conductivity in W/(m*K).

    Returns
    -------
    tuple[float, tuple[float, float]]
        Thermal and chemical Peclet numbers.
    """
    return Dimensionless(gas).pe(Uz, D, mu=mu, dab=Dab, k=k)


def Gr(
        gas: ct.Solution,
        D: float,
        Tw: float,
        mu: float | None = None,
        g: float = 9.81
    ) -> float:
    """ Compute Grashof number.

    Parameters
    ----------
    gas : ct.Solution
        Cantera gas solution.
    D : float
        Characteristic diameter in meters.
    Tw : float
        Wall temperature in kelvin.
    mu : float | None = None
        Dynamic viscosity in pascal second.
    g : float = 9.81
        Gravity acceleration in m/s^2.

    Returns
    -------
    float
        Grashof number.
    """
    return Dimensionless(gas).gr(D, Tw, mu=mu, g=g)


def Ra(
        gas: ct.Solution,
        D: float,
        Tw: float,
        mu: float | None = None,
        k: float | None = None,
        g: float = 9.81
    ) -> float:
    """ Compute Rayleigh number.

    Parameters
    ----------
    gas : ct.Solution
        Cantera gas solution.
    D : float
        Characteristic diameter in meters.
    Tw : float
        Wall temperature in kelvin.
    mu : float | None = None
        Dynamic viscosity in pascal second.
    k : float | None = None
        Thermal conductivity in W/(m*K).
    g : float = 9.81
        Gravity acceleration in m/s^2.

    Returns
    -------
    float
        Rayleigh number.
    """
    return Dimensionless(gas).ra(D, Tw, mu=mu, k=k, g=g)
