# -*- coding: utf-8 -*-

import cantera


def _get_mu(gas, mu):
    """ Helper to retrieve viscosity. """

    try:
        mu = gas.viscosity
    except Exception:
        assert mu is not None, 'Missing viscosity'
    return mu


def Re(gas, Uz, D, mu=None):
    """ Returns Reynolds number.

    Parameters
    ----------
    gas : cantera.Solution
        Gas phase with transport properties. If transport manager does not
        implement dynamic viscosity, this shall be provided through `mu`.
    U : float
        Fluid velocity in meters per second.
    D : float
        Problem characteristic diameter in meters.
    mu : float, optional
        Dynamic viscosity in pascal times second. Default is None.

    Raises
    ------
    AssertionError
        If gas has no transport manager and `mu` was not supplied.

    Returns
    -------
    float
        Reynolds number.
    """

    return gas.density * Uz * D / _get_mu(gas, mu)


def Pr(gas, mu=None, k=None):
    """ Returns Prandtl number.

    Parameters
    ----------
    gas : cantera.Solution
        Gas phase with transport properties. If transport manager does not
        implement dynamic viscosity, this shall be provided through `mu`.
    mu : float, optional
        Dynamic viscosity in pascal times second. Default is `None`.
    k : float, optional
        Conductivity in watts per meter times kelvin. Default is `None`.

    Raises
    ------
    AssertionError
        If gas has no transport manager and `mu` was not supplied.

    Returns
    -------
    float
        Prandtl number.
    """

    cp = gas.cp_mass

    def eucken():
        try:
            assert k is not None, 'Missing thermal conductivity'
            return cp * _get_mu(gas, mu) / k
        except Exception as err:
            raise NotImplementedError('Approximation not implemented')
            # FIXME this is wrong!
            # print(f'Using Eucken approximation : {err}')
            # return cp / (cp + 0.00125 * cantera.gas_constant)

    try:
        return cp * _get_mu(gas, mu) / gas.thermal_conductivity
    except Exception:
        return eucken()


def Sc(gas, mu=None, Dab=None):
    """ Returns Schmidt number.

    Parameters
    ----------
    gas : cantera.Solution
        Gas phase with transport properties. If transport manager does not
        implement dynamic viscosity, this shall be provided through `mu`.
    mu : float, optional
        Dynamic viscosity in pascal times second. Default is `None`.
    Dab : float, optional
        Species diffusivity in square meters per second. Default is `None`.

    Raises
    ------
    AssertionError
        If gas has no transport manager and `mu` or `Dab` were not supplied.

    Returns
    -------
    tuple
        Minimum and maximum Schmidt number.
    """

    try:
        kappa_min = min(gas.mix_diff_coeffs)
        kappa_max = max(gas.mix_diff_coeffs)
    except Exception:
        assert Dab is not None, 'Missing diffusivity'
        kappa_min = kappa_max = Dab

    nu = _get_mu(gas, mu) / gas.density
    return nu / kappa_max, nu / kappa_min


def Pe(gas, Uz, D, mu=None, Dab=None, k=None):
    """ Returns Peclet numbers.

    Parameters
    ----------
    gas : cantera.Solution
        Gas phase with transport properties. If transport manager does not
        implement dynamic viscosity, this shall be provided through `mu`.
    U : float
        Fluid velocity in meters per second.
    D : float
        Problem characteristic diameter in meters.
    mu : float, optional
        Dynamic viscosity in pascal times second. Default is `None`.
    Dab : float, optional
        Species diffusivity in square meters per second. Default is `None`.
    k : float, optional
        Conductivity in watts per meter times kelvin. Default is `None`.

    Returns
    -------
    tuple
        A tuple containing the thermal Peclet number and another tuple with
        the minimum and maximum chemical Peclet numbers.
    """

    Re_ev = Re(gas, Uz, D, mu=mu)
    Sc_ev = Sc(gas, mu=mu, Dab=Dab)
    Pe_th = Re_ev * Pr(gas, mu=mu, k=k)
    Pe_ch = (Re_ev * Sc_ev[0], Re_ev * Sc_ev[1])
    return Pe_th, Pe_ch


def Gr(gas, D, Tw, mu=None, g=9.81):
    """ Returns Grashof number.

    Parameters
    ----------
    D : float
        Problem characteristic diameter in meters.
    Tw : float
        Pipe wall temperature in kelvin.
    mu : float, optional
        Dynamic viscosity in pascal times second. Default is `None`.
    g : float, optional
        Gravity acceleration in meters per square second. Default is `9.81`.

    Raises
    ------
    AssertionError
        If gas has no transport manager and `mu` was not supplied.

    Returns
    -------
    float
        Grashof number.
    """

    beta = gas.thermal_expansion_coeff
    nu = _get_mu(gas, mu) / gas.density
    return beta * g * (Tw - gas.T) * D ** 3 / nu ** 2


def Ra(gas, D, Tw, mu=None, k=None, g=9.81):
    """ Returns Rayleigh number.

    Parameters
    ----------
    D : float
        Problem characteristic diameter in meters.
    Tw : float
        Pipe wall temperature in kelvin.
    mu : float, optional
        Dynamic viscosity in pascal times second. Default is `None`.
    k : float, optional
        Conductivity in watts per meter times kelvin. Default is `None`.
    g : float, optional
        Gravity acceleration in meters per square second. Default is `9.81`.

    Returns
    -------
    float
        Rayleigh number.
    """

    return Gr(gas, D, Tw, mu=mu, g=g) * Pr(gas, mu=mu, k=k)
