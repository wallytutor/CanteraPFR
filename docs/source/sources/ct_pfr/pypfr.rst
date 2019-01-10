Plug-Flow Reactor Model
=======================

A plug-flow reactor (PFR) model represents the limit of a tubular reactor for
which Péclet number in axial direction is too high so that mixing can be
neglected over flow direction (so that convective effects completely shadow
any diffusion of species). That implies segregated flow in which each volume
of fluid behaves as a closed system and travels in a single direction without
any mixing with its neighboring elements of fluid. The composition at each
coordinate of the reactor is homogeneous (a PFR can be seen under some cases
as a series of perfect-stirred reactors) so that any radial effects are
neglected. It is essentially a one-dimensional model and for most purposes is
treated on stationary state, since the non-mixing criterion makes the transient
approach meaningless. This package deals with PFR's that are both constant and
smoothly varying (**to be implemented**) in cross-sectional area (assumed
circular). Each of these cases can be treated according to slightly different
formulations, what includes:

1. isothermal PFR: is characterized by the absence of energy equation solution.
2. adiabatic PFR: heat released by reactions is not exchanged with environment.
3. imposed wall temperature PFR: simplified convection heat exchanged with environment.

**Currently no surface exchanges are implemented, that must be a next step.
This will require a way to generalize the interface to place a catalyst or
reacting interface inside the system, so that it can be usable in many cases.**

Numerical integration
---------------------

Solution of a PFR requires at least the coupled integration of mass, momentum,
and species conservation equations. For those cases for which temperature is
solved, energy conservation must be added to this set. These equations are then
coupled by means of ideal gas law. Following, equations solved by this package
are presented. For more details, please check reference [Ref.1]_.

Throughout the discussion the following symbols convention is adopted:

- :math:`\rho`: density of gas.
- :math:`\mu`: viscosity of gas.
- :math:`T`: temperature of gas.
- :math:`u`: velocity in axial direction.
- :math:`R`: radius of reactor.
- :math:`A`: cross-section of chamber.
- :math:`V`: volume of chamber.
- :math:`\mathcal{P}_{a}`: perimeter over cross-section :math:`A` of tube.
- :math:`\hat{h}`: global convective heat transfer coefficient.
- :math:`T_{w}`: local wall temperature.
- :math:`C_{p}`: mixture constant pressure heat capacity.
- :math:`h_{k}`: species :math:`k` enthalpy per unit mass.
- :math:`W`: molecular weight, with subscript to represent a species.
- :math:`\bar{W}`: mean molecular weight.
- :math:`Y`: mass fraction, with subscript to represent a species.
- :math:`\omega`: molar rate of production, with subscript to represent a species.
- :math:`CV`: control volume, used as integration limit.
- :math:`CS`: control surface, used as integration limit.

Mass conservation
-----------------

Mass continuity equation in the absence of wall (source) terms can be derived
by expanding steady state general compressible continuity equation
:math:`\nabla\cdot{}(\rho{}u)`. The full derivation is provided in [Ref.1]_ and
makes use of Reynolds transport theorem. The expansion of this expression over
:math:`z` axis gives the expression implemented in this package:

.. math::

        \rho \frac{\partial u}{\partial z} +
         u \frac{\partial \rho}{\partial z} = 0

Species conservation
--------------------

Species conservation starts by applying Reynolds transport theorem on the overall
system species mass conservation. By using divergence theorem (here specialized
in one dimension) this can be written as:

.. math::

        \int_{CV}\frac{\partial \rho{}Y_{k}}{\partial t}dV +
        \int_{CS}\rho{}Y_{k}u{}dA =
        \int_{CV}\omega_{k}W_{k}dV

Assuming steady state and no variations in channel:

.. math::

        \frac{\partial \rho{}uY_{k}}{\partial z} - \omega_{k}W_{k} = 0

The axial derivative can be simplified using the overall mass balance leading to
the expression that is implemented:

.. math::

        \rho{}u\frac{\partial Y_{k}}{\partial z} - \omega_{k}W_{k} = 0

Momentum conservation
---------------------

Momentum conservation is derived for a Newtonian fluid by taking into account
viscous effects (drag over channel walls) and pressure acting at inlet. Viscous
effects are introduced through a model equation as discussed in [Ref.1]_. The
equation is deriving by assuming a smooth pressure gradient and a friction factor
applicable for the laminar limit. Hydraulic diameter proposed in section 16.2 of
[Ref.1]_ is used to simplify the equation to its final form as:

.. math::

        \rho{}u\frac{\partial u}{\partial z} +
        \frac{\partial p}{\partial z} +
        8\frac{\mu u}{R^2} = 0

State equation
--------------

Coupling of state variables is made by ideal gas law, here implemented as:

.. math::

        p\bar{W} = \rho{}RT

Energy conservation
-------------------

Except for the isothermal case, energy equation is included in a general form as
follows. Notice that for adiabatic reactors :math:`\hat{h}=0`, what makes the
convective term disappear from the equation.

.. math::

    \rho{}u{}C_{p}\frac{\partial T}{\partial z} +
    \sum_{k}h_{k}\omega_{k}W_{k} -
    \mathcal{P}_{a}\hat{h}(T_{w}-T) = 0

Problem initialization
----------------------

Since we intend to solve the problem as a differential algebraic equations
system, the constraining equations for calculation of initial derivatives
is required. This is done by replacing the symbols in previous equations
by their respective initial states. In order to get the same number of
unknowns as equations, implicit total derivative of state equation is added to
the system. These equations form a set of linear equations that can be solved
by standard linear algebra methods.

.. math::

        \begin{align*}
          \omega_{k,0} W_{k}
            & = \rho_{0} u_{0} Y^{\prime}_{k,0} \\
          0 & = \rho_{0} u^{\prime}_{0} + u_{0} \rho^{\prime}_{0} \\
          -8\frac{\mu_{0} u_{0}}{R^2}
            & = \rho_{0} u_{0} u^{\prime}_{0} + p^{\prime}_{0} \\
          0 & = p_{0}\frac{\bar{W}^{2}_{0}}{W_{k,0}}Y^{\prime}_{k,0}
                + RT\rho^{\prime}_{0}
                - \bar{W}_{0}p^{\prime}_{0}
                + \rho_{0}RT^{\prime}_{0} \\
          -\sum_{k}h_{k,0}\omega_{k,0} W_{k}
          +\mathcal{P}_{a}\hat{h}(T_{w,0}-T_{0})
            & = \rho_{0} u_{0} C_{p,0} T^{\prime}_{0}
        \end{align*}

Notice that here equations were expressed in most general form considered in this
package. Temperature derivatives have to be eliminated for isothermal case and
convective heat transfer coefficient set to zero in adiabatic limit.

Reference API
-------------

.. module:: CanteraPFR.ct_pfr

.. autoclass:: PyPFR
    :members:
    :undoc-members:
    :show-inheritance:
