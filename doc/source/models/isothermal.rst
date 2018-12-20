Isothermal Plug-Flow Reactor Model
==================================

An isothermal plug-flow reactor (IPFR) is characterized by the absence of energy
equation solution in space and a segregated flow in which each volume of fluid
behaves as a closed system and travels in a single direction without any mixing
with its neighboring elements of fluid. In this sense, it is a system presenting
a Péclet number that is too high, so that convective effects completely shadow
any diffusion of species. The model hereafter described assumes, thus, an
one-dimensional flow in a channel of constant or slowly varying section. The
composition at each coordinate of the reactor is homogeneous (a PFR can be seen
under some circumstances as a series of perfect-stirred reactors) so that any
radial effects are neglected. Finally, a plug-flow model is essentially a steady
state model.

Solution of an IPFR require the coupled integration of mass, momentum, and species
conservation equations. This model implements only gas phase phenomena in such a
way that mass conservation is linked to density and velocity changes only. Species
mole/mass fractions can vary only due to homogeneous processes, which can lead to
changes in average molecular weight and impact directly the momentum and mass
equations. Momentum equation may take viscous effects into account and then be
able to properly represent the pressure loss along the reactor axis.

Model limitations/hypothesis:

1. Isothermal system (not suitable for flames).
2. Laminar flow in circular tubular reactor.
3. Constant reactor cross-section.
4. Homogeneous reactions only.
5. Ideal gas equation of state.

Functionalities to implement:

1. Model is not taking cross-section variations into account.
2. Generalize model by incorporating source (surface) terms.

For more details, please check reference [Ref.1]_.

Model equations
---------------

Throughout the discussion the following symbols convention is adopted:

- :math:`\rho`: density of gas.
- :math:`\mu`: viscosity of gas.
- :math:`u`: velocity in axial direction.
- :math:`R`: radius of reactor.
- :math:`V`: volume of chamber.
- :math:`A`: cross-section of chamber.
- :math:`W`: molecular weight, subscripted to represent a species.
- :math:`\bar{W}`: mean molecular weight.
- :math:`Y`: mass fraction, subscripted to represent a species.
- :math:`\omega`: molar rate of production, subscripted to represent a species.
- :math:`CV`: control volume, used as integration limit.
- :math:`CS`: control surface, used as integration limit.

Mass continuity equation in the absence of wall (source) terms can be derived
by expanding steady state general compressible continuity equation
:math:`\nabla\cdot{}(\rho{}u)`. The full derivation is provided in [Ref.1]_ and
makes use of Reynolds transport theorem. The expansion of this expression over
:math:`z` axis gives the expression implemented in this model:

.. math::

        \rho \frac{\partial u}{\partial z} +
         u \frac{\partial \rho}{\partial z} = 0

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

Coupling of state variables is made by ideal gas law, here implemented as:

.. math::

        p\bar{W} = \rho{}RT

Since we intend to solve the problem as a differential algebraic equations
system, the constraining equations for calculation of initial derivatives
need to be provided. This is done by replacing the symbols in previous equations
by their respective initial states. In order to get the same number of
unknowns as equations, implicit total derivative of state equation is added to
the system. These equations form a set of linear equations that can be solved
by standard linear algebra methods. 

.. math::

        \begin{align*}
          0 & = \rho_{0} u^{\prime}_{0} + u_{0} \rho^{\prime}_{0} \\
          0 & = \rho_{0} u_{0} Y^{\prime}_{k,0} - \omega_{k,0} W_{k} \\
          0 & = \rho_{0} u_{0} u^{\prime}_{0} + p^{\prime}_{0}
                + 8\frac{\mu_{0} u_{0}}{R^2} \\
          0 & = RT\rho^{\prime}_{0} + p_{0}\frac{\bar{W}^{2}_{0}}{W_{k,0}}
                - \bar{W}_{0}p^{\prime}_{0}
        \end{align*}

Module documentation
====================

.. autoclass:: CanteraPFR.models._isothermal.IsothermalPFR
   :members:

References
==========

.. [Ref.1] Chemically Reacting Flow, Theory and Practice. Robert J Kee,
      Michael E. Coltrin, Peter Glarborg. 2003.
