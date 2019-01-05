Dimensionless Numbers
=====================

.. module:: CanteraPFR.CanteraAux

When dealing with reactor models, it is always useful to be able to quickly
compute approximate dimensionless numbers for the studied case. In what follows
we provide a brief description of these quantities as implemented in this
package. Definitions might vary according to the author or field of application
and here we follow solely [Ref.2]_.

Reynolds Number
---------------

*This dimensionless group is named for Osborne Reynolds (1842-1912), professor
of engineering at the University of Manchester. He studied the laminar-turbulent
transition, turbulent heat transfer, and theory of lubrication* [Ref.2]_. In
general we denote Reynolds number by :math:`\mathrm{Re}` and it is used to
delineate flow regimes. For circular tubes it is defined as:

.. math::

    \mathrm{Re} = \frac{\rho \langle v_{z} \rangle D}{\mu}

where :math:`\langle v_{z} \rangle` is the average flow velocity in axial
direction and :math:`D` is the tube diameter. For values up 2100 the flow is
assumed laminar if steady state is established and density is constant. For
more, see [Ref.2]_, Chapter 2.

.. autofunction:: Re


Prandtl and Schmidt Numbers
---------------------------

*Ludwig Prandtl (1875-1953) (pronounced "Prahn-t'), who taught in Hannover and
Gottingen and later served as the Director of the Kaiser Wilhelm Institute for
Fluid Dynamics, was one of the people who shaped the future of his field at the
beginning of the twentieth century; he made contributions to turbulent flow and
heat transfer, but his development of the boundary-layer equations was his
crowning achievement* [Ref.2]_. The dimensionless quantity appear under two
forms of interest for the analysis of reactors: its thermal and its chemical
versions. In thermal version, this number compares the kinematic viscosity
:math:`nu` to the thermal diffusivity :math:`\alpha`, which is replaced by
species diffusivity in its chemical version, which is more often referred to as
Schmidt number. *Ernst Heinrich Wilhelm Schmidt (1892-1975), who taught at
the universities in Gdansk, Braunschweig, and Munich (where he was the
successor to Nusselt)* [Ref.2]_. The ratio :math:`\frac{\nu}{\alpha}` indicates
the relative ease of momentum and energy or species transport in flow systems.
This dimensionless ratio in thermal form is given by

.. math::

    \mathrm{Pr} = \frac{\nu}{\alpha} = \frac{C_{p} \mu}{k}

If transport properties for a gas are not available, thermal Prandtl number can
be estimated at low pressure and non-polar molecules mixtures with help of
Eucken formula as

.. math::

    \mathrm{Pr} = \frac{C_{p}}{C_{p} + \frac{5}{4}R}

Schmidt number range can be much broader than its thermal relative, Prandtl
number. This is given by the effects of cross-section and molar weight
determining mass diffusivity of gas species. Since this number is always defined
for a pair of species, here we will try to make it more general by using the
mixture averaged diffusion coefficients as provided by Cantera's `Transport`
method `mix_diff_coeffs`.

.. math::

    \mathrm{Sc}_{min,max} = \frac{\nu}{D_{min,max}}

For more, see [Ref.2]_, Chapter 9.

.. autofunction:: Pr

.. autofunction:: Sc


Péclet Number
-------------

*Jean-Claude-Eugene Peclet (pronounced "Pay-clay" with the second syllable
accented) (1793-1857) authored several books including one on heat conduction*
[Ref.2]_. This number is nothing more than the multiplication of Reynolds and
Prandtl or Schmidt numbers. By simplifying factors one easily determines that
it represents the ratio of convective by diffusive transport (thermal or
species). High :math:`\mathrm{Pe}` limit represents the plug-flow behavior.

.. math::

    \mathrm{Pe}_{th} = \mathrm{Re} \mathrm{Pr}\qquad
    \mathrm{Pe}_{ch} = \mathrm{Re} \mathrm{Sc}

.. autofunction:: Pe


Grashof Number
--------------

*Franz Grashof (1826-1893) (pronounced "Grahss-hoff). He was professor of
applied mechanics in Karlsruhe and one of the founders of the Verein Deutscher
Ingenieure in 1856* [Ref.2]_. The Grashof number is the characteristic group
occurring in analyses of free convection. It approximates the ratio of the
buoyancy to viscous force acting on a fluid.

**TODO: implement diffusional Gr. The diffusional Grashof number arises because
of the buoyant force caused by the concentration inhomogeneities.**

.. math::

    \mathrm{Gr} = \frac{g\beta l^{3} \Delta T}{\nu^{2}}

.. autofunction:: Gr


Rayleigh Number
---------------

Rayleigh number is the multiplication of Grashof and Prandtl numbers.
See Rayleigh_.

.. _Rayleigh: https://en.wikipedia.org/wiki/Rayleigh_number

.. autofunction:: Ra
