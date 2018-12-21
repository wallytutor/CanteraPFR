AdiabaticPFR Plug-Flow Reactor Model
====================================

Model limitations/hypothesis:


Functionalities to implement:


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

Module documentation
--------------------

.. autoclass:: CanteraPFR.models._hom_adiabatic.AdiabaticPFR
   :members:
