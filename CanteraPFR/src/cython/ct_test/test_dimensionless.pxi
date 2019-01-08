# -*- coding: utf-8 -*-


def test_basic():
    import cantera as ct
    from .ct_aux import Re, Pr, Sc, Pe, Gr, Ra

    print('\n\nTest basic functioning')

    Uz = 2.00
    Di = 0.05
    Tw = 1173.0

    for trans in ['mix', 'multi']:
        print(f'\nUsing transport model {trans}')
        gas = ct.Solution('gri30.cti', f'gri30_{trans}', transport=trans)
        gas.TPX = 298.15, 5000, "N2:0.64, C2H2:0.36"

        print(f'Re = {Re(gas, Uz, Di)}')
        print(f'Pr = {Pr(gas)}')
        print(f'Sc = {Sc(gas)}')
        print(f'Pe = {Pe(gas, Uz, Di)}')
        print(f'Gr = {Gr(gas, Di, Tw)}')
        print(f'Ra = {Ra(gas, Di, Tw)}')


def test_fallback():
    import cantera as ct
    from .ct_aux import Re, Pr, Sc, Pe, Gr, Ra

    print('\n\nTest function fallbacks')

    Uz = 2.00
    Di = 0.05
    Tw = 1173.0
    mu = 1.47e-05
    k = 0.025
    Dab = 0.00021

    gas = ct.Solution('CT-hydrocarbon-norinaga-2009-mech.cti')
    gas.TPX = 298.15, 5000, "N2:0.64, C2H2:0.36"

    print('\nTest Reynolds')
    print(f'Re = {Re(gas, Uz, Di, mu=mu)}')
    try:
        print(f'Re = {Re(gas, Uz, Di)}')
    except Exception as err:
        print(err)

    print('\nTest Prandtl')
    print(f'Pr = {Pr(gas, mu=mu, k=k)}')
    try:
        print(f'Pr = {Pr(gas, mu=mu)}')
    except Exception as err:
        print(err)

    try:
        print(f'Pr = {Pr(gas, k=k)}')
    except Exception as err:
        print(err)

    try:
        print(f'Pr = {Pr(gas)}')
    except Exception as err:
        print(err)

    print('\nTest Schmidt')
    print(f'Sc = {Sc(gas, mu=mu, Dab=Dab)}')
    try:
        print(f'Sc = {Sc(gas, mu=mu)}')
    except Exception as err:
        print(err)
    try:
        print(f'Sc = {Sc(gas, Dab=Dab)}')
    except Exception as err:
        print(err)
    try:
        print(f'Sc = {Sc(gas)}')
    except Exception as err:
        print(err)

    print('\nTest Peclet')
    print(f'Pe = {Pe(gas, Uz, Di, mu=mu, Dab=Dab, k=k)}')
    try:
        print(f'Pe = {Pe(gas, Uz, Di, Dab=Dab, k=k)}')
    except Exception as err:
        print(err)
    try:
        print(f'Pe = {Pe(gas, Uz, Di, mu=mu, k=k)}')
    except Exception as err:
        print(err)
    try:
        print(f'Pe = {Pe(gas, Uz, Di, mu=mu, Dab=Dab)}')
    except Exception as err:
        print(err)

    print('\nTest Grashof')
    print(f'Gr = {Gr(gas, Di, Tw, mu=mu)}')
    try:
        print(f'Gr = {Gr(gas, Di, Tw)}')
    except Exception as err:
        print(err)

    print('\nTest Rayleigh')
    print(f'Ra = {Ra(gas, Di, Tw, mu=mu, k=k)}')
    try:
        print(f'Ra = {Ra(gas, Di, Tw, mu=mu)}')
    except Exception as err:
        print(err)
    try:
        print(f'Ra = {Ra(gas, Di, Tw, k=k)}')
    except Exception as err:
        print(err)


def test_dimensionless():
    test_basic()
    test_fallback()
