# -*- coding: utf-8 -*-


MECH = """\
<?xml version="1.0"?>
<ctml>
    <validate reactions="yes" species="yes"/>
    <phase dim="3" id="{0}">
        <elementArray datasrc="elements.xml">
            {1}
        </elementArray>
        <speciesArray datasrc="{2}#species_data">
            {3} <skip element="undeclared"/>
        </speciesArray>
        <reactionArray datasrc="{2}#reaction_data">
            <skip species="undeclared" third_bodies="undeclared"/>
        </reactionArray>
        <thermo model="IdealGas"/>
        <kinetics model="GasKinetics"/>
        <transport model="{4}"/>
    </phase>
</ctml>
"""


def filter_mechanism(mech, species, transport='Multi', write=False,
                     idname=None, output=None, overwrite=False):
    """ Filter mechanism to contain only selected species.

    Read mechanism and remove reactions involving species not required by
    `species` list. Transport model if specified must be `'Multi'`.

    Note
    ----
    To date transport model `'Mix'` has not being implemented.

    Parameters
    ----------
    mech : path-like
        Mechanism file in CTI/XML format.
    species : list of str
        List of required species. Case must match with mechanism.
    transport : str, optional
        Transport model to use. Default is `'Multi'`.
    write : bool, optional
        If `True`, write XML filtering mechanism.
    idname : str, optional
        Name of new phase. Required if `write=True`. Default is `None`.
    output : path-like, optional
        Path to new mechanism. Required if `write=True`. Default is `None`.
    overwrite : bool, optional
        If `True` allows mechanism to be overwritten. Default is `False`.
    """

    import cantera as ct

    all_spec = ct.Species.listFromFile(mech)
    all_reac = ct.Reaction.listFromFile(mech)
    all_elem = [k for s in all_spec for k in s.composition.keys()]

    elem = list(set(all_elem))
    spec = [s for s in all_spec if s.name in species]
    reac = []

    if write:
        assert idname is not None, 'A phase `id` must be provided'
        assert output is not None, 'An output file must be provided'
        if not overwrite:
            from os.path import exists
            assert not exists(output), 'Output fils already exists'

        elem = ' '.join(sorted(elem))
        spec = ' '.join(sorted(species))

        with open(output, 'w') as writer:
            writer.write(MECH.format(idname, elem, mech, spec, transport))

        return ct.Solution(output)

    for r in all_reac:
        if not all(si in species for si in r.reactants):
            continue
        if not all(si in species for si in r.products):
            continue
        reac.append(r)

    gas_conf = dict(thermo='IdealGas',
                    kinetics='GasKinetics',
                    transport_model=transport,
                    species=spec,
                    reactions=reac)

    return ct.Solution(**gas_conf)
