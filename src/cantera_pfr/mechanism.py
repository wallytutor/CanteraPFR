# -*- coding: utf-8 -*-
from pathlib import Path
from typing import Sequence
import cantera as ct


MECH_XML_TEMPLATE = """\
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


class MechanismFilter:
    """ Filter kinetic mechanisms to selected species.

    Parameters
    ----------
    mech : Path | str
        Path to mechanism file in CTI, XML, or YAML format.

    Returns
    -------
    MechanismFilter
        Instance configured for filtering mechanism.
    """

    __slots__ = ('_mech',)

    def __init__(self, mech: Path | str) -> None:
        # XXX no documentation here, this is intentional
        self._mech = Path(mech)

    def filter(
        self,
        species_filter: Sequence[str] | None,
        transport: str = 'Multi',
        write: bool = False,
        idname: str | None = None,
        output: Path | str | None = None,
        overwrite: bool = False
    ) -> ct.Solution:
        """ Filter mechanism to contain only selected species.

        Parameters
        ----------
        species_filter : Sequence[str] | None
            Required species list. If None, original mechanism returned.
        transport : str = 'Multi'
            Transport model to use.
        write : bool = False
            Write filtered mechanism to XML file if True.
        idname : str | None = None
            Phase id name for output file.
        output : Path | str | None = None
            Output file path.
        overwrite : bool = False
            Allow existing output file overwrite.

        Returns
        -------
        ct.Solution
            Cantera solution object with filtered mechanism.
        """
        mech_str = str(self._mech)

        if species_filter is None:
            return ct.Solution(mech_str)

        all_spec_fn = getattr(
            ct.Species,
            'list_from_file',
            getattr(ct.Species, 'listFromFile', None)
        )
        all_reac_fn = getattr(
            ct.Reaction,
            'list_from_file',
            getattr(ct.Reaction, 'listFromFile', None)
        )

        all_spec = all_spec_fn(mech_str)
        all_reac = all_reac_fn(mech_str)
        all_elem = [k for s in all_spec for k in s.composition.keys()]

        elem = list(set(all_elem))
        spec = [s for s in all_spec if s.name in species_filter]
        reac = []

        if write:
            if idname is None:
                raise ValueError('A phase idname must be provided')

            if output is None:
                raise ValueError('An output file must be provided')

            out_path = Path(output)

            if not overwrite and out_path.exists():
                raise FileExistsError(f'Output file exists: {out_path}')

            elem_str = ' '.join(sorted(elem))
            spec_str = ' '.join(sorted(species_filter))
            content = MECH_XML_TEMPLATE.format(
                idname, elem_str, mech_str, spec_str, transport
            )
            out_path.write_text(content, encoding='utf-8')
            return ct.Solution(str(out_path))

        for r in all_reac:
            if not all(si in species_filter for si in r.reactants):
                continue

            if not all(si in species_filter for si in r.products):
                continue

            reac.append(r)

        return ct.Solution(
            thermo='IdealGas',
            kinetics='GasKinetics',
            transport_model=transport,
            species=spec,
            reactions=reac
        )


def filter_mechanism(
        mech: Path | str,
        species_filter: Sequence[str] | None,
        transport: str = 'Multi',
        write: bool = False,
        idname: str | None = None,
        output: Path | str | None = None,
        overwrite: bool = False,
        **kwargs: object
    ) -> ct.Solution:
    """ Filter kinetic mechanism to contain only selected species.

    Parameters
    ----------
    mech : Path | str
        Path to mechanism file in CTI, XML, or YAML format.
    species_filter : Sequence[str] | None
        List of required species.
    transport : str = 'Multi'
        Transport model to use.
    write : bool = False
        If True, write XML filtered mechanism.
    idname : str | None = None
        Name of new phase if write is True.
    output : Path | str | None = None
        Path to new mechanism if write is True.
    overwrite : bool = False
        Allow existing file overwrite.
    kwargs : object
        Additional unused options for compatibility.

    Returns
    -------
    ct.Solution
        Filtered Cantera solution object.
    """
    return MechanismFilter(mech).filter(
        species_filter,
        transport=transport,
        write=write,
        idname=idname,
        output=output,
        overwrite=overwrite
    )
