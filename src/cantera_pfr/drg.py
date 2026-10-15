# -*- coding: utf-8 -*-
import csv
from itertools import product
from pathlib import Path
import time
from typing import Any, Mapping, Sequence
import cantera as ct
from matplotlib import pyplot
import networkx as nx
import numpy as np
import pandas as pd
from .graph import DirectedGraph, analyse_graph
from .mechanism import filter_mechanism


class StateSpace:
    """ Represent initial state space for mechanism simplification.

    Parameters
    ----------
    t : list[float] | None = None
        List of temperatures in kelvin.
    p : list[float] | None = None
        List of pressures in pascal.
    x : list[Mapping[str, float] | str] | None = None
        List of initial gas compositions.

    Returns
    -------
    StateSpace
        Initialized state space container.
    """

    __slots__ = ('T', 'P', 'X')

    def __init__(
        self,
        t: list[float] | None = None,
        p: list[float] | None = None,
        x: list[Mapping[str, float] | str] | None = None
    ) -> None:
        # XXX no documentation here, this is intentional
        self.T = t
        self.P = p
        self.X = x


class SimplifySetup:
    """ Configuration for mechanism simplification workflow.

    Parameters
    ----------
    start_set : list[str] | None = None
        Species to retain in simplified mechanism.
    space : StateSpace | None = None
        State space conditions to scan.

    Returns
    -------
    SimplifySetup
        Initialized simplification configuration.
    """

    __slots__ = (
        'start_set',
        'plot_spec',
        'space',
        'integ',
        'nsimp',
        'threshold_min',
        'threshold_max',
        'species_filter',
        'dir_name'
    )

    def __init__(
        self,
        start_set: list[str] | None = None,
        space: StateSpace | None = None
    ) -> None:
        # XXX no documentation here, this is intentional
        self.start_set = start_set
        self.plot_spec: list[str] | None = None
        self.space = space
        self.integ: float | None = None
        self.nsimp: int | None = None
        self.threshold_min: float | None = None
        self.threshold_max: float | None = None
        self.species_filter: list[str] | None = None
        self.dir_name: Path | str | None = None

    def validate(self) -> None:
        """ Validate configuration settings.

        Returns
        -------
        None
            Returns when settings are valid.
        """
        if self.start_set is None:
            raise ValueError('Missing species start set')

        if self.threshold_max is None or self.threshold_min is None:
            raise ValueError('Missing threshold settings')

        if self.threshold_max <= self.threshold_min:
            msg = f'max {self.threshold_max} <= min {self.threshold_min}'
            raise ValueError(f'Threshold: {msg}')


class SimplifyDRG:
    """ Directed Relation Graph mechanism simplification.

    Parameters
    ----------
    graph : DirectedGraph
        Directed graph of kinetic mechanism.
    start_set : Sequence[str]
        Species to keep in simplified mechanism.

    Returns
    -------
    SimplifyDRG
        Configured simplification operator.
    """

    __slots__ = (
        '_graph',
        '_start',
        '_names',
        '_reactions',
        '_nu',
        '_nu_dict'
    )

    def __init__(
        self,
        graph: DirectedGraph,
        start_set: Sequence[str]
    ) -> None:
        # XXX no documentation here, this is intentional
        self._graph = graph
        self._start = list(set(start_set))

        nu_prods = graph.solution.product_stoich_coeffs()
        nu_reacs = graph.solution.reactant_stoich_coeffs()

        self._names = graph.solution.species_names
        self._reactions = graph.solution.reactions()
        self._nu = nu_prods - nu_reacs
        self._nu_dict = {
            s: self._nu[i, :] for i, s in enumerate(self._names)
        }

    def _update_graph(self, den_dict, rates):
        for i, (r, rt) in enumerate(zip(self._reactions, rates)):
            reac = list(r.reactants.keys())
            prod = list(r.products.keys())
            spec = list(set(reac + prod))

            def func(a, b, idx=i, rate=rt):
                return abs(self._nu_dict[a][idx] * rate)

            self._graph.update_edges(reac, spec, func)

            if r.reversible:
                prodonly = [x for x in prod if x not in reac]

                if prodonly:
                    self._graph.update_edges(prodonly, spec, func)

        self._graph.divide_edges(lambda a, b: den_dict[a], verbose=False)

    def _tag_species(self):
        h = nx.DiGraph()

        for spec in self._names:
            h.add_node(spec, mark=False, val=0.0)

            if spec in self._start:
                h.nodes[spec]['mark'] = True
                h.nodes[spec]['val'] = 1.0

        for a, b, prop in self._graph.sort_edges():
            h.add_edge(a, b, val=prop)

            if h.nodes[a]['mark'] and not h.nodes[b]['mark']:
                rab = h.edges[a, b]['val']
                tree = nx.dfs_tree(h, b)

                for node in tree.nodes():
                    if rab > h.nodes[node]['val']:
                        h.nodes[node]['val'] = rab

        sort_species = [(a, h.nodes[a]['val']) for a in h.nodes()]
        return sorted(sort_species, key=lambda p: p[1])

    def get_skeletal_mech(
        self,
        cond_list: Sequence[Any],
        threshold: float
    ) -> list[str]:
        """ Retrieve skeletal mechanism from sample states.

        Parameters
        ----------
        cond_list : Sequence[Any]
            Condition list (T, P, X, path) representing region.
        threshold : float
            Simplification threshold cutoff.

        Returns
        -------
        list[str]
            Retained species in skeletal mechanism.
        """
        keep_species: list[tuple[str, float]] = []

        for t_val, p_val, _, path in cond_list:
            data = pd.read_csv(Path(path))

            for _, row in data.iterrows():
                comp = {s: row[s] for s in self._names if s in row}
                self._graph.solution.TPX = t_val, p_val, comp
                rates = np.abs(self._graph.solution.net_rates_of_progress)
                den = np.abs(self._nu).dot(rates)
                den_dict = {s: d for s, d in zip(self._names, den)}

                self._graph.reinitialize_graph()
                self._update_graph(den_dict, rates)
                keep_species += self._tag_species()

        keep_dict: dict[str, float] = {}

        for k, v in keep_species:
            if k not in keep_dict:
                keep_dict[k] = v
            else:
                keep_dict[k] = max(v, keep_dict[k])

        retained = [
            k for k, val in keep_dict.items() if val >= threshold
        ]
        return list(set(retained))


class SimplifyDRGManager:
    """ Workflow manager for DRG mechanism simplification.

    Returns
    -------
    SimplifyDRGManager
        Manager instance.
    """

    __slots__ = ()

    @staticmethod
    def gen_name_TPXQ(
        mech: str,
        t: float,
        p: float,
        x: Mapping[str, float] | str,
        q: float = 0.0
    ) -> str:
        """ Generate standardized CSV sample file name.

        Parameters
        ----------
        mech : str
            Mechanism identifier name.
        t : float
            Temperature in kelvin.
        p : float
            Pressure in pascal.
        x : Mapping[str, float] | str
            Composition dict or label string.
        q : float = 0.0
            Inlet volume flow rate in SCCM.

        Returns
        -------
        str
            Standardized file name.
        """
        if isinstance(x, dict):
            comp_str = '_'.join(
                f'{k}_{v:.5e}' for k, v in x.items()
            )
        else:
            comp_str = str(x)

        name = f'{mech}_{q:.1f}sccm_{p:.1f}Pa_{t:.2f}K_{comp_str}'
        clean_name = name.replace('.', '_')
        return f'{clean_name}.csv'

    @staticmethod
    def gen_samples_psr(
        mech: Path | str,
        space: StateSpace,
        tend: float,
        overwrite: bool = False,
        **kwargs: object
    ) -> list[list[Any]]:
        """ Generate sample space using an isothermal PSR.

        Parameters
        ----------
        mech : Path | str
            Mechanism file path.
        space : StateSpace
            State space scan conditions.
        tend : float
            Integration end time in seconds.
        overwrite : bool = False
            Allow overwriting existing samples.
        kwargs : object
            Optional configuration arguments.

        Returns
        -------
        list[list[Any]]
            Condition tuples with file paths.
        """
        dir_name = Path(str(kwargs.get('dir_name', 'samples')))
        n_samples = int(kwargs.get('n_samples', 10))
        species = kwargs.get('species', None)
        species_seq = list(species) if species is not None else None
        flowrate = float(kwargs.get('flowrate', 0.0))

        dir_name.mkdir(parents=True, exist_ok=True)
        times = np.linspace(0, tend, n_samples)
        sol = filter_mechanism(mech, species_seq, **kwargs)
        rea = ct.IdealGasReactor(sol, energy='off')
        sim = ct.ReactorNet([rea])

        mech_name = Path(mech).name
        prod_tpx = product(
            space.T or [],
            space.P or [],
            space.X or []
        )
        head = ['T', 'P'] + sol.species_names
        cond_list: list[list[Any]] = []

        for t_val, p_val, x_val in prod_tpx:
            basename = SimplifyDRGManager.gen_name_TPXQ(
                mech_name, t_val, p_val, x_val, q=flowrate
            )
            sample_path = dir_name / basename
            cond_list.append([t_val, p_val, x_val, sample_path])

            if sample_path.exists() and not overwrite:
                continue

            sol.TPX = t_val, p_val, x_val
            rea.syncState()
            sim.set_initial_time(0.0)

            try:
                with sample_path.open('w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f)
                    writer.writerow(head)
                    writer.writerow([t_val, p_val] + list(sol.X))

                    for step_t in times[1:]:
                        sim.advance(step_t)
                        writer.writerow([t_val, p_val] + list(sol.X))
            except Exception as err:
                print(f'Error generating {sample_path}: {err}')

        return cond_list

    @staticmethod
    def manage(
        mech: Path | str,
        setup: SimplifySetup,
        overwrite: bool = False,
        **kwargs: object
    ) -> list[Path]:
        """ Manage DRG mechanism simplification workflow.

        Parameters
        ----------
        mech : Path | str
            Mechanism file path.
        setup : SimplifySetup
            Simplification configuration.
        overwrite : bool = False
            Allow overwriting existing outputs.
        kwargs : object
            Additional options.

        Returns
        -------
        list[Path]
            Generated skeletal mechanism species files.
        """
        t0 = time.time()
        setup.validate()
        base_dir = Path(setup.dir_name or 'simplification')
        base_dir.mkdir(parents=True, exist_ok=True)

        call_kwargs = dict(kwargs)
        call_kwargs['species'] = setup.species_filter
        call_kwargs['saveas'] = base_dir / 'graph' / 'analysis'
        call_kwargs['dir_name'] = base_dir / 'samples'

        dir_smp = base_dir / 'smp'
        dir_smp.mkdir(parents=True, exist_ok=True)

        dir_graph = base_dir / 'graph'
        dir_graph.mkdir(parents=True, exist_ok=True)

        if setup.space is None or setup.integ is None:
            raise ValueError('Setup space and integ cannot be None')

        cond_list = SimplifyDRGManager.gen_samples_psr(
            mech,
            setup.space,
            setup.integ,
            overwrite=overwrite,
            **call_kwargs
        )

        graph = analyse_graph(mech, DirectedGraph, **call_kwargs)

        if not isinstance(graph, DirectedGraph):
            raise TypeError('Expected DirectedGraph instance')

        if setup.start_set is None:
            raise ValueError('Setup start_set cannot be None')

        simplifier = SimplifyDRG(graph, setup.start_set)
        tmin = setup.threshold_min or 0.01
        tmax = setup.threshold_max or 0.20
        nsimp = setup.nsimp or 10
        threshold_list = np.linspace(tmin, tmax, nsimp)

        def simplify(threshold):
            smp = simplifier.get_skeletal_mech(cond_list, threshold)
            save_file = dir_smp / f'smp_{threshold:.6f}.txt'
            save_file.write_text(','.join(smp), encoding='utf-8')
            return len(smp)

        respath = dir_smp / 'residual.csv'
        residual_hist: list[list[float]] = [[], []]

        with respath.open('w', newline='', encoding='utf-8') as residuals:
            writer = csv.writer(residuals)

            for threshold in threshold_list:
                nspec = simplify(threshold)
                writer.writerow([threshold, nspec])
                residual_hist[0].append(float(threshold))
                residual_hist[1].append(float(nspec))

        figname = dir_smp / 'residual_hist.png'
        pyplot.close('all')
        pyplot.style.use('bmh')
        pyplot.plot(residual_hist[0], residual_hist[1])
        pyplot.xlabel('Simplification threshold, $\\varepsilon$')
        pyplot.ylabel('Number of residual species')
        pyplot.tight_layout()
        pyplot.savefig(figname, dpi=300)

        print(f'Simplification completed in {time.time() - t0:.2f} s')
        return sorted(dir_smp.glob('smp_*'))

    @staticmethod
    def test(
        mech: Path | str,
        setup: SimplifySetup,
        tpx: tuple[float, float, Any],
        all_smp: Sequence[Path | str],
        n_samples: int = 100,
        outfreq: int = 10,
        **kwargs: object
    ) -> None:
        """ Compare simplified mechanism solution against full mechanism.

        Parameters
        ----------
        mech : Path | str
            Mechanism file path.
        setup : SimplifySetup
            Simplification configuration.
        tpx : tuple[float, float, Any]
            Initial (T, P, X) conditions.
        all_smp : Sequence[Path | str]
            Skeletal mechanism species file paths.
        n_samples : int = 100
            Time evaluation step count.
        outfreq : int = 10
            Progress output step frequency.
        kwargs : object
            Additional options.

        Returns
        -------
        None
            Generates comparison plots and tables.
        """
        integ = setup.integ or 1.0
        times = np.linspace(0, integ, n_samples)
        base_dir = Path(setup.dir_name or 'simplification')
        outdir = base_dir / 'test'
        outdir.mkdir(parents=True, exist_ok=True)

        def run_mech(species_list, out_file):
            sol = filter_mechanism(mech, species_list, **kwargs)
            sol.TPX = tpx
            head = ['t', 'T', 'P'] + sol.species_names
            rea = ct.IdealGasReactor(sol, energy='off')
            sim = ct.ReactorNet([rea])
            sim.set_initial_time(0.0)

            with out_file.open('w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(head)
                writer.writerow([0.0, sol.T, sol.P] + list(sol.X))

                for i, step_t in enumerate(times[1:]):
                    sim.advance(step_t)
                    writer.writerow([step_t, sol.T, sol.P] + list(sol.X))

        ref_file = outdir / 'reference_solution.csv'
        run_mech(None, ref_file)
        plot_spec = setup.plot_spec or []
        data_src = pd.read_csv(ref_file, usecols=['t'] + plot_spec)

        for path in sorted(all_smp):
            smp_path = Path(path)
            basename = smp_path.stem
            smp_out = outdir / f'{basename}_smp.csv'
            species = smp_path.read_text(encoding='utf-8').strip().split(',')
            run_mech(species, smp_out)
            data_smp = pd.read_csv(smp_out, usecols=['t'] + plot_spec)

            for spec in plot_spec:
                save_fig = outdir / f'test_{spec}_{basename}.png'
                pyplot.close('all')
                pyplot.style.use('bmh')
                pyplot.plot(
                    data_smp['t'], data_smp[spec], 'o-', label='Simplified'
                )
                pyplot.plot(
                    data_src['t'], data_src[spec], label='Original'
                )
                pyplot.xlabel('Time (s)')
                pyplot.ylabel(f'Mole fraction of {spec}')
                pyplot.legend()
                pyplot.tight_layout()
                pyplot.savefig(save_fig, dpi=300)
