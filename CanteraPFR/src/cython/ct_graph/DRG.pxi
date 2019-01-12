# -*- coding: utf-8 -*-
from .ct_aux import filter_mechanism

import os
import csv
import glob
import time
import numpy
import cantera
import networkx
from itertools import product
from matplotlib import pyplot
from pandas import read_csv


class StateSpace(object):
    """ Represent the initial state space for mechanism simplification.

    The lists provided by this data class are provided an internal product for
    generation of all possible combinations of initial states.
    """

    T: list = None
    P: list = None
    X: list = None

    def __add__(self, other):
        """ Add generated iterator.

        Note
        ----
        This method does not add the initial lists, but the iterators created
        by this and another instance. This is useful when another region of
        sample states is desired but one does not wish to cross with the other
        conditions provided by the current instance.
        """
        pass

    def __iter__(self):
        pass


class SimplifySetup(object):
    start_set: list = None
    plot_spec: list = None
    space: StateSpace = None
    integ: float = None
    nsimp: int = None
    threshold_min: float = None
    threshold_max: float = None
    species_filter: list = None
    dir_name: str = None

    def validate(self):
        assert self.start_set is not None,\
            'Missing species start set'
        assert self.threshold_max > self.threshold_min,\
            f'Threshold: max {self.threshold_max} <= min {self.threshold_min}'


class SimplifyDRG(object):
    """ Apply DRG approach to mechanism simplification.

    Reference
    ---------
    Tianfeng Lu and Chung K. Law. On the applicability of directed relation
    graphs to the reduction of reaction mechanism. Combustion and Flame.
    vol. 146, no. 3, pp. 472-483, 2006.

    Parameters
    ----------
    graph : DirectedGraph
        Directed graph of kinetic mechanism.
    start_set: list
        List of species to keep in simplified version.
    """

    def __init__(self, graph, start_set):
        self._graph = graph
        self._start = list(set(start_set))

        nu_prods = graph.solution.product_stoich_coeffs()
        nu_reacs = graph.solution.reactant_stoich_coeffs()

        self._names = graph.solution.species_names
        self._reactions = graph.solution.reactions()

        self._nu = nu_prods - nu_reacs
        self._nu_dict = {s: n for s, n in zip(self._names, self._nu)}

    def _update_graph(self, den_dict, rates):
        """ Update relational graph edges weights.

        TODO
        ----
        This method has much to be optimized!

        Parameters
        ----------
        den_dict: dict
            Dictionary of denominators for species `abs(del_nu * omega)`.
        rates: array-like
            Reaction rates array.
        """

        for i, (r, rt) in enumerate(zip(self._reactions, rates)):
           reac = list(r.reactants.keys())
           prod = list(r.products.keys())
           spec = list(set(reac + prod))

           def func(A, B):
               return abs(self._nu_dict[A][i] * rt)

           self._graph.update_edges(reac, spec, func)

           if r.reversible:
               # TODO think about the physical meaning of this!
               prodonly = [x for x in prod if x not in reac]
               if prodonly:
                   self._graph.update_edges(prodonly, spec, func)

        # FIXME understand why so many zero denominators!
        self._graph.divide_edges(lambda A, B: den_dict[A], verbose=False)

    def _tag_species(self):
        """ Perform DFS with species tagging. """

        h = networkx.DiGraph()

        for spec in self._names:
            h.add_node(spec, mark=False, val=0.0)

            if spec in self._start:
                h.nodes[spec]['mark'] = True
                h.nodes[spec]['val'] = 1.0

        for A, B, prop in self._graph.sort_edges():
            h.add_edge(A, B, val=prop)

            #TODO Check this algorithm
            if h.nodes[A]['mark'] and not h.node[B]['mark']:
                rAB = h.edges[A, B]['val']
                tree = networkx.dfs_tree(h, B)

                for node in tree.nodes():
                    if rAB > h.nodes[node]['val']: # not in original
                        h.nodes[node]['val'] = rAB

        sort_species = [(A, h.nodes[A]['val']) for A in h.nodes()]
        return sorted(sort_species, key=lambda p: p[1])

    def get_skeletal_mech(self, cond_list, threshold):
        """ Retrieve skeletal mechanism from samples state space.

        This routine produces skeletal mechanisms that are robust within a state
        space represented by points in the files listed in the conditions list.

        TODO
        ----
        This method has much to be optimized!

        Parameters
        ----------
        cond_list : list
            List of files containing state space points representing the region
            where integration of simplified mechanism is expected to be of good
            precision. The format of this list must be the one returned by
            :func:`gen_samples_psr`.
        threshold : float
            Simplification threshold. See references.

        Returns
        -------
        list
            List of species retained with use of given threshold.
        """

        keep_species = []
        for T, P, X, path in cond_list:
            print(f'\n Simplifying from state {path}')
            data = read_csv(path)

            for nrow, row in data.iterrows():
                print(f' Currently on row {nrow}')

                self._graph.solution.TPX = T, P, dict(row[self._names])
                rates = numpy.abs(self._graph.solution.net_rates_of_progress)

                den = numpy.abs(self._nu).dot(rates)
                den_dict = {s: d for s, d in zip(self._names, den)}

                self._graph.reinitialize_graph()
                self._update_graph(den_dict, rates)
                keep_species += self._tag_species()

        keep_species = sorted(keep_species, key=lambda x: x[0])
        keep_dict = {}

        for k, v in keep_species:
            if k not in keep_dict:
                keep_dict[k] = v
            else:
                keep_dict[k] = max(v, keep_dict[k])

        keep_species = []
        for key in keep_dict.keys():
            if keep_dict[key] >= threshold:
                keep_species.append(key)

        return list(set(keep_species))


class SimplifyDRGManager(object):
    """ Provides a workflow for DRG simplification.

    """

    @staticmethod
    def gen_name_TPXQ(mech, T, P, X, Q=0.0):
        """ Generate a standard file name.

        This function is intended to generate a CSV file name in terms of the
        operating conditions of a gas phase perfect stirred reactor. The
        parameters taken into account comprise the pressure, temperature, inlet
        composition, and flow rate.

        Parameters
        ----------
        mech : str
            Mechanism name or another identifier.
        T : float
            Characteristic temperature in kelvin.
        P : float
            Operating pressure in pascal.
        X : dict or str
            Dictionary of species and respective mole fractions. If many species
            are used in intial state, use a composition identifier instead.
        Q : float, optional
            Inlet volume flow rate in SCCM. Default is 0.0.

        Returns
        -------
        str
            Standardized CSV file name.
        """

        if isinstance(X, dict):
            X = '_'.join(map(lambda a: f'{a[0]:s}_{a[1]:.5e}', X.items()))

        name = f'{mech}_{Q:.1f}sccm_{P:.1f}Pa_{T:.2f}K_{X}'
        return f'{name.replace(".", "_")}.csv'

    @staticmethod
    def gen_samples_psr(mech, space, tend, overwrite=False, **kwargs):
        """ Generate sample space using an isothermal PSR.

        TODO
        ----
        Replace CSV writer by a cantera writer or a DataFrame.
        Consider generating a single file instead.

        Parameters
        ----------
        mech : str
            Mechanism file or path.
        space : StateSpace
            Object describing set of conditions to scan.
        tend : float
            Integration time to sample from.
        overwrite : bool, optional
            If `True` allows samples to be overwritten. Default is `True`.
            FIXME cannot have overwrite because of kwargs!!
        Returns
        -------
        list
            List of lists with conditions and corresponding file name.
        """

        dir_name = kwargs.get('dir_name', 'samples')
        n_samples = kwargs.get('n_samples', 10)
        species = kwargs.get('species', None)
        flowrate = kwargs.get('flowrate', 0.0)

        if not os.path.exists(dir_name):
            os.makedirs(dir_name)

        times = numpy.linspace(0, tend, n_samples)
        sol = filter_mechanism(mech, species, **kwargs)
        rea = cantera.IdealGasReactor(sol, energy='off')
        sim = cantera.ReactorNet([rea])

        mech_name = os.path.basename(mech)
        prod_TPX = product(space.T, space.P, space.X)
        head = ['T', 'P'] + sol.species_names
        cond_list = []

        for T, P, X in prod_TPX:
            basename = SimplifyDRGManager.gen_name_TPXQ(mech_name, T, P, X,
                                                        Q=flowrate)
            basename = os.path.join(dir_name, basename)
            cond_list.append([T, P, X, basename])

            print(f'\n Working on space:\n {basename}\n')
            if os.path.exists(basename) and not overwrite:
                print(' Sample already exists, skipping.')
                continue

            sol.TPX = T, P, X
            rea.syncState()
            sim.set_initial_time(0.0)

            try:
                with open(basename, 'w') as outfile:
                    writer = csv.writer(outfile)
                    writer.writerow(head)
                    writer.writerow([T, P] + list(sol.X))

                    for t in times[1:]:
                        print(f' Step at {t:.3e} s')
                        sim.advance(t)
                        writer.writerow([T, P] + list(sol.X))
            except (Exception) as err:
                print(f'\n {err}:\n while generating:\n {basename}')

        return cond_list

    @staticmethod
    def manage(mech, setup, overwrite=False, **kwargs):
        """ Manage DRG mechanism simplification.

        TODO
        ----
        Save `setup` to a dictionary to be able to check if existing
        simplifications are not being repeated and thus avoid running again
        during post-processing.

        Parameters
        ----------
        mech : str
            Mechanism file or path.
        setup : SimplifySetup
            Structure containing simplification setup.
        """

        t0 = time.time()
        setup.validate()

        if not os.path.exists(setup.dir_name):
            os.makedirs(setup.dir_name)

        if 'species' in kwargs:
            print('Overwritting species list from keyword arguments')

        if 'saveas' in kwargs:
            print('Overwritting save name from keyword arguments')

        if 'dir_name' in kwargs:
            print('Overwritting directory name from keyword arguments')

        # TODO kwargs parsed to filter_mechanism by analyse_graph.
        #transport='Multi', write=False,
        #idname=None, output=None, overwrite=False)

        kwargs['species'] = setup.species_filter
        kwargs['saveas'] = os.path.join(setup.dir_name, 'graph', 'analysis')
        kwargs['dir_name'] = os.path.join(setup.dir_name, 'samples')

        dir_smp = os.path.join(setup.dir_name, 'smp')
        if not os.path.exists(dir_smp):
            os.makedirs(dir_smp)

        # FIXME this should be in analyse_graph
        dir_graph = os.path.join(setup.dir_name, 'graph')
        if not os.path.exists(dir_graph):
            os.makedirs(dir_graph)

        cond_list = SimplifyDRGManager.gen_samples_psr(mech, setup.space,
                                                       setup.integ,
                                                       overwrite=overwrite,
                                                       **kwargs)

        graph = analyse_graph(mech, DirectedGraph, **kwargs)
        simplifier = SimplifyDRG(graph, setup.start_set)

        tmin = setup.threshold_min
        tmax = setup.threshold_max
        threshold_list = numpy.linspace(tmin, tmax, setup.nsimp)

        def simplify(threshold):
            smp = simplifier.get_skeletal_mech(cond_list, threshold)
            saveas = os.path.join(dir_smp, f'smp_{threshold:.6f}.txt')

            # TODO transform all this in JSON dict!
            with open(saveas, 'w') as writer:
                writer.write(','.join(smp))

            return len(smp)

        # TODO it would be much faster if threshold loop was done inside the
        # simplification class, but new data structures need to be developped.
        respath = os.path.join(dir_smp, 'residual.csv')
        residual_hist = [[], []]
        with open(respath, 'w') as residuals:
            writer = csv.writer(residuals)

            for threshold in threshold_list:
                print(f'\n Working with threshold {threshold:.4f}')
                nspec = simplify(threshold)
                writer.writerow([threshold, nspec])
                residual_hist[0].append(threshold)
                residual_hist[1].append(nspec)

        figname = os.path.join(dir_smp, 'residual_hist.png')

        pyplot.close('all')
        pyplot.style.use('bmh')
        pyplot.plot(residual_hist[0], residual_hist[1])
        pyplot.xlabel('Simplication threshold, $\\varepsilon$ ')
        pyplot.ylabel('Number of residual species')
        pyplot.tight_layout()
        pyplot.savefig(figname, dpi=300)

        print(f'Simplification took {time.time()-t0} s')
        return glob.glob(os.path.join(dir_smp, r'smp_*'))

    @staticmethod
    def test(mech, setup, TPX, all_smp, n_samples=100, outfreq=10, **kwargs):
        """ Compares solution of simplified mechanism to its full counterpart.

        Parameters
        ----------
        """

        times = numpy.linspace(0, setup.integ, n_samples)
        outdir = os.path.join(setup.dir_name, 'test')
        if not os.path.exists(outdir):
            os.makedirs(outdir)

        def run_mech(species, basename):
            sol = filter_mechanism(mech, species, **kwargs)
            sol.TPX = TPX

            head = ['t', 'T', 'P'] + sol.species_names
            rea = cantera.IdealGasReactor(sol, energy='off')
            sim = cantera.ReactorNet([rea])
            sim.set_initial_time(0.0)

            try:
                with open(basename, 'w') as outfile:
                    writer = csv.writer(outfile)
                    writer.writerow(head)
                    writer.writerow([0, sol.T, sol.P] + list(sol.X))

                    for i, t in enumerate(times[1:]):
                        if not i % outfreq:
                            print(f' Step at {t:.3e} s')
                        sim.advance(t)
                        writer.writerow([t, sol.T, sol.P] + list(sol.X))

            except (Exception) as err:
                print(f'\n {err}:\n while generating:\n {basename}')

        basename_src = os.path.join(outdir, f'reference_solution.csv')
        run_mech(None, basename_src)
        data_src = read_csv(basename_src, usecols=['t'] + setup.plot_spec)

        for path in sorted(all_smp):
            print(f'\n Testing {path}')
            basename = '_'.join(os.path.basename(path).split('.')[:-1])
            basename_smp = os.path.join(outdir, f'{basename}_smp.csv')

            # TODO fix this when using JSON!
            with open(path, 'r') as reader:
                species = reader.read().split(',')

            run_mech(species, basename_smp)
            data_smp = read_csv(basename_smp, usecols=['t'] + setup.plot_spec)

            for spec in setup.plot_spec:
                saveas = os.path.join(outdir, f'test_{spec}_{basename}.png')

                pyplot.close('all')
                pyplot.style.use('bmh')
                pyplot.plot(data_smp['t'], data_smp[spec], 'o-', label='Simplified')
                pyplot.plot(data_src['t'], data_src[spec], label='Original')
                pyplot.xlabel('Time (s)')
                pyplot.ylabel(f'Mole fraction of {spec}')
                pyplot.legend()
                pyplot.tight_layout()
                pyplot.savefig(saveas, dpi=300)
