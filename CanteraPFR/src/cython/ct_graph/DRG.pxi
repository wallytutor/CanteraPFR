# -*- coding: utf-8 -*-
from .ct_aux import plot_adjmatrix
from .ct_aux import plot_deghist
from .ct_aux import filter_mechanism

import os
import itertools
import numpy
import networkx
from itertools import combinations
from cantera import Solution

class BaseGraph(networkx.Graph):
    """ Base class for chemical graphs.

    This class should be initialized in derived classes by calling method
    `_base_init` after inheritance initilization. Its interface is described
    here instead of the `__init__` one.

    Parameters
    ----------
    sol : cantera.Solution
        Solution object representing the kinetics mechanism.
    fullfmt : bool
        If `True`, print all edges in mechanism
    """

    def __init__(self):
        super(BaseGraph, self).__init__()

    def _base_init(self, sol, fullfmt):
        self._sol = sol
        self._fullfmt = fullfmt
        self.add_nodes_from(self._sol.species_names)

    def __str__(self):
        """ Formatted string representation. """

        basefmt = (f' {self.__class__.__name__} object\n'
                   f' -- Size   : {self.size()}\n'
                   f' -- Order  : {self.order()}\n')

        def more(edges):
            more = map(lambda e: f' {e[0]:>20s} - {e[1]:s}\n', edges)
            # extra = [' {:>20s} - {:s}\n'.format(*e) for e in self.edges()]
            return ''.join(more)

        return (basefmt + '' if not self._fullfmt else more(self.edges()))

    def adjacency_matrix(self):
        """ Returns the adjacency matrix.

        Returns
        -------
        numpy.array
            The adjacency matrix of interacting nodes.
        """

        return networkx.to_numpy_matrix(self)

    def degree_histogram(self):
        """ Return of nodes degrees.

        Returns
        -------
        list
            List indexed by degree with the number of its occurences.
        """

        return networkx.degree_histogram(self)

    def ordered_degree(self):
        """ Returns of nodes ordered by degree.

        Returns
        -------
        list
            List of pairs <node, degree>, ordered by degree.
        """

        baselist = map(lambda v: (v, self.degree(v)), self)
        # baselist = [(v, self.degree(v)) for v in self]
        return sorted(baselist, key=lambda x: x[1], reverse=True)

    @property
    def density(self):
        """ Graph density. """

        return networkx.density(self)

    def plot_adjmatrix(self, saveas, overwrite=True, **kwargs):
        """ Interface to :meth:`CanteraPFR.ct_aux.plot_adjmatrix`. """

        plot_adjmatrix(self, saveas, overwrite=overwrite, **kwargs)

    def plot_deghist(self, saveas, overwrite=True, **kwargs):
        """ Interface to :meth:`CanteraPFR.ct_aux.plot_deghist`. """

        plot_deghist(self, saveas, overwrite=overwrite, **kwargs)


class DirectedGraph(networkx.DiGraph, BaseGraph):
    """ Generate species directed relational graph.

    An edge is added for every pair of interacting species in a reaction. If
    :math:`A\\rightarrow{}B`, then a directed edge is added between these.

    Note
    ----
        No self loops :math:`(A\\rightarrow{}A)` are generated in the graph.

    Parameters
    ----------
    sol : cantera.Solution
        Solution object representing the kinetics mechanism.
    fullfmt : bool, optional
        If `True`, print all edges in mechanism. Default is `False`.
    """

    def __init__(self, sol, fullfmt=False):
        # networkx.DiGraph.__init__(self)
        # BaseGraph.__init__(self, sol, fullfmt)
        super(DirectedGraph, self).__init__()
        self._base_init(sol, fullfmt)

        for r in self._sol.reactions():
            reac = list(r.reactants.keys())
            prod = list(r.products.keys())
            spec = reac + prod

            self._add_edges(reac, spec)
            if r.reversible:
                self._add_edges(prod, spec)

        self = networkx.freeze(self)

    def _add_edges(self, src, dst):
        """ Add edges from src to dst without self loops on nodes. """

        # TODO improve with itertools (or other).
        src = numpy.array(src)
        dst = numpy.array(dst)

        for B in src:
            for A in dst[dst != B]:
                self.add_edge(A, B, weight=0.0)

    def update_edges(self, src, dst, func):
        """ Update edges values.

        Parameters
        ----------
        src : list
            List of source nodes names.
        dst : list
            List of destination nodes names.
        func : function
            Function used to increment edges `weight` property. This function
            takes the names of nodes `A` and `B` as arguments.
        """

        # TODO initialize edges visit here and remove reinitialize_graph.
        # TODO improve with itertools (or other).
        src = numpy.array(src)
        dst = numpy.array(dst)

        for B in src:
            for A in dst[dst != B]:
                self[A][B]['weight'] += func(A, B)

    def divide_edges(self, func):
        """ Divide edges by function.

        Parameters
        ----------
        func : function
            Function used to increment edges `weight` property. This function
            takes the names of nodes `A` and `B` as arguments.
        """

        for A, B in networkx.edges_iter(self):
            val = func(A, B)

            if val == 0.0:
                print(f'Warning {A}->{B} has zero denominator. Setting 1.')
                self[A][B]['val'] = 1.0
                continue

            ratio = self[A][B]['val'] / val
            self[A][B]['val'] = min(ratio, 1.0)

    def sort_edges(self, func, w='weight', reverse=True):
        """ Sort graph edges by key using a function.

        Parameters
        ----------
        func : function
            Function to provide ordering.
        reverse : bool, optional
            If `True`, returns in reverse order.

        Returns
        -------
        list
            List of edges ordered accoding to `func`.
        """

        def mapper(A, B):
            return (A, B, self[A][B][w])

        # Check if sorted accepts map directly.
        # baselist = [(A, B, {w: self[A][B][w]}) for A, B in networkx.edges_iter(self)]
        baselist = list(map(mapper, networkx.edges_iter(self)))
        return sorted(baselist, key=func, reverse=reverse)

    def ordered_in_degree(self):
        """ Returns list of tuples (node, degree) ordered by degree.

        Returns
        -------
        list
            List of tuples of species and respective degree.
        """

        # baselist = [(v, self.in_degree(v)) for v in self]
        baselist = map(lambda v: (v, self.in_degree(v)), self)
        return sorted(baselist, key=lambda x: x[1], reverse=True)

    def ordered_out_degree(self):
        """ Returns list of tuples (node, degree) ordered by degree.

        Returns
        -------
        list
            List of tuples of species and respective degree.
        """

        # baselist = [(v, self.out_degree(v)) for v in self]
        baselist = map(lambda v: (v, self.out_degree(v)), self)
        return sorted(baselist, key=lambda x: x[1], reverse=True)

    def reinitialize_graph(self):
        """ Set the value of property `weight` of all edges to zero. """

        for A, B in networkx.edges_iter(self):
            self[A][B]['edges'] = 0.0

    def complexes(self):
        """ Return complexes for building Feinberg–Horn–Jackson graph. """

        raise NotImplementedError('Not yet implemented')

    def plot_adjmatrix(self, saveas, overwrite=True, **kwargs):
        print('plot_adjmatrix is not available for DirectedGraph')


class UndirectedGraph(BaseGraph):
    """ Generate species undirected relational graph.

    An edge is added for every pair of species in a reaction.

    Note
    ----
        No self loops :math:`(A\\rightarrow{}A)` are generated in the graph.

    Parameters
    ----------
    sol : cantera.Solution
        Solution object representing the kinetics mechanism.
    fullfmt : bool, optional
        If `True`, print all edges in mechanism. Default is `False`.
    """

    def __init__(self, sol, fullfmt=True):
        super(UndirectedGraph, self).__init__()
        self._base_init(sol, fullfmt)

        for r in self._sol.reactions():
            reac = r.reactants.keys()
            prod = r.products.keys()
            spec = list(reac) + list(prod)

            # for p, q in list(combinations(spec, 2)):
            for p, q in combinations(spec, 2):
                if (p, q) not in self and (p != q):
                    self.add_edge(p, q)

        self = networkx.freeze(self)


def analyse_graph(mech, graph_class, **kwargs):
    """ Automatic analysis of species graph.

    This method provides a preliminary mechanism species graph analysis.
    It is able to filter species and reactions from a species file and
    perform analysis of simpler versions of a mechanism.  An input file
    and graph type are required inputs.  Optional arguments are the base
    name to save the resulting files and the list of species to filter
    the mechanism.  Graph species connection degree output is supplied
    as well as a basic report with the number of nodes, edges and density.
    Plots adjacency matrix graph histogram.

    TODO
    ----
    Document keyword arguments.

    Parameters
    ----------
    mech: str
        Mechanism file path in Cantera's format.
    graph_class: BaseGraph
        Graph class to be used.

    Returns
    -------
    BaseGraph
        Return generated graph with derived type `graph_class`.
    """

    print(f' Starting {graph_class.__name__} analysis approach')

    saveas = kwargs.get('saveas', f'results_{graph_class.__name__}')
    species = kwargs.get('species', None)

    sol = filter_mechanism(mech, species, **kwargs)
    mgraph = graph_class(sol)

    if isinstance(mgraph, DirectedGraph):
        with open(f'{saveas}-inward.txt', 'w') as writer:
            for v, d in mgraph.ordered_in_degree():
                writer.write(f'{v:16s} {d:4d}\n')

        with open(f'{saveas}-outward.txt', 'w') as writer:
            for v, d in mgraph.ordered_out_degree():
                writer.write(f'{v:16s} {d:4d}\n')

    with open(f'{saveas}-global.txt', 'w') as writer:
        for v, d in mgraph.ordered_degree():
            writer.write(f'{v:16s} {d:4d}\n')

    with open(f'{saveas}-report.txt', 'w') as writer:
        basestr = (f'Species graph contains {mgraph.size()} edges\n'
                   f'Species graph contains {mgraph.order()} nodes\n'
                   f'Species graph density  {mgraph.density}\n')
        writer.write(basestr)

    mgraph.plot_deghist(f'{saveas}-histogram.png')
    mgraph.plot_adjmatrix(f'{saveas}-adjacency.png')

    return mgraph






#         def gen_nameTPXQ(mech, T, P, X, Q=0.0):
#             """ Generate a standard file name.
#
#             This function is intended to generate a csv file name in terms of the
#             operating conditions of a gas phase chemical reactor.  The parameters
#             taken into account comprise the pressure, temperature, flow rate and
#             inlet composition.
#
#             Parameters
#             ----------
#             mech : :obj:`str`
#                 Mechanism name.
#             T : :obj:`float`
#                 Characteristic temperature in :math:`[K]`.
#             P : :obj:`float`
#                 Operating pressure in :math:`[Pa]`.
#             X : :obj:`dict` of :obj:`str`: :obj:`float`
#                 Dictionary of species and respective mole fractions.
#             Q : :obj:`float`, optional
#                 Inlet volume flow rate in :math:`[cm^{3}\\,min^{-1}]`.
#                 Default is 0.0.
#
#             Returns
#             -------
#             :obj:`str`
#                 Standard csv file name.
#             """
#             nQPT = '{:04d}sccm-{:06d}Pa-{:04d}K'.format(int(Q), int(P), int(T))
#             X = ['{:s}_{:.6e}'.format(k, v) for k, v in X.items()]
#             X = '-'.join(sorted([x.replace('.', '_') for x in X]))
#             return '{}-{}-{}.csv'.format(mech.split('.')[0], nQPT, X)
#
#
#         def generate_samples(mechanism, space, tend, **kwargs):
#             """ Generate sample space from IdealGasConsPressureReactor
#                 simulations carried with energy equation off.
#
#                 mechanism: cti mechanism file on cantera's path.
#                 space: space state to scan
#
#                 returns: list of lists of conditions used and files generated.
#             """
#             saveat = kwargs.get('saveat', 'generated-samples')
#             samples = kwargs.get('samples', 10)
#             species_list = kwargs.get('species_list', None)
#
#             if not os.path.exists(saveat):
#                 os.mkdir(saveat)
#
#             sol = mechanism_from_species(mechanism, species_list)
#             T, P, X = space['T'], space['P'], space['X']
#             rea = ct.IdealGasReactor(sol, energy='off')
#             sim = ct.ReactorNet([rea])
#             condlist = []
#
#             for TPX in list(itertools.product(T, P, X)):
#                 T, P, X = TPX
#                 sol.TPX = T, P, X
#                 rea.syncState()
#                 sim.set_initial_time(0.0)
#                 basename = gen_nameTPXQ(mechanism, P, T, X, Q=0.0)
#                 basename = os.path.join(saveat, 'sample-{}'.format(basename))
#                 condlist.append([T, P, basename])
#
#                 print('\n  Working on space:\n  %s\n' % basename)
#                 if os.path.isfile(basename):
#                     print('  Sample already exists, skipping.')
#                     continue
#
#                 try:
#                     outfile = open(basename, 'w')
#                     writer = csv.writer(outfile)
#                     writer.writerow(sol.species_names)
#
#                     dt = t/samples
#                     for i in range(1, samples+1):
#                         print('  Step at %.3e s' %(i*dt))
#                         sim.advance(i*dt)
#                         writer.writerow(sol.X)
#                 except IOError:
#                     print('\n  Cannot write results:\n  %s\n' % basename)
#                     pass
#                 finally:
#                     outfile.close()
#             return condlist
#
#         #class lu_law_simplify:
#         #    def __init__(self, digraph, start_set):
#         #        """ Reduce kinect mechanism contained in directed graph.
#         #
#         #            digraph: the mechanism species directed graph.
#         #            start_set: species to keep in mechanism.
#         #        """
#         #        self.sol = digraph.sol
#         #        self.gph = digraph
#         #        self.set = start_set
#         #
#         #        self.names = self.sol.species_names
#         #
#         #        nu_prods = self.sol.product_stoich_coeffs()
#         #        nu_reacs = self.sol.reactant_stoich_coeffs()
#         #        self.nu = np.array(nu_prods-nu_reacs)
#         #        self.nu_dict = {spec: self.nu[i] for i,spec in enumerate(self.names)}
#         #
#         #        self.reactions = self.sol.reactions()
#         #
#         #
#         #    def get_skeletal_mechanism(self,  cond_list, threshold=0.2):
#         #        """ For given conditions list, return a list of species to
#         #            retain in simplified mechanism.
#         #
#         #            cond_list: list of conditions (T, P, sample file).
#         #            threshold: cut-off value for mechanism reduction.
#         #
#         #            returns: set of species to keep in reduced mechanism.
#         #        """
#         #        keep_species = []
#         #        for T,P,f in cond_list:
#         #            print('\n  Simplifying from state %s' % f)
#         #            data = pd.read_csv(f, header=0, names=self.names)
#         #
#         #            for nrow,row in enumerate(data.iterrows()):
#         #                print('  Currently on row %d' % nrow)
#         #
#         #                X = [dict(x) for x in row[1:]][0] #TODO take care, workaround
#         #                self.sol.TPX = T, P, X
#         #
#         #                rates = np.abs(self.sol.net_rates_of_progress)
#         #
#         #                den = np.abs(self.nu).dot(rates)
#         #                den_dict = {spec: den[i] for i,spec in enumerate(self.names)}
#         #
#         #                self.gph.reinitialize_graph()
#         #                self._updategraph(den_dict, rates)
#         #
#         #                func = lambda x: self.gph.digraph[x[0]][x[1]]['val']
#         #                sort_ratio_g = self.gph.sort_edges(func)
#         #
#         #                h = nx.DiGraph()
#         #
#         #                for spec in self.names:
#         #                    h.add_node(spec, {'mark': False, 'val':0.0})
#         #
#         #                    if spec in self.set:
#         #                        h.node[spec]['mark'] = True
#         #                        h.node[spec]['val'] = 1.0
#         #
#         #                for t in sort_ratio_g:
#         #                    A, B, prop = t
#         #                    h.add_edge(A, B, prop)
#         #
#         #                    #TODO Check this algorithm
#         #                    if h.node[A]['mark'] and not h.node[B]['mark']:
#         #                        rAB = h[A][B]['val']
#         #                        tree = nx.dfs_tree(h, B)
#         #
#         #                        for node in nx.nodes_iter(tree):
#         #                            if rAB > h.node[node]['val']: # not in original
#         #                                h.node[node]['val'] = rAB
#         #
#         #                sort_species = [(A, h.node[A]['val']) for A in nx.nodes_iter(h)]
#         #                sort_species = sorted(sort_species, key=lambda p: p[1])
#         #                keep_species += sort_species
#         #
#         #        keep_species = sorted(keep_species, key=lambda x: x[0])
#         #        keep_dict = {}
#         #        for k,v in keep_species:
#         #            if k not in keep_dict:
#         #                keep_dict[k] = v
#         #            else:
#         #                keep_dict[k] = max(v, keep_dict[k])
#         #
#         #        keep_species = []
#         #        for key in keep_dict.keys():
#         #            if keep_dict[key] >= threshold:
#         #                keep_species.append(key)
#         #
#         #        return list(set(keep_species))
#         #
#         #
#         #    def _updategraph(self, den_dict, rates):
#         #        """ Update edges values from rates.
#         #
#         #            den_dict: dictionary of denominators for species abs(del_nu*omega).
#         #            rates: reaction rates omega.
#         #        """
#         #        for i,(r,rt) in enumerate(zip(self.reactions,rates)):
#         #            reac = list(r.reactants.keys())
#         #            prod = list(r.products.keys())
#         #            spec = list(set(reac+prod))
#         #
#         #            func = lambda A, B: abs(self.nu_dict[A][i]*rt)
#         #            self.gph.update_edges(reac, spec, func)
#         #
#         #            if r.reversible:
#         #                prodonly = [x for x in prod if x not in reac] #TODO think about
#         #                if prodonly:
#         #                    self.gph.update_edges(prodonly, spec, func)
#         #
#         #        func = lambda A, B: den_dict[A]
#         #        self.gph.divide_nodes(func)
#         #
#         ## ****************************************************************************
#         ## function simplify_mechanism
#         ##TODO allow to write cti file from here!
#         ## ****************************************************************************
#         #
#         #def simplify_mechanism(mechanism, setup, **kwargs):
#         #    """ Use reduce class to produce a simplified reaction mechanism.
#         #
#         #        mechanism: name o cti file containing mechanism to reduce.
#         #        digraph: a directed graph of the mechanism.
#         #        cond_list: conditions to explore in reduction.
#         #        start_set: species to keep anyways in mechanism. TODO check!
#         #        threshold_list: list of cut-off parameters to explore reduction.
#         #        save: existing directory to save results.
#         #
#         #    """
#         #    # Parse kwargs
#         #    saveat    = xu._get_kwarg('saveat',  kwargs, 'simplified-mechanism')
#         #    thr_max   = xu._get_kwarg('thr_max', kwargs, 0.2)
#         #    thr_num   = xu._get_kwarg('thr_max', kwargs, 20)
#         #    gonly     = xu._get_kwarg('gonly',     kwargs, None)
#         #    speciesfi = xu._get_kwarg('speciesfi', kwargs, None)
#         #
#         #    # Deal with folder
#         #    xu._get_outfolder(saveat)
#         #
#         #    # Deal with setup
#         #    msg = '\n  Cannot continue without `start_set`'
#         #    start_set = setup['start_set'] if 'start_set' in setup else sys.exit(msg)
#         #
#         #    msg = '\n  Cannot continue without `cond_list`'
#         #    cond_list = setup['cond_list'] if 'cond_list' in setup else sys.exit(msg)
#         #
#         #    thr_max = setup['thr_max'] if 'thr_max' in setup else thr_max
#         #    thr_num = setup['thr_num'] if 'thr_num' in setup else thr_num
#         #
#         #    if 'digraph' in setup:
#         #        digraph = setup['digraph']
#         #    else:
#         #        analysis_setup = {'gonly': gonly, 'speciesfi': speciesfi}
#         #        digraph = analyse_graph(mechanism, 'digraph', **analysis_setup)
#         #
#         #    # Create simplification manager
#         #    name = digraph.name
#         #    red = lu_law_simplify(digraph, start_set)
#         #
#         #    # Auxiliary private function
#         #    def _simp(thr, fname):
#         #        smp = red.get_skeletal_mechanism(cond_list, threshold=thr)
#         #
#         #        try:
#         #            specfile = open(saveas, 'w')
#         #            specfile.write(' '.join([str(s)+' 'for s in smp]))
#         #        except IOError:
#         #            print('\n  Cannot write %s' % fname)
#         #            raise SystemExit
#         #        finally:
#         #            specfile.close()
#         #
#         #        return len(smp)
#         #
#         #    # Loop over threshold list
#         #    # TODO it would be much faster if threshold loop was done inside the
#         #    # simplification class, but new data structures need to be developped.
#         #    resname = name+'-residual.csv'
#         #    resname = os.path.join(saveat, resname)
#         #
#         #    if not os.path.isfile(resname):
#         #        try:
#         #            outfile = open(resname, 'w')
#         #            writer = csv.writer(outfile)
#         #
#         #            for thr in np.linspace(0.0, thr_max, thr_num):
#         #                print('\n  Working with threshold %.4f' % thr)
#         #                fname = name+'-threshold-'+str('%.5f' % thr).replace('.','_')+'.txt'
#         #                saveas = os.path.join(saveat, fname)
#         #
#         #                if os.path.isfile(saveas):
#         #                    print('  Reduction already exists, skipping.')
#         #                else:
#         #                    nspec = _simp(thr, fname)
#         #                    writer.writerow([thr, nspec])
#         #        except IOError:
#         #            print('\n  Cannot write %s' % resname)
#         #            raise SystemExit
#         #        finally:
#         #            outfile.close()
#         #
#         #    # Plot residual species chart
#         #    plot_residual_species(resname)
#         #
#         #    return {'dir': saveat, 'mech': mechanism}
#         #
#         ## ****************************************************************************
#         ## function test_mechanism
#         ## ****************************************************************************
#         #
#         #def test_mechanism(solution, TPX, interval, **kwargs):
#         #    """ Test reduced mechanisms.
#         #    """
#         #    # Parse kwargs
#         #    #gonly     = xu._get_kwarg('gonly',     kwargs, None)
#         #    speciesfi = xu._get_kwarg('speciesfi', kwargs, None)
#         #    plotspec  = xu._get_kwarg('plotspec',  kwargs, None)
#         #
#         #    # Deal with arguments
#         #    msg = '\n  Missing solution directory'
#         #    directory = solution['dir'] if 'dir' in solution else sys.exit(msg)
#         #
#         #    msg = '\n  Missing solution mechanism'
#         #    mechanism = solution['mech'] if 'mech' in solution else sys.exit(msg)
#         #
#         #    # Iterate over simplifications
#         #    plt.clf()
#         #    filelist = sorted(glob.glob(os.path.join(directory,'*threshold*.txt')))
#         #    for filename in filelist:
#         #        print('\n  Testing %s' % filename)
#         #
#         #        outfilename = os.path.join(filename.split('.')[0]+'-sol.csv')
#         #        if os.path.isfile(outfilename):
#         #            print('  File already tested, skipping.')
#         #        else:
#         #            new = load_mechanism(mechanism, filename, speciesfi=speciesfi)
#         #            new.TPX = TPX
#         #            header = ['time']+new.species_names
#         #
#         #            r = ct.IdealGasConstPressureReactor(new, energy='off')
#         #            sim = ct.ReactorNet([r])
#         #
#         #            try:
#         #                print('\n  Saving %s' % outfilename)
#         #                outfile = open(outfilename, 'w')
#         #                writer = csv.writer(outfile)
#         #                writer.writerow(header)
#         #
#         #                for next_time in interval:
#         #                    sim.advance(next_time)
#         #                    writer.writerow([next_time, *new.X])
#         #            except IOError:
#         #                print('\n  Cannot write %s' % outfilename)
#         #                raise SystemExit
#         #            finally:
#         #                outfile.close()
#         #
#         #            if plotspec is not None:
#         #                # Plot simplification comparison
#         #                data = pd.read_csv(outfilename, header=0, names=header)
#         #                for spec in plotspec:
#         #                    plt.plot(data['time'], data[spec])
#         #
#         #    plt.savefig('test.png', bbox_inches='tight')
