# # Recycle XCantera Analyse.py here
#
# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Wed Sep  6 06:23:35 2017
#
# @author: dalmazsi1
# """
#
# import os
# import itertools
# import numpy as np
# import networkx as nx
# import matplotlib.pyplot as plt
# from matplotlib.ticker import MaxNLocator
#
# plt.style.use('seaborn')
# # TODO add frame to plots.
# # TODO document kwargs.
# # TODO improve descriptions after updating simplification routines.
# # TODO correct 'Returns' in functions.
#
#
# class BaseGraph(nx.Graph):
#     def __init__(self):
#         """ A base classe providing graphics utilities for derived classes.
#
#             It can be also used standalone to plot existing networkx graphs.
#         """
#         super(BaseGraph, self).__init__()
#
#     def __str__(self):
#         """ Formatted string representation. """
#         size, order = self.size(), self.order()
#         basefmt = ' {} object\n -- Size   : {}\n -- Order  : {}\n'
#         basefmt = basefmt.format(self.__class__.__name__, size, order)
#         if self._fullfmt:
#             for edge in self.edges():
#                 basefmt += ' {:>20s} - {:s}\n'.format(*edge)
#         return basefmt
#
#     def adjacency_matrix(self):
#         """ Returns a numpy array from the adjacency matrix.
#
#         Returns
#         -------
#         A :obj:`numpy.array`.
#         """
#         return nx.to_numpy_matrix(self)
#
#     def degree_histogram(self):
#         """ Return a list of the frequency of each degree value.
#
#         Returns
#         -------
#         A :obj:`list` of :obj:`int`.
#         """
#         return nx.degree_histogram(self)
#
#     def ordered_degree(self):
#         """ Returns list of tuples (node, degree) ordered by degree.
#
#         Returns
#         -------
#         A :obj:`list` of :obj:`tuples` of (:obj:`str`, :obj:`int`).
#         """
#         baselist = [(v, self.degree(v)) for v in self]
#         return sorted(baselist, key=lambda x: x[1], reverse=True)
#
#     @property
#     def density(self):
#         """ :obj:`float` : Graph density. """
#         return nx.density(self)
#
#     @classmethod
#     def plot_adjmatrix(cls, saveas, graph=None, **kwargs):
#         """ Plot the graph's adjacency matrix.
#
#         Parameters
#         ----------
#         saveas : :obj:`str`
#             Path to save file. Returns `None` if file exists.
#         graph : :obj:`BaseGraph`
#             Source graph.
#         """
#         if os.path.exists(saveas):
#             print(' Graph cannot overwrite: {}'.format(saveas))
#             return
#
#         if graph is None:
#             if hasattr(cls, '_graph'):
#                 graph = cls._graph
#             else:
#                 print(' Class not set for use with graphical output')
#                 return
#
#         figsize = kwargs.get('figsize', (6, 6))
#         cmap = kwargs.get('cmap', 'Greys')
#         dpi = kwargs.get('dpi', 300)
#
#         plt.clf()
#         plt.close('all')
#         fig, ax = plt.subplots(figsize=figsize)
#         ax.imshow(nx.to_numpy_matrix(graph), cmap=cmap)
#         ax.axes.get_xaxis().set_visible(False)
#         ax.axes.get_yaxis().set_visible(False)
#         fig.subplots_adjust(top=1.0, bottom=0.0, right=1.0, left=0.0)
#         plt.savefig(saveas, dpi=dpi)
#         plt.close('all')
#
#     @classmethod
#     def plot_histogram(cls, saveas, graph=None, **kwargs):
#         """ Plot the histogram of graph degrees.
#
#         Parameters
#         ----------
#         saveas : :obj:`str`
#             Path to save file. Returns `None` if file exists.
#         graph : :obj:`BaseGraph`
#             Source graph.
#         """
#         if os.path.exists(saveas):
#             print(' Graph cannot overwrite: {}'.format(saveas))
#             return
#
#         if graph is None:
#             if hasattr(cls, '_graph'):
#                 graph = cls._graph
#             else:
#                 print(' Class not set for use with graphical output')
#                 return
#
#         figsize = kwargs.get('figsize', (6, 6))
#         color = kwargs.get('color', 'k')
#         xlim = kwargs.get('xlim', None)
#         dpi = kwargs.get('dpi', 300)
#
#         plt.clf()
#         plt.close('all')
#         deghis = nx.degree_histogram(graph)
#         fig, ax = plt.subplots(figsize=figsize)
#         ax.bar(np.arange(0, len(deghis), 1), deghis, color=color)
#         ax.set_xlabel('Degree')
#         ax.set_ylabel('Counts')
#         if type(xlim) is not dict:
#             ax.set_xlim(xlim)
#         else:
#             ax.set_xlim(**xlim)
#         ax.xaxis.set_major_locator(MaxNLocator(integer=True))
#         ax.yaxis.set_major_locator(MaxNLocator(integer=True))
#         plt.savefig(saveas, dpi=dpi, )
#         plt.close('all')
#
#
# class DGraph(nx.DiGraph, BaseGraph):
#     def __init__(self, sol, fullfmt=True):
#         """ Generate species directed relational graph.
#
#         An edge is added for every pair of self influencing species in a
#         reaction.
#
#         Note
#         ----
#             No self loops :math:`(A\\rightarrow{}A)` are generated in the
#             graph.
#
#         Parameters
#         ----------
#         sol : :obj:`cantera.Solution`
#             Solution object representing the kinetics mechanism.
#         fullfmt : :obj:`bool`
#             If `True`, print all edges in mechanism formatted string.
#             Default is `True`.
#         """
#         super(DGraph, self).__init__()
#         self._sol, self._fullfmt = sol, fullfmt
#         self._build_graph()
#
#     def _build_graph(self):
#         """ Connect species from reaction mechanism. """
#         # Each species represents a node.
#         self.add_nodes_from(self._sol.species_names)
#
#         # Pairwise combine species from reactions.
#         for r in self._sol.reactions():
#             reac = list(r.reactants.keys())
#             prod = list(r.products.keys())
#             spec = reac + prod
#
#             self._add_edges(reac, spec)
#             if r.reversible:
#                 self._add_edges(prod, spec)
#
#         # Do not allow graph to be modified.
#         self = nx.freeze(self)
#
#         # To be able to call graphical output.
#         self.__class__._graph = self
#
#     def _add_edges(self, src, dst):
#         """ Add edges from src to dst without self loops on nodes. """
#         # TODO improve with itertools (or other).
#         src, dst = np.array(src), np.array(dst)
#         for B in src:
#             for A in dst[dst != B]:
#                 self.add_edge(A, B, {'val': 0.0})
#
#     def update_edges(self, src, dst, func):
#         """ Update edges values.
#
#         Parameters
#         ----------
#         src : :obj:`list` of :obj:`str`
#             List of source nodes names.
#         dst : :obj:`list` of :obj:`str`
#             List of destination nodes names.
#         func : :obj:`function`
#             Function used to increment nodes `val` property.
#         """
#         # TODO initialize edges visit here and remove reinitialize_graph.
#         src, dst = np.array(src), np.array(dst)
#         for B in src:
#             for A in dst[dst != B]:
#                 self[A][B]['val'] += func(A, B)
#
#     def divide_nodes(self, func):
#         """ Divide nodes by function.
#
#         Parameters
#         ----------
#         func : :obj:`function`
#             Function used to divide nodes `val` property.
#         """
#         for A, B in nx.edges_iter(self):
#             val = func(A, B)
#             if val == 0.0:
#                 print('Warning {:s} has denominator={:.2f}'.format(A, val))
#                 self[A][B]['val'] = 1.0
#             else:
#                 ratio = self[A][B]['val'] / val
#                 self[A][B]['val'] = min(ratio, 1.0)
#
#     def sort_edges(self, func, reverse=True):
#         """ Sort graph edges by key using a function.
#
#         Parameters
#         ----------
#         func : :obj:`function`
#             Function to provide ordering.
#         reverse : :obj:`bool`, optional
#             If `True`, returns in reverse order.
#
#         Returns
#         -------
#         A :obj:`list` of :obj:`str`.
#         """
#         G, k = self, 'val'
#         baselist = [(A, B, {k: G[A][B][k]}) for A, B in nx.edges_iter(G)]
#         return sorted(baselist, key=func, reverse=reverse)
#
#     def ordered_in_degree(self):
#         """ Returns list of tuples (node, degree) ordered by degree.
#
#         Returns
#         -------
#         A :obj:`list` of :obj:`tuples` of (:obj:`str`, :obj:`int`).
#         """
#         baselist = [(v, self.in_degree(v)) for v in self]
#         return sorted(baselist, key=lambda x: x[1], reverse=True)
#
#     def ordered_out_degree(self):
#         """ Returns list of tuples (node, degree) ordered by degree.
#
#         Returns
#         -------
#         A :obj:`list` of :obj:`tuples` of (:obj:`str`, :obj:`int`).
#         """
#         baselist = [(v, self.out_degree(v)) for v in self]
#         return sorted(baselist, key=lambda x: x[1], reverse=True)
#
#     def reinitialize_graph(self):
#         """ Set the value of property `val` of all edges to zero. """
#         for A, B in nx.edges_iter(self):
#             self[A][B]['val'] = 0.0
#
#     def complexes(self):
#         """ Return complexes for building Feinberg–Horn–Jackson graph. """
#         raise NotImplemented('Not yet implemented')
#
#
# class UGraph(BaseGraph):
#     def __init__(self, sol, fullfmt=True):
#         """ Generate species undirected relational graph.
#
#         An edge is added for every pair of species in a reaction.
#
#         Note
#         ----
#             No self loops :math:`(A\\rightarrow{}A)` are generated in the
#             graph.
#
#         Parameters
#         ----------
#         sol : :obj:`cantera.Solution`
#             Solution object representing the kinetics mechanism.
#         fullfmt : :obj:`bool`, optional
#             If `True`, print all edges in mechanism formatted string.
#             Default is `True`.
#         """
#         super(UGraph, self).__init__()
#         self._sol, self._fullfmt = sol, fullfmt
#         self._build_graph()
#
#     def _build_graph(self):
#         """ Connect species from reaction mechanism. """
#         # Each species represents a node.
#         self.add_nodes_from(self._sol.species_names)
#
#         # Pairwise combine species from reactions.
#         for r in self._sol.reactions():
#             spec = list(r.reactants.keys()) + list(r.products.keys())
#             for p, q in list(itertools.combinations(spec, 2)):
#                 if (p, q) not in self and (p != q):
#                     self.add_edge(p, q)
#
#         # Do not allow graph to be modified.
#         self = nx.freeze(self)
#
#         # To be able to call graphical output.
#         self.__class__._graph = self
#
#         #!/usr/bin/env python3
#         # -*- coding: utf-8 -*-
#         """
#         Created on Fri Sep  8 06:43:56 2017
#
#         @author: dalmazsi1
#         """
#
#         import os
#         import itertools
#         import cantera as ct
#         from . import UGraph
#         from . import DGraph
#         from . import ChemUtils
#
#
#         def mechanism_from_species(mechanism, species, transport=None):
#             """ Filter mechanism from species list in file or list.
#
#             A solution object is filtered from a list of species, which can be stored
#             in a file or in a Python :obj:`list` of :obj:`str`.  If the object type of
#             `species` is :obj:`str`, the method will try to read the species list from
#             a space separated file.  If this fails, it will try to convert the string
#             into a list by splitting at white spaces. If a :obj:`list` is provided, it
#             is directly employed by the filter routine.
#
#             Parameters
#             ----------
#             mechanism : :obj:`str`
#                 Mechanism file path in Cantera's format.
#             species : :obj:`str`, :obj:`list` or other iterable
#                 File path or list of species to keep.
#
#             Returns
#             -------
#             :obj:`cantera.Solution`
#                 Solution object.
#             """
#             try:
#                 if type(species) is str:
#                     gonly = open(species, 'r').read().split(sep=' ')
#                 else:
#                     gonly = species
#             except (IOError, Exception) as err:
#                 print(err)
#                 if type(species) is str:
#                     gonly = species.split(sep=' ')
#             return ChemUtils.filterSpecies(mechanism, gonly, transport=transport)
#
#
#         def analyse_graph(mechanism, grapht, **kwargs):
#             """ Automatic analysis of species graph.
#
#             This method provides a preliminary mechanism species graph analysis.
#             It is able to filter species and reactions from a species file and
#             perform analysis of simpler versions of a mechanism.  An input file
#             and graph type are required inputs.  Optional arguments are the base
#             name to save the resulting files and the list of species to filter
#             the mechanism.  Graph species connection degree output is supplied
#             as well as a basic report with the number of nodes, edges and density.
#             Plots adjacency matrix graph histogram.
#
#             Parameters
#             ----------
#             mechanism: :obj:`str`
#                 Mechanism file path in Cantera's format.
#             grapht: :obj:`str`
#                 Graph type, 'UGraph' or 'DGraph'.
#
#             Returns
#             -------
#             :obj:`grapht`
#                 Return generated graph depending on the type name `grapht`.
#             """
#             print(' Starting {} analysis approach'.format(grapht))
#             saveas = kwargs.get('saveas', 'results-{}'.format(grapht))
#             species_list = kwargs.get('species_list', None)
#
#             if species_list is None:
#                 sol = ct.Solution(mechanism)
#             else:
#                 sol = mechanism_from_species(mechanism, species_list)
#
#             if grapht not in ['UGraph', 'DGraph']:
#                 raise ValueError('Unexpected graph type {}'.format(grapht))
#             elif grapht == 'UGraph':
#                 mgraph = UGraph(sol)
#             elif grapht == 'DGraph':
#                 mgraph = DGraph(sol)
#                 ordidg = mgraph.ordered_in_degree()
#                 fname = '{}-inward.txt'.format(saveas)
#                 with open(fname, 'w') as writer:
#                     for v, d in ordidg:
#                         writer.write('{:16s} {:4d}\n'.format(v, d))
#                 if writer:
#                     writer.close()
#                 ordodg = mgraph.ordered_out_degree()
#                 fname = '{}-outward.txt'.format(saveas)
#                 with open(fname, 'w') as writer:
#                     for v, d in ordodg:
#                         writer.write('{:16s} {:4d}\n'.format(v, d))
#                 if writer:
#                     writer.close()
#
#             orddeg = mgraph.ordered_degree()
#             fname = '{}-global.txt'.format(saveas)
#             with open(fname, 'w') as writer:
#                 for v, d in orddeg:
#                     writer.write('{:16s} {:4d}\n'.format(v, d))
#             if writer:
#                 writer.close()
#
#             fname = '{}-report.txt'.format(saveas)
#             with open(fname, 'w') as writer:
#                 basestr = ('Species graph contains {} edges\n'
#                            'Species graph contains {} nodes\n'
#                            'Species graph density  {}\n')
#                 basestr = basestr.format(mgraph.size(), mgraph.order(), mgraph.density)
#                 writer.write(basestr)
#             if writer:
#                 writer.close()
#
#             mgraph.plot_histogram('{}-histogram.png'.format(saveas))
#             mgraph.plot_adjmatrix('{}-adjacency.png'.format(saveas))
#             return mgraph
#
#
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
