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

        adjmat = networkx.to_numpy_matrix(self)
        plot_adjmatrix(adjmat, saveas, overwrite=overwrite, **kwargs)

    def plot_deghist(self, saveas, overwrite=True, **kwargs):
        """ Interface to :meth:`CanteraPFR.ct_aux.plot_deghist`. """

        deghist = networkx.degree_histogram(self)
        plot_deghist(deghist, saveas, overwrite=overwrite, **kwargs)

    @property
    def solution(self):
        """ Access to graph's `cantera.Solution` object. """

        return self._sol


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
                self.add_edge(A, B, val=0.0)

    def update_edges(self, src, dst, func):
        """ Update edges values.

        Parameters
        ----------
        src : list
            List of source nodes names.
        dst : list
            List of destination nodes names.
        func : function
            Function used to increment edges `val` property. This function
            takes the names of nodes `A` and `B` as arguments.
        """

        # TODO initialize edges visit here and remove reinitialize_graph.
        # TODO improve with itertools (or other).
        src = numpy.array(src)
        dst = numpy.array(dst)

        for B in src:
            for A in dst[dst != B]:
                self.edges[A, B]['val'] += func(A, B)

    def divide_edges(self, func, verbose=False):
        """ Divide edges by function.

        Parameters
        ----------
        func : function
            Function used to increment edges `val` property. This function
            takes the names of nodes `A` and `B` as arguments.
        """

        for A, B in self.edges():
            val = func(A, B)

            if val == 0.0 and verbose:
                print(f'Warning {A}->{B} has zero denominator. Setting 1.')
                self.edges[A, B]['val'] = 1.0
                continue

            ratio = self.edges[A, B]['val'] / val
            self.edges[A, B]['val'] = min(ratio, 1.0)

    def sort_edges(self, func=None, w='val', reverse=True):
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

        def mapper(pars):
            A, B = pars
            return (A, B, self.edges[A, B][w])

        if func is None:
            def func(pars):
                return pars[2]

        baselist = list(map(mapper, self.edges()))
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
        """ Set the value of property `val` of all edges to zero. """

        for A, B in self.edges():
            self.edges[A, B]['val'] = 0.0

    def complexes(self):
        """ Return complexes for building Feinberg–Horn–Jackson graph. """

        raise NotImplementedError('Not yet implemented')

    # def plot_adjmatrix(self, saveas, overwrite=True, **kwargs):
    #     print('plot_adjmatrix is not available for DirectedGraph')


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
