# -*- coding: utf-8 -*-
from itertools import combinations
from pathlib import Path
from typing import Any, Callable, Sequence
import cantera as ct
import networkx as nx
import numpy as np
from .mechanism import filter_mechanism
from .plot import plot_adjmatrix, plot_deghist


class BaseGraph(nx.Graph):
    """ Base class for chemical species reaction graphs.

    Parameters
    ----------
    sol : ct.Solution
        Solution object representing the kinetics mechanism.
    fullfmt : bool = False
        Print all edges in string representation if True.

    Returns
    -------
    BaseGraph
        Initialized graph instance.
    """

    __slots__ = ('_sol', '_fullfmt')

    def __init__(self, sol: ct.Solution, fullfmt: bool = False) -> None:
        # XXX no documentation here, this is intentional
        super().__init__()
        self._sol = sol
        self._fullfmt = fullfmt
        self.add_nodes_from(self._sol.species_names)

    def __str__(self) -> str:
        basefmt = (
            f' {self.__class__.__name__} object\n'
            f' -- Size   : {self.size()}\n'
            f' -- Order  : {self.order()}\n'
        )

        if not self._fullfmt:
            return basefmt

        more = map(lambda e: f' {e[0]:>20s} - {e[1]:s}\n', self.edges())
        return basefmt + ''.join(more)

    def adjacency_matrix(self) -> np.ndarray:
        """ Return adjacency matrix.

        Returns
        -------
        np.ndarray
            Adjacency matrix of interacting nodes.
        """
        return nx.to_numpy_array(self)

    def degree_histogram(self) -> list[int]:
        """ Return histogram of node degrees.

        Returns
        -------
        list[int]
            List indexed by degree with counts of its occurrences.
        """
        return nx.degree_histogram(self)

    def ordered_degree(self) -> list[tuple[str, int]]:
        """ Return nodes ordered by degree.

        Returns
        -------
        list[tuple[str, int]]
            Pairs of (node, degree) ordered descending.
        """
        baselist = [(v, self.degree(v)) for v in self]
        return sorted(baselist, key=lambda x: x[1], reverse=True)

    @property
    def density(self) -> float:
        """ Return graph density.

        Returns
        -------
        float
            Graph density.
        """
        return float(nx.density(self))

    def plot_adjmatrix(
        self,
        saveas: Path | str,
        overwrite: bool = True,
        **kwargs: object
    ) -> None:
        """ Plot graph adjacency matrix.

        Parameters
        ----------
        saveas : Path | str
            Destination image file path.
        overwrite : bool = True
            Allow overwriting existing output file.
        kwargs : object
            Additional plotting options.

        Returns
        -------
        None
            Generates plot file on disk.
        """
        adjmat = nx.to_numpy_array(self)
        plot_adjmatrix(adjmat, saveas, overwrite=overwrite, **kwargs)

    def plot_deghist(
        self,
        saveas: Path | str,
        overwrite: bool = True,
        **kwargs: object
    ) -> None:
        """ Plot node degree histogram.

        Parameters
        ----------
        saveas : Path | str
            Destination image file path.
        overwrite : bool = True
            Allow overwriting existing output file.
        kwargs : object
            Additional plotting options.

        Returns
        -------
        None
            Generates plot file on disk.
        """
        deghist = nx.degree_histogram(self)
        plot_deghist(deghist, saveas, overwrite=overwrite, **kwargs)

    @property
    def solution(self) -> ct.Solution:
        """ Return associated Cantera solution.

        Returns
        -------
        ct.Solution
            Associated kinetics solution.
        """
        return self._sol


class DirectedGraph(nx.DiGraph, BaseGraph):
    """ Generate species directed relational graph.

    Parameters
    ----------
    sol : ct.Solution
        Solution object representing the kinetics mechanism.
    fullfmt : bool = False
        Print all edges in string representation if True.

    Returns
    -------
    DirectedGraph
        Initialized directed relational graph.
    """

    __slots__ = ()

    def __init__(self, sol: ct.Solution, fullfmt: bool = False) -> None:
        # XXX no documentation here, this is intentional
        super().__init__()
        self._sol = sol
        self._fullfmt = fullfmt
        self.add_nodes_from(self._sol.species_names)

        for r in self._sol.reactions():
            reac = list(r.reactants.keys())
            prod = list(r.products.keys())
            spec = reac + prod

            self._add_edges(reac, spec)

            if r.reversible:
                self._add_edges(prod, spec)

    def _add_edges(self, src, dst):
        src_arr = np.array(src)
        dst_arr = np.array(dst)

        for b in src_arr:
            for a in dst_arr[dst_arr != b]:
                self.add_edge(a, b, val=0.0)

    def update_edges(
        self,
        src: Sequence[str],
        dst: Sequence[str],
        func: Callable[[str, str], float]
    ) -> None:
        """ Update edge weights.

        Parameters
        ----------
        src : Sequence[str]
            Source node names.
        dst : Sequence[str]
            Destination node names.
        func : Callable[[str, str], float]
            Weight increment function (A, B) -> delta.

        Returns
        -------
        None
            Updates edge weights in-place.
        """
        src_arr = np.array(src)
        dst_arr = np.array(dst)

        for b in src_arr:
            for a in dst_arr[dst_arr != b]:
                self.edges[a, b]['val'] += func(a, b)

    def divide_edges(
        self,
        func: Callable[[str, str], float],
        verbose: bool = False
    ) -> None:
        """ Divide edge weights by function values.

        Parameters
        ----------
        func : Callable[[str, str], float]
            Denominator function taking nodes A and B.
        verbose : bool = False
            Print warning messages if True.

        Returns
        -------
        None
            Divides edge weights in-place.
        """
        for a, b in self.edges():
            val = func(a, b)

            if val == 0.0 and verbose:
                print(f'Warning {a}->{b} has zero denominator. Setting 1.')
                self.edges[a, b]['val'] = 1.0
                continue

            ratio = self.edges[a, b]['val'] / val if val != 0.0 else 1.0
            self.edges[a, b]['val'] = min(ratio, 1.0)

    def sort_edges(
        self,
        func: Callable[[tuple[str, str, float]], Any] | None = None,
        w: str = 'val',
        reverse: bool = True
    ) -> list[tuple[str, str, float]]:
        """ Sort graph edges by property.

        Parameters
        ----------
        func : Callable[[tuple[str, str, float]], Any] | None = None
            Sorting key function.
        w : str = 'val'
            Edge attribute name.
        reverse : bool = True
            Sort descending if True.

        Returns
        -------
        list[tuple[str, str, float]]
            List of sorted edges (A, B, weight).
        """
        def mapper(edge):
            a, b = edge
            return (a, b, self.edges[a, b][w])

        key_func = func if func is not None else lambda pars: pars[2]
        baselist = list(map(mapper, self.edges()))
        return sorted(baselist, key=key_func, reverse=reverse)

    def ordered_in_degree(self) -> list[tuple[str, int]]:
        """ Return nodes ordered by in-degree.

        Returns
        -------
        list[tuple[str, int]]
            Pairs of (node, in-degree) ordered descending.
        """
        baselist = [(v, self.in_degree(v)) for v in self]
        return sorted(baselist, key=lambda x: x[1], reverse=True)

    def ordered_out_degree(self) -> list[tuple[str, int]]:
        """ Return nodes ordered by out-degree.

        Returns
        -------
        list[tuple[str, int]]
            Pairs of (node, out-degree) ordered descending.
        """
        baselist = [(v, self.out_degree(v)) for v in self]
        return sorted(baselist, key=lambda x: x[1], reverse=True)

    def reinitialize_graph(self) -> None:
        """ Reset all edge values to zero.

        Returns
        -------
        None
            Resets edge attributes in-place.
        """
        for a, b in self.edges():
            self.edges[a, b]['val'] = 0.0


class UndirectedGraph(BaseGraph):
    """ Generate species undirected relational graph.

    Parameters
    ----------
    sol : ct.Solution
        Solution object representing the kinetics mechanism.
    fullfmt : bool = True
        Print all edges in string representation if True.

    Returns
    -------
    UndirectedGraph
        Initialized undirected graph.
    """

    __slots__ = ()

    def __init__(self, sol: ct.Solution, fullfmt: bool = True) -> None:
        # XXX no documentation here, this is intentional
        super().__init__(sol, fullfmt=fullfmt)

        for r in self._sol.reactions():
            reac = r.reactants.keys()
            prod = r.products.keys()
            spec = list(reac) + list(prod)

            for p, q in combinations(spec, 2):
                if (p, q) not in self and (p != q):
                    self.add_edge(p, q)


def analyse_graph(
        mech: Path | str,
        graph_class: type[BaseGraph],
        **kwargs: object
    ) -> BaseGraph:
    """ Perform species graph analysis.

    Parameters
    ----------
    mech : Path | str
        Mechanism file path.
    graph_class : type[BaseGraph]
        Graph class to construct.
    kwargs : object
        Optional species filtering and plotting arguments.

    Returns
    -------
    BaseGraph
        Constructed and analyzed graph instance.
    """
    saveas = kwargs.get('saveas', f'results_{graph_class.__name__}')
    species = kwargs.get('species', None)
    species_seq = list(species) if species is not None else None

    sol = filter_mechanism(mech, species_seq, **kwargs)
    mgraph = graph_class(sol)
    save_path = Path(str(saveas))
    save_stem = str(save_path)

    if isinstance(mgraph, DirectedGraph):
        with Path(f'{save_stem}-inward.txt').open('w', encoding='utf-8') as w:
            for v, d in mgraph.ordered_in_degree():
                w.write(f'{v:16s} {d:4d}\n')

        with Path(f'{save_stem}-outward.txt').open('w', encoding='utf-8') as w:
            for v, d in mgraph.ordered_out_degree():
                w.write(f'{v:16s} {d:4d}\n')

    with Path(f'{save_stem}-global.txt').open('w', encoding='utf-8') as w:
        for v, d in mgraph.ordered_degree():
            w.write(f'{v:16s} {d:4d}\n')

    with Path(f'{save_stem}-report.txt').open('w', encoding='utf-8') as w:
        basestr = (
            f'Species graph contains {mgraph.size()} edges\n'
            f'Species graph contains {mgraph.order()} nodes\n'
            f'Species graph density  {mgraph.density}\n'
        )
        w.write(basestr)

    mgraph.plot_deghist(Path(f'{save_stem}-histogram.png'))
    mgraph.plot_adjmatrix(Path(f'{save_stem}-adjacency.png'))
    return mgraph
