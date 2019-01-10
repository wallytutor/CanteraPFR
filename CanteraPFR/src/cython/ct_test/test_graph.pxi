# -*- coding: utf-8 -*-


def test_analyse_graph():
    from .ct_graph import analyse_graph
    from .ct_graph import DirectedGraph
    from .ct_graph import UndirectedGraph

    mech = 'CT-hydrocarbon-norinaga-2009-mech.cti'

    graph_class = DirectedGraph
    g = analyse_graph(mech, graph_class)
    #print(g)

    graph_class = UndirectedGraph
    g = analyse_graph(mech, graph_class)
    #print(g)
