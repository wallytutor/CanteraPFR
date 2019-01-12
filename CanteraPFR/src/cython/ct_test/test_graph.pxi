# -*- coding: utf-8 -*-


def test_analyse_graph():
    from .ct_graph import analyse_graph
    from .ct_graph import DirectedGraph
    from .ct_graph import UndirectedGraph

    mech = 'CT-hydrocarbon-norinaga-2009-mech.cti'

    graph_class = DirectedGraph
    analyse_graph(mech, graph_class)

    graph_class = UndirectedGraph
    analyse_graph(mech, graph_class)


def test_drg_simplification():
    from .ct_graph import StateSpace
    from .ct_graph import SimplifySetup
    from .ct_graph import SimplifyDRGManager

    mech = 'CT-hydrocarbon-norinaga-2009-mech.cti'

    space = StateSpace()
    space.T = [1073, 1173]
    space.P = [1000, 5000]
    space.X = [{'N2': 0.64, 'C2H2': 0.36}]

    setup = SimplifySetup()
    setup.start_set = ['N2', 'C2H2']
    setup.plot_spec = ['C2H2', 'H2']
    setup.space = space
    setup.integ = 5.0
    setup.nsimp = 10
    setup.threshold_min = 0.01
    setup.threshold_max = 0.20
    setup.species_filter = None
    setup.dir_name = 'test_simplification'

    TPX = 1123, 3000, space.X[0]

    mngr = SimplifyDRGManager()
    all_smp = mngr.manage(mech, setup)
    mngr.test(mech, setup, TPX, all_smp, n_samples=50, transport=None)


def test_simplification():
    test_drg_simplification()
