try:
    import networkx as nx
except ModuleNotFoundError:
    import os
    os.system('pip install networkx[default]')
    
    import networkx as nx

import communityDetection
import simulation

def main():

    graph_path = 'graphs\\LFR\\50N\\graph_data_0.graphml'
    G = nx.read_graphml(graph_path)

    com = communityDetection.communityDetection(G)
    partition = com.louvain()

    steps = -1
    prob_e = 0.4
    prob_i = 0.4
    prob_r = 0.4
    prob_q = 0.3
    prob_v = 0.2
    infection_duration = 7 
    susceptible_again = False
    nat_immunity = -1


    sim = simulation.Simulation(G, partition,
                              steps  = steps,
                              prob_e = prob_e,
                              prob_i = prob_i,
                              prob_r = prob_r,
                              prob_q = prob_q,
                              prob_v = prob_v,
                              infection_duration = infection_duration,
                              susceptible_again  = susceptible_again,
                              nat_immunity       = nat_immunity)

    sim.run_scenraio(100, 2)
    sim.scenraio_results(100, "main_results")
    
main()