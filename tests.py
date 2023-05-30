try:
    import networkx as nx
except ModuleNotFoundError:
    import os
    os.system('pip install networkx[default]')
    
    import networkx as nx

import communityDetection
import simulation
import time

try:
    import numpy as np
except ModuleNotFoundError:
    import os
    os.system('pip install numpy')
    import numpy as np

try:
    from scipy.stats import entropy

except ModuleNotFoundError:
    import os 
    os.system('pip isntall scipy')
    from scipy.stats import entropy

try:
    from sklearn.metrics import normalized_mutual_info_score, mutual_info_score
    from sklearn.metrics.cluster import contingency_matrix
    from sklearn.metrics import rand_score

except ModuleNotFoundError:
    import os
    os.system('pip install scikit-learn')

    from sklearn.metrics import normalized_mutual_info_score, mutual_info_score
    from sklearn.metrics.cluster import contingency_matrix
    from sklearn.metrics import rand_score

from networkx.algorithms.community import label_propagation_communities, greedy_modularity_communities


'''
Fucntion that generates a new LFR benchmark graph. 

params:
size - number of nodes the graph should have

return:
G - networkX intance of the LFR graph with communities as a new node attribute
'''
def generateLFR_graph(size):
    
    #benchmark parameters
    tau1 = 3
    tau2 = 1.2
    mu = 0.3
    G = nx.LFR_benchmark_graph(size, tau1, tau2, mu, average_degree=5, min_community=10, seed=10)

    #Create list with communties
    communtites = []
    for i in G.nodes(data=True):
        if list(i[1]['community']) not in communtites:
            communtites.append(list(i[1]['community']))

    #Label nodes based on the intex in the 'communites' list
    for node in G.nodes(data=True):
        for i in range(len(communtites)):
            if node[0] in communtites[i]:
                node[1]['community'] = i


    return G


'''
This procedure creates LFR-bechmark graphs and saves them in the '.\graphs\LFR' location basesd on their size

params:
sizes   - a list containing the number of nodes for each graph
numbers - the number of graphs to be created with the number of nodes idicated by the elements of list [sizes]
'''
def generateLFR_data(sizes, numbers):

    path = ".\\graphs\\LFR"
    
    for size in sizes:
        #Create folder path 
        folder = "\\" + str(size)+"N"

        #Generate the number of graphs for 'size'
        for i in range(numbers):
            G = generateLFR_graph(size)

            #Create save path
            file = "\\graph_data_"+str(i)+".graphml"
            saveFile = path+folder+file
            #Write graph in new 'graphml' file. If the folder does note exists it is created now.
            try:
                nx.write_graphml(G, saveFile)
            except FileNotFoundError as e:
                import os
                os.makedirs(path+"\\{}N".format(size))
                nx.write_graphml(G, saveFile)

'''
This procedure writes the partition dictionary for the nodes in graph G, based on the information in [l_com]

params:
G      - networkX instance of graph G
part   - partition diictionary to be writen 
l_com  - list of communities represent as a collection of nodes
'''
def extract_community_label(G, part, l_com):
    for node in G.nodes:
        for index in range(len(l_com)):
            if node in l_com[index]:
                part[node] = index

#====================== Community Detection ====================#

#Compute the purity score for the labels 'y_true' and 'y_pred'
def purity_score(y_true, y_pred):
    contingency = contingency_matrix(y_true, y_pred)
    return np.sum(np.amax(contingency, axis=0)) / np.sum(contingency) 

#Compute the vi scoore for the labels 'y_true' and 'y_pred'
def vi_score(y_true, y_pred):
    mi = mutual_info_score(y_true, y_pred)
    h_true = entropy(np.bincount(y_true))
    h_pred = entropy(np.bincount(y_pred))
    vi = h_true + h_pred - 2*mi

    return vi

#Returns list of labels discovered by nx.label_propagation_communities algorithm
def label_propagation(G):
    com = label_propagation_communities(G)
    part = {node: -1 for node in G.nodes}
    
    l_com = list(com)
    
    extract_community_label(G, part, l_com)
    labels = [i for i in part.values()]
    
    return labels

#Returns list of labels discovered by nx.greedy_modularity_communitis algorithm
def greedy_modularity(G):
    
    com = greedy_modularity_communities(G)
    part = {node: -1 for node in G.nodes}
    l_com = list(list(i) for i in com)

    
    extract_community_label(G, part, l_com)
    labels = [i for i in part.values()]
    
    return labels


'''
Function that runs all the algorithms using the graphs from 'path'
!! The functoin uses the list [sizes] and [graph] to format the path string. To change 
it modify the line 127 to bring it to the correct foramt !!

The function prints the following results :
- Detection time
- Normalized Mutual Information 
- Rand Score
- Purity Score
- VI
'''
def communityDetection_AccuracyTests(path):
    
    # Graph size and index in LFR folder 
    sizes = [50,60,70,80,90,100]
    graph = [1]
    test_count = 0

    #Run tests
    for size in sizes:
        print("====== SIZE {} ========".format(size))

        for i in graph:
            graph_path = path.format(size= size, i= i)
            test_count+=1
            
            print("Test: "+ str(test_count))
            print("---- louvain ----")
            G = nx.read_graphml(graph_path)
            print(G)

            #Get the acutal list of labels for all communities in Graph G 
            y_true = [i[1]["community"]for i in G.nodes(data=True)]
            comm = communityDetection.communityDetection(G)

            #Print the actual number of communities
            print("Actual Communties: ", len(set(y_true)))
            
            #Parition the graph and record the runtime
            start =time.time()
            partitoin = comm.louvain()
            end  = time.time() - start

            #Get the list of labels predicted by the project's implementation
            y_pred =[i for i in partitoin.values()]

            #Print results for this project's implementaiton
            print("Detection time :", end,"\n")
            print("Detected Communities: ", len(set(y_pred)))
            print("Rand Score: ", rand_score(y_true, y_pred))
            print("NMI: ",normalized_mutual_info_score(y_true, y_pred))
            print("Purity Score: ",purity_score(y_true, y_pred))
            print("VI: ", vi_score(y_true, y_pred))

            
            #Print result for the label propagation algorithm 
            print("\n---- label_propagation_mx --- \n")
            y_pred_test = label_propagation(G)
            
            print("Detected Communities: ", len(set(y_pred_test)))
            print("Rand Score: ", rand_score(y_true, y_pred_test))
            print("NMI: ",normalized_mutual_info_score(y_true, y_pred_test))
            print("Purity Score: ",purity_score(y_true, y_pred))
            print("VI: ", vi_score(y_true, y_pred_test))
            
            
            #Print result for the greedy modularity algorithm 
            print("\n---- greedy_modularity_mx --- \n")
            y_pred_test = greedy_modularity(G)
            
            print("Detected Communities: ", len(set(y_pred_test)))
            print("Rand Score: ", rand_score(y_true, y_pred_test))
            print("NMI: ",normalized_mutual_info_score(y_true, y_pred_test))
            print("Purity Score: ",purity_score(y_true, y_pred))
            print("VI: ", vi_score(y_true, y_pred_test))
            print("\n\n")




#========================= Simulation ==========================#

'''
This function runs the scenario function of the Simulation class for all the scenarios
in list [scenarios]. Uses the community detection method implemented in this project. 

Prints the result metrics implemented in the Simulation class.
'''

def run_viral_simulation():

    #Load graph and find the communities
    G = nx.read_graphml("graphs\LFR\\50N\graph_data_1.graphml")
    G = nx.convert_node_labels_to_integers(G)
    com = communityDetection.communityDetection(G)
    partition = com.louvain()

    #Initialsie the simulation class
    sim = simulation.Simulation(G, partition)

    scenario_list = [0,1,2,3]
    times = 1000

    #Run simulation for each memeber of list [scenario_list] for the value of times.  
    for i in scenario_list:
        print("=========SCENARIO ",i,"========\n")

        sim.run_scenraio(times, i)
        sim.scenraio_results(times, "Scenario " + str(i))
        print("\n\n")
    




def main():
    
    path = "graphs\\LFR\\{size}N\\graph_data_{i}.graphml"

    communityDetection_AccuracyTests(path)
    run_viral_simulation()

main()
