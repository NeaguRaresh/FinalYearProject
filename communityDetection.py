import networkx as nx


class communityDetection:

    '''
    Community Detection Class:

    It implements the Louvain Method of maximsing the modularity score.

    Params:
    G - networkX instance of a graph

    Returns:
    partittion - dictonary with the node as Key and coresponding community as Value
    '''
    def __init__(self, G):
        
        self.G = nx.convert_node_labels_to_integers(G)
        
        #Initialise the partition dict with each node in its own community
        self.partition = {i: i for i in self.G.nodes()}

        #Pre calculate the degree of each node to speed up the 'get_modularity' function
        self.degrees = {i[0]: i[1] for i in self.G.degree()}

        self.savePath = "results\\community_graphs\\"
        

    '''
    This function calculates the modularity based on the function

    Q = l/m - ((k_i / 2m) * (d / 2m) - (k_tot /2m)^2 - (k_i / 2m)^2

    where

    m - number of edges in the graph G
    l - number of edges in one community
    
    d     - degree of nodes in a community 
    k_i   - degree of nodes in one community including the degree of adjacent nodes
    k_tot - degree of node in the whole graph
    ---------------
    Params:
    parittion - dictionary representing the current partition, where they Key is the node and the Value is the community

    Returns:
    Q - float representing the modularity of current partition
    '''
    def get_modularity(self, partition):
    
        m = self.G.number_of_edges()
        q = 0

        #Repeat for all current communitites
        for i in set(partition.values()):

            #Create the subgraph coresponding to the community 'i'
            nodes = [n for n in partition if partition[n] == i]
            subgraph = self.G.subgraph(nodes)
            l = subgraph.number_of_edges()

            #Calculate the value of q for this community
            d = sum(dict(self.G.degree(subgraph)).values())
            
            k_i = sum([self.degrees[node] for node in nodes])
            
            k_tot = sum(self.degrees.values())
            q += (l / m) - ((k_i / (2 * m)) * (d / (2 * m)))
            q -= ((k_tot / (2 * m)) ** 2) - ((k_i / (2 * m)) ** 2)

        #Return the total value of modularity q for the partition 'partition'
        return q


    '''
    This function merges the communitites based on the maximum value of modularity

    G         - networkX instance of graph G
    partition - current partitioning of the graph G

    returns:
    partition - dictionary representing the improved partitioning of the graph
    '''
    def merge_communities(self, G, partition):
        
        nodes = list(partition.keys())
        
        #Compute modularity for the current partitioning
        q = self.get_modularity( partition)
        
        #Repeat until no improvements can be made
        while True:
            improvement = False

            #Iterate over all the node and it's neighbors to find the best community for the node 'i'
            for i in G.nodes:
                for j in G.neighbors(i):

                    node_i = nodes[i]
                    node_j = nodes[j]

                    #If node_i and node_j are not part of the same partition check if their communties can be merged              
                    if partition[node_i] != partition[node_j]:
                        
                        #Assign the node_i to the community of node j
                        new_partition = partition.copy()
                        new_partition[node_i] = new_partition[node_j]

                        #Check if the new modularity is better than the one before                        
                        new_q = self.get_modularity( new_partition)
                        
                        if new_q > q:
                            #Save new partition and modularity value
                            partition = new_partition
                            q = new_q
                            improvement = True
            
            if not improvement:
                break
        return partition


    '''
    This function calls the 'merge_communities' until the maximum value of modularity, and the 
    best partitioning is found. 

    returns:
    partition - dictionary representing the best partititioning of the graph, where they Key is the 
                node and the Value is the community
    '''
    def louvain(self):
       

        
        partition = self.partition.copy()

        while True:
            #Find improved partition
            new_partition = self.merge_communities(self.G, partition)

            #If the last partition is equal to the new one, the best solutoin was found and loop end
            if sorted(new_partition.items()) == sorted(partition.items()):
                break

            #Save new communitites dict
            partition = new_partition
            

        # return the final partition 
        self.partition = partition
        return partition


    '''
    Creates new graphml file with the same number of nodes and edges as self.G, but 
    adds the 'community' node attribute representing the community label each node is part of.

    Saves the graph in the self.savePath location
    '''
    def createCommunityLabelGraph(self):

        
        #Create copy of graph G 
        g = nx.Graph()
        g.add_nodes_from(self.partition.keys())
        g.add_edges_from(self.G.edges)
        
        #Create the new node attributes
        value_list = [i for i in self.partition.values()]
        node_attrs = {node:{'community': value_list[node]} for node in g.nodes}

        #Add new attributes to the nodes
        nx.set_node_attributes(g, node_attrs)

        #Save file
        nx.write_graphml(g, self.savePath+"Community_Labels_Graph.graphml")
        
    '''
    Creates new graphml file with a graph representing the communities and the relationships between them.
    For each node it creates an attribute called 'population' that represents the number of members each 
    community has. Edge weight of edge (u,v) represents the number of members' relationiships from community 'u'
    to community 'v'.


    Saves the graph in the self.savePath location
    '''
    def createCommunityGraph(self):
        
        
        #Create community graph
        g = nx.Graph()
        g.add_nodes_from(self.partition.values())

        #Adds new node attributes to the graph
        value_list =[i for i in self.partition.values()]
        node_attrs = {node: {'population': value_list.count(value_list[node])} for node in g.nodes}
        nx.set_node_attributes(g, node_attrs)
        
            
        #Add weighted edges to the graph.
        for i in self.G.edges:
            
            #Avoid self loops
            if self.partition[i[0]] == self.partition[i[1]]:
                continue

            #Add edge between communities
            g.add_edge(self.partition[i[0]],self.partition[i[1]])
            
            #Update the 'weight' attribute of the edge
            if 'weight' in g[self.partition[i[0]]][self.partition[i[1]]]:
                g[self.partition[i[0]]][self.partition[i[1]]]['weight'] += 1
            else:
                g[self.partition[i[0]]][self.partition[i[1]]]['weight'] = 1
        
        #Save file
        nx.write_graphml(g, self.savePath+"\\Community_Graph_.graphml")

