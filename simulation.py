from random import choice, random
try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    import os
    os.system('pip install matplotlib')

try:
    import networkx as nx
except ModuleNotFoundError:
    import os
    os.system('pip install networkx[default]')



class Simulation:


    '''
    This class represents a SEIR (Susceptible - Exposed - Infected - Recovered) viral spread model. Nodes of graph 'G' are categorized into one of the four
    categories and are allowed to change state only if the requrirements are met. This class has multiple data strucutre created to keep track of the nodes state at all times
    during the simulation and scenario runs. When itializing, the user can change all of the following simulation parameters:

    G         - NetworkX instacnce of grpah G representing the network the virus will act upon
    steps     - number of days each simulation runs, if value is -1 simulation continues unill no more nodes can be infected
    partitoin - dictionary representing the communities nodes belong to. Keys reprsent the node and the values their respective communities
    
    prob_e - probability of exposure given contact
    prob_i - probability of getting the disease given exposure
    prob_r - probability to recover from the disease
    prob_v - probability that a node will accept the vaccine
    prob_q - probability that a node will respect the quarantine

    infection_duration - number of days a node has to be infected before it can recover
    susceptible_again  - control variable that allows nodes to become susceptible again after they recovered
    nat_immunity       - number of days a node stays in the recovered state if 'susceptible_again' is active 
    '''
    def __init__(self, G, partition, steps = -1, prob_e = 0.4, prob_i = 0.4, prob_r = 0.3, prob_q = 0.3, prob_v = 0.2, infection_duration = 7, susceptible_again = False, nat_immunity = -1):

        #Simulation graph and its community structure
        self.G = nx.convert_node_labels_to_integers(G)



        self.partition = partition

        #Number of repretitions, -1 means the simulation runs until no more nodes can be infected
        self.steps = steps

        #Values for the simulation parameters
        self.prob_e = prob_e
        self.prob_i = prob_i
        self.prob_r = prob_r
        self.prob_q = prob_q
        self.prob_v = prob_v

        #Control variable - Dictates if the nodes can become susceptible again
        self.susceptible_again=susceptible_again

        #Duration each node stays infected before it can recover
        self.infection_duration = infection_duration

        #Duration for which nodes stay in recover state if 'susceptible_again' is True
        self.nat_immunity = nat_immunity
        
        self.refresh_states()


    '''
    This method resets the metric data structures that keep track of nodes states 
    and other information about the simulation to defalt values

    t                - tracks the number of days the simulation runs
    timekeeping_list - dictionary that records the time each node enters a new state
    states           - dictionary that tracks the nodes in each state at time 't'
    infector_dict    - dictionary that keeps track who infected each node
    '''
    def refresh_states(self):
        self.timekeeping_list = {node: {"EXPOSED_TIME": 0, 
                                     "INFECTED_TIME": 0,
                                     "RECOVERED_TIME": 0,
                                     "QUARANTINE_TIME": 0,
                                     "VACCINATED_TIME": 0} for node in self.G.nodes} 
      
        self.t = 0

        self.states = {self.t: {"SUSCEPTIBLE": set(),
                                "EXPOSED": set(),
                                "INFECTED": set(),
                                "RECOVERED": set(),
                                "QUARANTINED" : set(),
                                "VACCINATED" : set()},
                        }
            
        self.infector_dict = {node: -1 for node in self.G} 


    '''
    This method runs a "scenario". A scenario is a repretition of the same simulation for a determined number of times.

    params:
    times               - number of repetitions for this scenario
    protective_measures - the indicator for the protective measure to be used in the scenario

    modifies the [stats_scenario] list to keep track of each simulation metric data structures 
    '''
    def run_scenraio(self, times, protective_measures):

        time = 0
        
        self.stats_scenario = []
        
        while time < times:

            self.run_simulation(protective_measures)
            self.stats_scenario.append((self.t, self.states, self.infector_dict))

            self.refresh_states()
            time+=1


    '''
    This method begins a new simulation. 

    protective_measures - the indicator for the protective measures to be used in the simulation
            0 - No protective measures
            1 - Vaccination campaing
            2 - Quaraninte orders
            3 - Vaccination and Quaraninte 
    '''
    def run_simulation(self, protecitve_measures):
       
        #Choose the first infected node at random
        pacient_zero = choice(list(self.G.nodes()))
        
        #Make all nodes susceptible
        self.states[self.t]["SUSCEPTIBLE"].update(set([node for node in self.G.nodes() if node != pacient_zero]))

        #Add pacient_zero to metrics 
        self.states[self.t]["INFECTED"].add(pacient_zero)
        self.timekeeping_list[pacient_zero]["INFECTED"] = self.t
        
        #pacient zero is recored as infecting itself
        self.infector_dict[pacient_zero] = pacient_zero


        while self.states[self.t]["INFECTED"] or self.states[self.t]["EXPOSED"]:
        
            #Update time and metric data structures
            self.t +=1
            self.add_simulation_step()

            #infected nodes expose suscetible nodes            
            self.exposure_step()
        
            #Exposed nodes develop the virus
            self.infection_step()

            #Infected nodes recover
            self.recovery_step()


            # Include Preventive measure steps
            if protecitve_measures == 1 or protecitve_measures == 3:
                self.vaccination_step()

            if protecitve_measures == 2 or protecitve_measures == 3:
                self.quarantine_step()

            #Check if the number of steps for this simulation have been reached
            if self.steps < self.t and self.steps > -1:
                break

    '''
    In this method the SUSCEPTIBLE nodes are EXPOSED to the virus. A node can become inefcted if they are not VACCINATED.

    The longer a node is infected the more infectious it is. Each day the node is 'deltaT' / 10 more infectious, where deltaT 
    is the number of days the node was infected.   
    '''
    def exposure_step(self):
        
        
        for infected in self.states[self.t-1]["INFECTED"]:
            
            #If node is qurantined it cannot infect other nodes
            if infected in self.states[self.t-1]["QUARANTINED"]:
                continue
            
            #Find the time this node was infected for
            deltaT = self.t - self.timekeeping_list[infected]["INFECTED_TIME"]    
            
            #Node is meeting its neighbours
            for neighbour in self.G.neighbors(infected):

                #If a neighbour is vaccinated it cannot contact the virus
                if neighbour in self.states[self.t]["VACCINATED"]:
                    continue
                
                #If neigbour can contact the virus
                if neighbour in self.states[self.t]["SUSCEPTIBLE"]:
                    
                    # Determine if the neighbor contacts the virus based on the time it was infected
                    if random() < self.prob_e + (deltaT / 10):
                        self.states[self.t]["EXPOSED"].add(neighbour)
                        self.states[self.t]["SUSCEPTIBLE"].remove(neighbour)
                        self.timekeeping_list[neighbour]["EXPOSED_TIME"] = self.t-1
                        self.infector_dict[neighbour] = infected
                        
    '''
    In this method the EXPOSED nodes become INFECTED. The longer a node is EXPOSED the less likely it is to be come INFECTED.

    If a  node is in EXPOSED state for more than 3 days it becomes SUSCEPTIVLE. If the number of days is less than 3, the node
    has a chance of 'prob_i' * (1 / deltaT), reducing the chances to become infected. 
    '''         
    def infection_step(self):
        
        
        for exposed in self.states[self.t]["EXPOSED"]:
            deltaT = self.t - self.timekeeping_list[exposed]["EXPOSED_TIME"]

            #If the number of days is less than 3
            if deltaT < 3: 

                #If deltaT is 0 the probability of infection is not affected
                if deltaT == 0:

                    if random() < self.prob_i:
                        self.states[self.t]["INFECTED"].add(exposed)
                        self.timekeeping_list[exposed]["INFECTED_TIME"] = self.t
                #If detltaT is greather than 0, the probability of infection decreases
                elif random() < self.prob_i * (1 / deltaT):
                    self.states[self.t]["INFECTED"].add(exposed)
                    self.timekeeping_list[exposed]["INFECTED_TIME"] = self.t
            #Node beocme Susceptible after 3 days 
            else:   
                self.states[self.t]["SUSCEPTIBLE"].add(exposed)

        #Update State lists
        self.states[self.t]["EXPOSED"] -= self.states[self.t]["EXPOSED"].intersection(self.states[self.t]["INFECTED"])
        self.states[self.t]["EXPOSED"] -= self.states[self.t]["EXPOSED"].intersection(self.states[self.t]["SUSCEPTIBLE"])
        

    '''
    In this method INFECTED nodes become RECOVERED. a node is infected for a minimum of 'infection_duration' days. After that
    number of days have passed the nodes have a probability of 'prob_r' to recover. 

    If the 'susceptible_again' variable is True, RECOVERED nodes return tu SUSCEPTIBLE state if they have been RECOVERED more than
    'nat_immunity' value. 
    '''
    def recovery_step(self):
   
        
        for infected in self.states[self.t]["INFECTED"]:  
            #If note is infected for less than 'infection_duration' 
            if self.t - self.timekeeping_list[infected]["INFECTED_TIME"] < self.infection_duration:
                continue
            if random() < self.prob_r:
                self.states[self.t]["RECOVERED"].add(infected)
                self.timekeeping_list[infected]["RECOVERED_TIME"] = self.t
        
        self.states[self.t]["INFECTED"] -= self.states[self.t]["INFECTED"].intersection(self.states[self.t]["RECOVERED"])
        self.states[self.t]["QUARANTINED"] -= self.states[self.t]["QUARANTINED"].intersection(self.states[self.t]["RECOVERED"])

        if self.susceptible_again:
            for recovered in self.states[self.t]["RECOVERED"]:
                if self.t - self.timekeeping_list[infected]["RECOVERED_TIME"] < self.nat_immunity:
                    self.states[self.t]["SUSCEPTIBLE"].add(recovered)

        self.states[self.t]["RECOVERED"] -= self.states[self.t]["RECOVERED"].intersection(self.states[self.t]["SUSCEPTIBLE"])


    '''
    This method adds a new entry in the {states} dictionary
    '''
    def add_simulation_step(self):
        self.states.setdefault(self.t, {})
        self.states[self.t].setdefault("EXPOSED", self.states[self.t-1]["EXPOSED"].copy())
        self.states[self.t].setdefault("INFECTED", self.states[self.t-1]["INFECTED"].copy())
        self.states[self.t].setdefault("RECOVERED", self.states[self.t-1]["RECOVERED"].copy())
        self.states[self.t].setdefault("SUSCEPTIBLE", self.states[self.t-1]["SUSCEPTIBLE"].copy())
        self.states[self.t].setdefault("VACCINATED", self.states[self.t-1]["VACCINATED"].copy())
        self.states[self.t].setdefault("QUARANTINED", self.states[self.t-1]["QUARANTINED"].copy())


    '''
    In this method INFECTED nodes change state to QUARANTINED. The proability of an INFECTED node to become QUARANTINED is 'prob_q'.

    A QUARANTINED node cannot infect SUSCEPTIBLE nodes anymore.
    '''
    def quarantine_step(self):
        
        for node in self.states[self.t]["INFECTED"]:
            if random() < self.prob_q:
                self.states[self.t]["QUARANTINED"].add(node)


    '''
    In this method SUSCEPTIBLE nodes change state to VACCINATED. The proability of an SUSCEPTIBLE node to become VACCIANTED is 'prob_v'

    A VACCINATED node cannot contact the virus anymore.
    '''
    def vaccination_step(self):
        for node in self.states[self.t]["SUSCEPTIBLE"]:
            if random() < self.prob_v:
                self.states[self.t]["VACCINATED"].add(node)

        self.states[self.t]["SUSCEPTIBLE"] -= self.states[self.t]["VACCINATED"]

############## Graphic Representation & Numerical data ##################

    '''
    This method extracts the number of nodes in each state for 't' steps of the simulation

    params:
    t      - number of simulation steps
    states - dictionary of shape states[ time ][ node_state ]
    show   - True for data to be ploted on graph, False otherwise

    returns:
    touple - touple of lists with the number of node is each state 
            (Susceptible, Exposed, Infected, Recovered, Vaccinated, Quarantined)
    '''
    def simulationGraphic(self, states, t, show=False):
       
        #Creates list containing the number of nodes in a state at anytime in the simulation
        sus = [len(states[t]['SUSCEPTIBLE']) for t in range(t)]
        exp = [len(states[t]['EXPOSED']) for t in range(t)]
        inf = [len(states[t]['INFECTED']) for t in range(t)]
        rec = [len(states[t]['RECOVERED']) for t in range(t)]
        vac = [len(states[t]["VACCINATED"]) for t in range(t)]
        qur = [len(states[t]["QUARANTINED"]) for t in range(t)]

        #Plot Graphs
        if show:
            plt.plot(sus, label = "SUSCEPTIBLE")
            plt.plot(exp, label = "EXPOSED", color="orange")
            plt.plot(inf, label = "INFECTED", color="red")
            plt.plot(rec, label = "RECOVERED", color="green")
            plt.plot(vac, label="VACCINATED", color = "Gray")
            plt.plot(qur, label = "QUARANTINED", color = "Purple")
        
            plt.xlabel("Time")
            plt.ylabel("Number of nodes")      
            plt.legend()
            plt.title("Low prob_e - High prob_i")
            plt.savefig("results\\simulation_results\\"+"indiv_sim_prob_e_h.png")
            plt.show()  

        #Return results
        return (sus, exp, inf, rec)

    '''
    This method extracts the number of nodes in each state for every community.

    params:
    t         - number of simulation steps
    states    - dictionary of shape states[ time ][ node_state ], representing the states all nodes are in at time 't'
    partition - dictionary with keys representing the node and 
                values the community each nodes is part of
    show      - True for data to be ploted on graph, False otherwise

    returns:
    list - list of touples containing the nodes in each state for every community. 
    '''
    def simulationCommunityGraphic(self, states, partition, t, show=True):

        #Infection in each comunity at time t
        results = []
        
        for t in range(t):
            
            # Record the number of infected nodes in each community at time t
            inf={part: 0 for part in set(partition.values())}
            exp ={part: 0 for part in set(partition.values())}
            sus={part: 0 for part in set(partition.values())}
            rec ={part: 0 for part in set(partition.values())}
            qur={part: 0 for part in set(partition.values())}
            vac ={part: 0 for part in set(partition.values())}
            
            for node in states[t]["INFECTED"]: 
                inf[partition[node]]+=1
            
            #Record the number of exposed node s in each community at time t
            for node in states[t]["EXPOSED"]:
                exp[partition[node]]+=1

            for node in states[t]["RECOVERED"]:
                sus[partition[node]]+=1

            for node in states[t]["SUSCEPTIBLE"]:
                rec[partition[node]]+=1
            
            for node in states[t]["QUARANTINED"]:
                qur[partition[node]]+=1

            for node in states[t]["VACCINATED"]:
                vac[partition[node]]+=1

            result = (sus, exp, inf, rec, qur, vac)
            results.append(result)
            
        return results

    '''
    This method extract numerical data from the simulation. This method calaculates the
    infectious rate, reproduction number and the most infectious node

    returns:
    touple - touple containing the (infectious rate, index of the most infectious node, reproduction number)
    '''
    def numericalData(self):

        #Infection rate 
        inf_rate = self.infection_rate()

        #Super Spreader
        index = self.superSpreader() 
        
        #Reproduction number
        r0_num = self.reproductionNumber()


        return (inf_rate, index, r0_num) 

    '''
    This method finds the most infectious node in the network, called Super Spreader

    params:
    infector_dict - dictionary that keeps track who infected each node
    
    returns:
    int - the index of the most infection nodes in the {infector_dict}
    '''
    def superSpreader(self, infector_dict):
        scores = [0] * len(infector_dict.keys())
        for node in infector_dict.keys():
            scores[infector_dict[node]] += 1

        
        return self.find_max_value_index(scores)



    '''
    This method calculates the infectious rate values for the range [day1, 't']. The infectious rate is defined as:
    inf_rate = (new recoreded cases) / (total number of nodes that could become infected) * 100

    params:
    t         - number of simulation steps
    states    - dictionary of shape states[ time ][ node_state ], representing the states all nodes are in at time 't'

    list - list of size 't' containing all the infectious rate values for the simulation
    '''
    def infection_rate(self, states, t):
        inf_rate = []
        
        for t in range(1, t):
            #If no new nodes are infected the infection rate is 0
            if len(states[t-1]["INFECTED"]) > len(states[t]["INFECTED"]):
                inf_rate.append(0)
            else:
                rate = len(states[t]["INFECTED"]) - len(states[t-1]["INFECTED"]) 

                #Avoid divison by 0
                try:
                    #Calculate the rate value for time 't'    
                    rate /= len(states[t]["EXPOSED"]) + len(states[t]["SUSCEPTIBLE"]) 
                except ZeroDivisionError as e:
                    pass

            #Add rate value of time t into the inf_rate list
            inf_rate.append(abs(rate)*100)

        return inf_rate   


    '''
    This method caculates the reproduction number (R0). Reproduction number is determined by the probability of exposure, infection, 
    infection duration and the proportion of exposed nodes in the whole network. 

    Rt = ( EXPOSED_NODES / n ) * prob_i * prob_e * infection_duration 

    params:
    t         - number of simulation steps
    states    - dictionary of shape states[ time ][ node_state ], representing the states all nodes are in at time 't'

    returns:
    list - list representing the value of reproduction nuumber for each day in the range(0, 't')
    '''
    def reproductionNumber(self, states, t, show=True):
        R0 = self.prob_i * self.prob_e * self.infection_duration

        Rt = []
        for t in range(t):
            sus_prop = len(states[t]["EXPOSED"]) / len(self.G)
            Rt.append(sus_prop*R0)

        if show:
            plt.plot(Rt, label="effective reproduction number")
            plt.legend()
            plt.show()

        return Rt
    
    '''
    This method calculates the average number of nodes in a particular state for each community. This method is used for 
    the average number of infected nodes, but it can be modifiied to return values for each step. This method should be used to 
    find the avereage number over a scenario, meaning a colection of simulations.

    params:
    t              - number of simulation steps
    community_data - dictonary whose keys are the simulation number and values are the touples returned by 'simulationCommunityGraphic()'

    return:
    avg_inf - average number of nodes in a paprticular state for each community
    '''
    def communityAverage(self, times, community_data):

        total_infected = {part: 0 for part in set(self.partition.values())}
        num_steps = 0
        avg_inf_sim = []
        for time in range(times):
        # Iterate over the results and accumulate the total number of infected nodes for each community
            for result in community_data[time]:
                sus, exp, inf, rec, qur, vac = result
                for part in total_infected:
                    total_infected[part] += inf[part]
                num_steps += 1

            avg_inf_sim.append({part: total_infected[part] / num_steps for part in total_infected})

        avg_inf = {part: sum(avg_inf_sim[s][part] for s in range(times)) / times for part in total_infected}
        return avg_inf

    '''
    This method calculates and displays the statisical data extracted from all the scenarios run. This method finds:
            - average simulation time
            - average infectious rate
            - average reproduction number
            - most common super-spreader
            - the community super-spreader belongs to 
            - average number of infecred nodes in each community

    This method saves the data into a .txt file in 'results\simuluation_results\\numerical_data' folder

    params:
    times - number of simulations to include in the scenario reults
    title - title to use to save the data.
    '''
    
    def scenraio_results(self, times, title):
        
        avg_time = 0
        saveFile = open("results\\simulation_results\\numerical_data\\"+title+".txt", "w")

        for i in range(times):
            avg_time += self.stats_scenario[i][0]

        avg_time/= times
        print("Average Time: ", avg_time)
        saveFile.write("Avereage Time: " + str(avg_time) + "\n")
        avg_inf = []
        for i in range(times):
            avg_inf.append(sum(self.infection_rate(self.stats_scenario[i][1], self.stats_scenario[i][0]))/self.stats_scenario[i][0])
        
        
        saveFile.write("Average Infections Rate: " + str(sum(avg_inf) / times) + "\n")
        print("Average Infections Rate: ",sum(avg_inf) / times)

        avg_r0 = []
        for i in range(times):
            avg_r0.append(sum(self.reproductionNumber(self.stats_scenario[i][1], self.stats_scenario[i][0], False)))

        saveFile.write("AverageR0: " + str(sum(avg_r0) / times) + "\n")
        print("Average R0: ", sum(avg_r0) /times)

      

        nodes_list = []
        for i in range(times):
            nodes_list.append(self.superSpreader(self.stats_scenario[i][2]))
        superSpreader = nodes_list[self.find_max_value_index(nodes_list)]
        
        saveFile.write("Super-Spreader :" + str(superSpreader) + "\n")
        print("Super Spreader: ", superSpreader)
        
        saveFile.write("Community of SuperSpreader: " + str( self.partition[superSpreader])  + "\n")
        print("Community of SuperSpreader: ", self.partition[superSpreader])

        community_data = []
        for i in range(times):
            community_data.append(self.simulationCommunityGraphic(self.stats_scenario[i][1], self.partition, self.stats_scenario[i][0], False))

        print("\nAverage number of nodes by Communiteis:")
        saveFile.write("\nAverage number of nodes by Communiteis:\n")
        inf_avg = self.communityAverage(times, community_data)
        
        for k in inf_avg.keys():
            saveFile.write("\tCommunity" + str(k) + " : " + str(inf_avg[k]) + "\n")
            print("\tCommunity", str(k), " : ", str(inf_avg[k]))

        plt.bar([str(i) for i in inf_avg.keys()], inf_avg.values(), align='center')
        plt.title(title)
        plt.xlabel("Community label")
        plt.ylabel("% of infected nodes")
        plt.savefig(".\\results\\simulation_results\\"+title+".png")
        plt.show()

        saveFile.close()
        

        
    #Helper function that finds the index of the max value in list [scores]
    def find_max_value_index(self, scores):
        max = -1
        index = -1
        for i in range(len(scores)):
            if scores[i] > max:
                max = scores[i]
                index = i
        return index
        
