# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:23:36 2019

@author: QUA042
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

# Status=0     => Susceptible
# Status=1-5   => Exposed
# Status=10    => Reported Infectious
# Status=-1    => Killed by management

size=10000   #size of the World
N=100      # number total of farms
i=1        # number of infected farms at the beginning, (N-i) are susceptible   

A = 100   #The world is divided into A² cells
a = size/A  #Each cell is a square of a side = size / A

#parameters : we choose these parameters from the article (Sellman) in the case of sheeps
trans_1=0.00083
trans_2=0.49
susc_1=1
susc_2=0.2
alpha=1
beta=1600
gamma=4.6
radius = 20
radius_2 = 20
radius_3 = 100
max_animals = 300   # max_animals, used in case of a random world generated
budget = 3      # daily budget of 3 farms

import numpy as np
from math import sqrt
import math
import pylab
import random
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import copy


# Each farm is defined by its position in the World, the number of animals it contains and its status    
class farm:
    def __init__(self, x, y, animals, status):
        self.x = x
        self.y = y
        self.animals = animals
        self.status = status
    
    def Become_Exposed(self):
        self.status=1
        
    def Transmissibility(self):
        return trans_1*(self.animals**trans_2)
    
    def Susceptibility(self):
        return susc_1*(self.animals**susc_2)
    
# Calculate the distance between two farms
def d(f1,f2):
    return sqrt((f1.x-f2.x)**2+(f1.y-f2.y)**2)
    
# K function, used in the disease model
def K(d):
    return alpha/(1+((d/beta)**gamma))

# Function that kills a farm (index of the farm in the World is taken as a parameter)
def Kill(W, elem):
    W[elem].status=-1
    return W
    
# Return the distance to the nearest farm to a given farm f and a list of farms cell. Useful especially when cell doesn't contain f       
def distance_nearest_farm(f, cell):
    m = len(cell)
    min=d(f, cell[0])
    for k in range(m):
        if d(f, cell[k])<min:
            min=d(f, cell[k])
    return min
    

########### These functions were used to try to implement the Conditional Entry Algorithm ##############
# Function that determines the number of the cell than contains a given farm
def cell_of_farm(farm):
    k = np.floor(farm.x*A/size)
    p = np.floor(farm.y*A/size)
    return p*A+k+1

# Function that tranforms the World in a list of its cells
def World_in_cells(World):
    W_cells=[]
    for k in range(1,A*A+1):
        cell_k = []
        for f in World:
            if cell_of_farm(f)==k:
                cell_k.append(f)
        W_cells.append(cell_k)
    return W_cells



# Probability that at least one of the susceptible farms in a cell become infected by an infectious farm
def v(farm, cell):
    if cell==[]:
        return 0
    else:
        T_i = farm.Transmissibility()
        S_max=0
        for f in cell:
            if f.Susceptibility() > S_max :
                S_max = f.Susceptibility()
        return 1-math.exp(-T_i*S_max*K(distance_nearest_farm(farm, cell)))
    

# Probability that at least one of the susceptible farms in a cell become infected by an infectious farm
def w(farm, cell):
    J=0
    for f in cell:
        if f.status==0:
            J+=1
    return 1-(1-v(farm,cell))**J


# Probability that an infectious farm (1) infects a susceptible one (2)
def p(farm_1,farm_2):
    return (1-math.exp(-farm_1.Transmissibility()*farm_2.Susceptibility()*K(d(farm_1,farm_2))))

# Function that displays the state of each farm of the World
def aff(W):
    L=[]
    for f in W:
        L.append(f.status)
    return L

# Return the index of a given farm in a World 
def index(farm, W):
    p=0
    for k in range(len (W)):
        if W[k]==farm:
            p=k
    return p

def takeSecond(elem):
    return elem[1]

# Function that removes an element from a list
def remove(list, element):
    new_list=[]
    for e in list:
        if e!=element:
            new_list.append(e)
    return new_list


## To create a World, we have different functions : a Random_World function that gives a set of farms with random positions
## For our simulation, we use a World extracted from a real data from a csv file
    
# Generate randomly a World of farms with N farms (N-i are susceptible, i are infectious)
def Random_World(N,i):
    np.random.seed(1)
    randoms_x=size * np.random.random_sample((N,))
    np.random.seed(5)
    randoms_y=size * np.random.random_sample((N,))
    Random_World=[]
    for k in range(N-i):
        x_k = randoms_x[k]
        y_k = randoms_y[k]
        random_number = max_animals
        f=farm(x_k,y_k,random_number,0)
        Random_World.append(f)
    for k in range(N-i,N):
        x_k = randoms_x[k]
        y_k = randoms_y[k]
        rand_number = max_animals
        g=farm(x_k,y_k,rand_number,1)
        Random_World.append(g)
    return Random_World

# Construction of the World used for the simulation :
data = np.loadtxt('Data.csv', delimiter=",")
total_data=[]
total_size=len(data)
for ff in range(total_size):
    f = farm(data[ff][2],data[ff][3],data[ff][4],0)
    total_data.append(f)

# In the simulation, we build W which is a part of the data total_data
# We choose a part in 'the middle' of the World
W=[]
for f in total_data:
    if 650000<f.x<660000 and 750000<f.y<760000:
        W.append(f)
size_world=len(W)
# Then, we choose randomly one farm that is going to be infectious
infectious_farm = random.randint(0,size_world)
W[infectious_farm].status=1
    

# Calculate the probability that a farm become infectious by a set of farms 
def ProbaInfection(farm, set):
    product=1
    for infecting in set:
        product = product * (1-p(farm, infecting))
    return 1-product

# Function that spreads the disease among all the farms of W using the REAL state of the World
def Simple_Algo(W):
    infectious = Infectious(W)
    susceptibles = Susceptible_Farms(W)
    for i in infectious:
        for j in susceptibles:
            random_ij = np.random.random()
            if random_ij < p(i,j):
                j.status=1
    return W

# Function that spreads the disease among all the farms of W using the REPORTED state of the World
def Algo_Reported(W):
    reported_infectious = Infectious_Reported(W)
    susceptibles = Susceptible_Farms(W)
    for i in reported_infectious:
        for j in susceptibles:
            random_ij = np.random.random()
            if random_ij < p(i,j):
                j.status=1
    return W


# Returns the farms that are around a given farm and diven radius. We choose only the non-killed farms (status>=0)
def Alive_Ring(W, farm, radius):
    ring=[]
    for f in W:
        if (d(f, farm)<radius and farm.status>=0):
            ring.append(f)
    return ring


#############################################   MANAGEMENT STRATEGIES   #############################################################


# Update the World without any management
def Iteration_0(W):
# We increment status of all infectious
    for f in W:
        if f.status>0:
            f.status=f.status+1
    infectious = Infectious(W)
    susceptibles = Susceptible_Farms(W)
    for i in infectious:
        for j in susceptibles:
            random_ij = np.random.random()
            if random_ij < p(i,j):
                j.status=1
    return W    

## Strategy 1 :
    
def Iteration_1(W):
    reported_infectious = Infectious_Reported(W)
    remaining_budget = budget
    
    # We increment status of all infectious
    for f in W:
        if f.status>0:
            f.status=f.status+1
                
    # We update the World with the Conditional Entry Algorithm that spreads the disease
    W=Simple_Algo(W)
    
    # We kill firstly all the reported infectious
    for inf in reported_infectious:
        if remaining_budget>0:
            inf.status=-1
            remaining_budget=remaining_budget-1
        

    if remaining_budget==0:
        return W
    
    else:
    # If we have left budget after killing all reported infectious, we kill alive farms by ranking them with probability of infection
        alive = Alive_Farms(W)
        Probas=[]
        for f in alive :
            proba_f = ProbaInfection(f, reported_infectious)
            Probas.append((f, proba_f))
            Probas_sorted = sorted(Probas, key=takeSecond, reverse=True)
        if len(reported_infectious)>0:
            for k in range(remaining_budget):
                if (remaining_budget>0 and len(Probas_sorted)>remaining_budget):
                    Probas_sorted[k][0].status=-1
                    remaining_budget=remaining_budget-1
    return W




    
## Strategy 3 :
    
def Iteration_3(W, reported_killed):
    reported_infectious = Infectious_Reported(W)
    remaining_budget = budget
    updated_reported_killed = reported_killed
    # We increment status of all infectious
    for f in W:
        if f.status>0:
            f.status=f.status+1
                
    # We update the World
    W=Simple_Algo(W)

    for inf in reported_infectious:
        if remaining_budget>0:
            inf.status=-1
            remaining_budget=remaining_budget-1
            updated_reported_killed.append(inf)
    
            
    size = len(updated_reported_killed)
    if (size>0):
        for j in range(size):
            reported_j = updated_reported_killed[size-j-1]
            ring_j = Alive_Ring(W, reported_j, radius)   # ring contains all the alive farms around last_killed
            for farm_ring in ring_j:
                if (remaining_budget>0):
                    farm_ring.status=-1
                    remaining_budget=remaining_budget-1
    return W, updated_reported_killed

    

## Strategy 2 :
    
# In this strategy, we need to rank susceptibles upon their propensity to inflict further damage. To this purpose, we consider
# firstly all susceptibles that are close to an infectious farm. We determine it by this function :
def CloserSusceptibles(W, radius_2):
    CloserSusceptibles=[]
    infectious = Infectious(W)
    susceptibles = Susceptible_Farms(W)
    for inf in infectious:
        for susc in susceptibles:
            if d(inf, susc)<radius_2:
                CloserSusceptibles.append(susc)
    return CloserSusceptibles
    
# After the choice of a susceptible, we simulate the state of the world during T timesteps, with this susceptible and without it.
# At each timestep, we kill all reported infectious farms. We calculate then the difference between number of killed farms in the
# two cases :


def Simulate_World(W, sus, T):
    W_with_s = copy.deepcopy(W)
    W_without_s = copy.deepcopy(W)
    W_without_s[sus].status=-1
    
    for k in range(T):
        W_with_s = Simple_Algo(W_with_s)
    
    for k in range(T):
        W_without_s = Simple_Algo(W_without_s)
    
    return 1+Nb_Culled_Farms(W_with_s)-Nb_Culled_Farms(W_without_s)


def Rank_Susceptibles(W, T, radius_2):
    closer = CloserSusceptibles(W, radius_2)
    new_closer = []
    for susc in closer:
        ind = index(susc, W)
        difference = Simulate_World(W, ind, T)
        new_closer.append((susc, difference))
    closer_ranked = sorted(new_closer, key=takeSecond, reverse=True)
    return closer_ranked

    
def Iteration_2(W, T):
    reported_infectious = Infectious_Reported(W)
    remaining_budget = budget
    
# We increment status of all infectious
    for f in W:
        if f.status>0:
            f.status=f.status+1
                
# We update the World with the Conditional Entry Algorithm that spreads the disease
    W=Simple_Algo(W)
# We kill firstly all the reported infectious
    for inf in reported_infectious:
        if remaining_budget>0:
            inf.status=-1
            remaining_budget=remaining_budget-1
       
    left_budget = remaining_budget
    
    if(Suc_Farms(W)==0):
        return W
    
# If we have left budget after killing all reported infectious :
    elif(Suc_Farms(W)>0 and remaining_budget>0):
        ranking = Rank_Susceptibles(W, T, radius_2)
        for i in range(left_budget):
            if (remaining_budget>0):
                ranking[i][0].status=-1
                remaining_budget=remaining_budget-1
    return W



##################################################   GRAPHICS   #############################################################

# Function that draw the whole World W
def Draw_World(World):
    X1_draw=[] #susceptibles
    Y1_draw=[]
    X2_draw=[] # infectious
    Y2_draw=[]
    X3_draw=[] # killed
    Y3_draw=[]
    for farm in World:
        if farm.status==0:
            X1_draw.append(farm.x)
            Y1_draw.append(farm.y)
            plt.scatter(X1_draw,Y1_draw, c='green', s=5)
        elif farm.status>0:
            X2_draw.append(farm.x)
            Y2_draw.append(farm.y)
            plt.scatter(X2_draw,Y2_draw, c='red', s=5)
        elif farm.status==-1:
            X3_draw.append(farm.x)
            Y3_draw.append(farm.y)
            plt.scatter(X3_draw,Y3_draw, c='black', s=5)


# Function that draw the susceptible farms
def Draw_Susceptibles(World):
    X_s=[]
    Y_s=[]
    for farm in World:
        if farm.status==0:
            X_s.append(farm.x)
            Y_s.append(farm.y)
        plt.scatter(X_s,Y_s, c='chartreuse', s=5)


# Function that draw the infectious farms
def Draw_Infectious(World):
    X_i=[]
    Y_i=[]
    for farm in World:
        if 0<farm.status<13:
            X_i.append(farm.x)
            Y_i.append(farm.y)
        plt.scatter(X_i,Y_i, c='red', s=5)

# Function that draw the killed farms
def Draw_Culled(World):
    X_c=[]
    Y_c=[]
    for farm in World:
        if farm.status==-1:
            X_c.append(farm.x)
            Y_c.append(farm.y)
        plt.scatter(X_c,Y_c, c='black', s=5)

# Function that returns a list of killed farms
def Culled_Farms(World):
    c=[]
    for f in World:
        if f.status==-1:
            c.append(f)
    return c

# Function that gives the number of killed farms
def Nb_Culled_Farms(World):
    c=0
    for f in World:
        if f.status==-1:
            c+=1
    return c

# Function that gives the number of killed animals
def Killed_Animals(World):
    total_killed=0
    culled = Culled_Farms(World)
    for c in culled:
        total_killed = total_killed + c.animals
    return total_killed

# Function that gives the number of kinfectious farms
def Infectious_Farms(World):
    i=0
    for f in World:
        if f.status>0:
            i+=1
    return i

# Function that gives the number of susceptible farms
def Suc_Farms(World):
    s=0
    for f in World:
        if f.status==0:
            s+=1
    return s

# Function that gives the number of kreported farms
def Reported_Infectious_Farms(World):
    i=0
    for f in World:
        if f.status>9:
            i+=1
    return i

# Function that returns a list of infectious farms
def Infectious(World):
    S=[]
    for f in World:
        if f.status>0:
            S.append(f)
    return S

# Function that returns a list of reported infectious farms
def Infectious_Reported(World):
    S=[]
    for f in World:
        if f.status>9:
            S.append(f)
    return S

# Function that returns a list of infectious farms
def Susceptible_Farms(W):
    size=len(W)
    S=[]
    for k in range(size):
        if W[k].status==0:
            S.append(W[k])
    return S   

# Function that returns a list of alive farms
def Alive_Farms(W):
    size=len(W)
    S=[]
    for k in range(size):
        if W[k].status>=0:
            S.append(W[k])
    return S  



strategy=3 

Initial_World = W
Current_World = Initial_World

#def Iteration(W, strategy):
#    if strategy==1:
#        return Iteration_1(W)
#    elif strategy==2:
#        return Iteration_2(W, 5)
#    else:
#        return Iteration_3(W)


def Iterate_World(strategy, p):
    Current_World = Initial_World
    if p==0: 
        return Initial_World
    else: 
        if strategy==1:
            Current_World = Iteration_1(Current_World)
        elif strategy==2:
            Current_World = Iteration_2(Current_World)
        else:
            Current_World = Iteration_3(Current_World)
    return Current_World
 

# Number of animals killed per day :
def Several_Worlds(N_Sim, strategy):
    Initial_World=W
    List_of_lists=[]
    for k in range(N_Sim):
        World_k= copy.deepcopy(Initial_World)  # We create a copy of the initial World for the simulation
        list_of_killed_k=[]
        day_k=0
        inf_k=Infectious_Farms(World_k)
        rep_k=Reported_Infectious_Farms(World_k)
        reported_killed=[]
        while(rep_k==0):   # Loop while to initialize the World; we start the management at the 1st infectious reported
            day_k = day_k+1
            World_k=Iteration_0(World_k)
            list_of_killed_k.append([day_k, 0])
            inf_k=Infectious_Farms(World_k)              # update the value of the infectious farms
            rep_k=Reported_Infectious_Farms(World_k)     # update the value of the reported infectious farms
        while(inf_k>0):
            day_k = day_k+1
            World_k, reported_killed=Iteration_3(World_k, reported_killed)
            inf_k=Infectious_Farms(World_k)
            rep_k=Reported_Infectious_Farms(World_k)
            Killed_k = Nb_Culled_Farms(World_k)
            list_of_killed_k.append([day_k, Killed_k])
#            print(day_k)
#            print(Infectious_Farms(World_k))
        List_of_lists.append(list_of_killed_k)
    return List_of_lists

def ListToArray(list):
    width=len(list)
    # length is the maximum of the set {len(list[k])}
    length=len(list[0])
    for p in range(len(list)):
        if length<len(list[p]):
            length=len(list[p])
    array=np.zeros((width, length))
    for i in range(width):
        for j in range(len(list[i])):
            array[i][j]=list[i][j][1]
        for j in range(len(list[i]),length):
            array[i][j]=list[i][len(list[i])-1][1]
    return array

# Function that gives 
def T_end(W1, strategy):
    World= copy.deepcopy(W1)   # We create a copy of the initial World for the simulation
    inf=Infectious_Farms(World)
    day=0
    reported_killed=[]
    while(inf>0):
        World, reported_killed=Iteration_3(World, reported_killed)
        inf=Infectious_Farms(World)     # update the value of the infectious farms
        day = day+1
    return day
 
def T_ends(N_Sim, strategy):
    T_ends=[]
    for k in range(N_Sim):
        W_k=copy.deepcopy(W)
        T_end_k = T_end(W_k, strategy)
        T_ends.append(T_end_k)
    return T_ends
        
N2_Sim=100
L1 = T_ends(N2_Sim, 1)
L3 = T_ends(N2_Sim, 3)
plt.scatter(['Strategy 1']*N2_Sim, L1, c='blue', s=5)
plt.scatter(['Strategy 3']*N2_Sim, L3, c='red', s=5)
plt.title('T_end values in days using ' + str(N2_Sim) + ' simulations')






N_Sim=100

#Plot of strategy 1
L_1=Several_Worlds(N_Sim, 1)
array_1 = ListToArray(L_1)
array_to_plot_1 = np.zeros(len(array_1[0]))
for m in range(len(array_1[0])):
    Sum_m=0
    for j in range(len(array_1)):
        Sum_m=Sum_m + array_1[j][m]
    array_to_plot_1[m] = Sum_m / len(array_1)
pylab.plot(array_to_plot_1, label='Strategy 1')
pylab.legend(loc='upper left')
        
#Plot of strategy 3
L_3=Several_Worlds(N_Sim, 3)
array_3 = ListToArray(L_3)
array_to_plot_3 = np.zeros(len(array_3[0]))
for m in range(len(array_3[0])):
    Sum_m=0
    for j in range(len(array_3)):
        Sum_m=Sum_m + array_3[j][m]
    array_to_plot_3[m] = Sum_m / len(array_3)
pylab.plot(array_to_plot_3, label='Strategy 3')
pylab.legend(loc='upper left')


pylab.xlabel('Duration (days)')
pylab.ylabel('Cumulated number of killed farms')
#



