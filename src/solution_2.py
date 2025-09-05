import gurobipy as gp
import data as dt
import numpy as np


# data class 
data = dt.Data("C:\\Users\\jackg\\OneDrive\\Documents\\MATH3205\\Team Project\\2023.0223\\data\\instances\\district_06-01.txt")

# preprocessing
print("starting preprocessing")
patterns_of_size = {i:set() for i in range(1,9)}
patterns_of_size[1].add(frozenset([(0,0),]) )
for size in range(2,9):
    for pattern in patterns_of_size[size-1]:
        for x,y in pattern:
            for dx,dy in [(1,0),(0,1),(-1,0),(0,-1)]:
                if (x+dx,y+dy) not in pattern:
                    patterns_of_size[size].add(pattern.union(frozenset([(x+dx,y+dy),])))
# patterns from root of size 8 that intersect point
patterns_hit_point8 = {}
patterns_hit_point6 = {}
for pattern in patterns_of_size[8]:
    for x,y in pattern:
        if (x,y) not in patterns_hit_point8:
            patterns_hit_point8[(x,y)] = set()
        patterns_hit_point8[(x,y)].add(pattern)
        
for pattern in patterns_of_size[6]:
    for x,y in pattern:
        if (x,y) not in patterns_hit_point6:
            patterns_hit_point6[(x,y)] = set()
        patterns_hit_point6[(x,y)].add(pattern)
        
# functions                          
def get_patterns(root,point,size):
    """ gets a set of patterns from a generated from a root point
    that intersect a particular point"""
    ref_point = point - root
    if size == 6:
        return patterns_hit_point6[ref_point]
    elif size == 8:
        return patterns_hit_point8[ref_point]
    
def cost_of_patterns(patterns):
    """ a fixed cost of choosing a pattern from some combinor"""
    cost_of_patterns = {}
    for pattern in patterns:
        cost_of_patterns[pattern] = sum(sum(abs(point)*data.edge_len) for point in pattern)
    return cost_of_patterns
    
def discard_patterns(patterns:set,invalid_points:frozenset):
    """ discards patterns that intersect with invalid points"""
    for pattern in patterns:
        if not invalid_points.isdisjoint(pattern):
            patterns.discard(pattern)
                    
def Neighbors(v,data: dt.Data,V:set):
    """ get neighbors of nodes every node v """
    VNeighbors = []
    for v in V:
        for a in [[-1, 0], [1, 0], [0, 1], [0, -1]]:
            loc = (data.location[v] + np.array(a)) % np.array([data.row_num, data.col_num])
            if data.graph[tuple(loc)] in V:
                VNeighbors.append(data.graph[tuple(loc)])
    return VNeighbors

    
#Sets and Subsets
V = range(len(data.location)) # set of all nodes for PVA`s, Combiners , Inverters and service paths
VL = data.way_set # subset of V for each service way l, VL[l] for l in LD = range(10)
VD = data.district # subset of V for each district d , VD[d] for d in D
VN = [Neighbors(v,data,set(V)) for v in V] # subset of V for each neighbor v, VN[v]
D = range(data.district_num) # set of all districts 
L = range(data.way_num) # set of all candidate service paths 
LD = [[l for l in range(data.way_num) if data.way_cover[l][d] == 1] for d in D] # set of service paths that pass through district l
T = range(0,2) # set of combiner box types 

# Data
path_cost_per_meter = 46.312
V_dist = data.distance # distance between v1 in V and v2 in V so V_dist[v1][v2]
V_loc = data.location # coords for v in V so V_loc[v]
loc_V = {V_loc[v]:v for v in V}

L_cost_path = [length*path_cost_per_meter for length in data.way_len] # cost of each service path l in L
IC_cable_cost = 10.904 # cable cost for inverter to combiner per meter of cable
CP_cable_cost = 3.937 # cable cost for combiner to PVA per meter of cable
CT_capacity = [6,8] # capacity of combiner Box of type t in T
CT_cost = [94.476,110.222] # cost of combiner box of type t in T
Cost_of_patterns = cost_of_patterns(patterns_of_size[6]) | cost_of_patterns(patterns_of_size[8]) # cost of choosing pattern p for 
# some combiner box, Cost_of_patterns[p]
"""BMP1 benders master problem  is where do we put the the service paths
such that each service path passes through a district """



"""
BMP = gp.Model()
Z = {l:BMP.addVar(vtype = gp.GRB.BINARY) for l in L}
B = {(l,d): BMP.addVar(vtype= gp.GRB.BINARY) for d in D for l in LD[d]}
Theta = {(l,d): BMP.addVar(lb = 0,vtype= gp.GRB.CONTINUOUS) for d in D for l in LD[d]}
# cost of service paths chosen plus lowerbound of subproblem cost
BMP.setObjective(gp.quicksum(L_cost_path[l]*Z[l] for l in L) + gp.quicksum(Theta[l,d] for d in D for l in LD[d]),gp.GRB.MINIMIZE)
# for each district look at service paths that pass through that district and choose 1 of those service paths to put the invertor
choose_service_path_in_D = {d:BMP.addConstr(gp.quicksum(B[l,d] for l in LD[d]) == 1)  for d in D }
# if you choose some service path in district k then you must choose that service path overall
link_Z_and_B = {(l,d):BMP.addConstr(Z[l] >= B[l,d]) for d in D for l in LD[d]}
# for each district there must be some service path that passes through it
choose_service_path_in_D = {d:BMP.addConstr(gp.quicksum(Z[l] for l in LD[d]) >= 1)  for d in D }
"""
def subproblem(nodes,inverter,service_path,combiners):
    """
    nodes is the set of pva`s in the district for this subproblem with the inverter and
    service path taken out
    combiners a set containing type and fixed position
    inverter a fixed position
    service path a set of fixed positions
    
    so this problem is basicly solving how to assign pva`s to fixed combinor boxes 
    given a fixed inverter and fixed service path. this is done by generating all 
    patterns from each combiner box such that a pva is overlapped by exactly 1 pattern.
    """
    SP = gp.Model()
    X = {(c,p): SP.addVar(vtype = gp.GRB.BINARY) for c in combiners for p in patterns_of_size[c[0]]}
    SP.setObjective(gp.quicksum(Cost_of_patterns[p]*X[c,p] for c in combiners for p in patterns_of_size[c[0]]))
    
    Node_coverage = {node:SP.addConstr(
                    gp.quicksum(
                    gp.quicksum(X[c,p] for p in get_patterns(c[1],V_loc[node],c[0])) 
                    for c in combiners) == 1) 
                    for node in nodes}
    No_combiner_overlap = {c1:SP.addConstr(
                        gp.quicksum(X[c1,p] for c2 in combiners for p in get_patterns(c1[1],c2[1],c1[0]) if c1 != c2) == 0) 
                        for c1 in combiners}
    No_inverter_overlap = {c:SP.addConstr(
                        gp.quicksum(X[c,p] for p in get_patterns(c[1],V_loc[inverter],c[0])) == 0) 
                        for c in combiners}
    No_service_path_overlap = {(c,node):SP.addConstr(
                        gp.quicksum(X[c,p] for p in get_patterns(c[1],V_loc[node],c[0])) == 0) 
                        for c in combiners for node in service_path} 
    
        







