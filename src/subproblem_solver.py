import gurobipy as gp
import numpy as np
from gurobipy import GRB
import math
import time
import random as rd

import data as dt


class DistrictSolver:
    def __init__(self, data: dt.Data, district_id: int):
        self.data = data
        self.district = district_id
        self.nodes = data.district[district_id]
        self.node_num = len(self.nodes)
        self.way = [i for i in range(data.way_num) if data.way_cover[i][district_id] == 1]
        self.lb = {node: None for node in self.nodes}
        self.cover_range = {
            i: [j for j in self.nodes if sum(abs(data.location[i] - data.location[j])) < max(data.capacity)]
            for i in self.nodes
        }
        self.adjacent = {i: [] for i in self.nodes}
        for i in self.nodes:
            for a in [[-1, 0], [1, 0], [0, 1], [0, -1]]:
                loc = (data.location[i] + np.array(a)) % np.array([data.row_num, data.col_num])
                if data.graph[tuple(loc)] in self.nodes:
                    self.adjacent[i].append(data.graph[tuple(loc)])

    def __relaxation(self, inverter: int):
        nodes = list(self.nodes)
        nodes.remove(inverter)
        cover_range = {node: self.cover_range[node][:] for node in nodes}
        for node in nodes:
            if inverter in cover_range[node]:
                cover_range[node].remove(inverter)
        types = self.data.capacity

        model = gp.Model()
        combiner = model.addVars(nodes, types, vtype=GRB.CONTINUOUS, ub=1.0)
        secondary = model.addVars(nodes, nodes, vtype=GRB.CONTINUOUS, ub=1.0)

        model.addConstrs(
            secondary[i, i] == gp.quicksum(combiner[i, t] for t in self.data.capacity) for i in nodes
        )
        model.addConstrs(gp.quicksum(secondary[i, j] for j in nodes) == 1 for i in nodes)
        model.addConstrs(
            gp.quicksum(secondary[j, i] for j in nodes)
            <= gp.quicksum(t * combiner[i, t] for t in self.data.capacity)
            for i in nodes
        )
        model.addConstrs(
            secondary[i, j] == 0 for i in nodes for j in nodes if j not in cover_range[i]
        )
        model.addConstrs(
            secondary[j, i] <= secondary[i, i] for i in nodes for j in nodes
        )
        model.addConstr(
            gp.quicksum(t * combiner[i, t] for t in self.data.capacity for i in nodes) >= len(nodes)
        )
        model.addConstr(
            gp.quicksum(combiner[i, t] for i in nodes for t in self.data.capacity)
            >= math.ceil(len(nodes) / max(self.data.capacity))
        )

        p_cost = self.data.primary_cost * gp.quicksum(
            self.data.distance[i, inverter] * gp.quicksum(combiner[i, t] for t in types)
            for i in nodes
        )
        s_cost = self.data.secondary_cost * gp.quicksum(
            self.data.distance[i, j] * secondary[i, j] for i in nodes for j in cover_range[i]
        )
        fix_cost = gp.quicksum(
            gp.quicksum(self.data.price[t] * combiner[i, t] for t in self.data.capacity)
            for i in nodes
        )

        model.Params.OutputFlag = 0
        model.setObjective(p_cost + s_cost + fix_cost, GRB.MINIMIZE)
        model.optimize()

        return model.objVal, model.Runtime

    def relax_all(self):
        total_time = 0
        for node in self.nodes:
            obj, runtime = self.__relaxation(node)
            self.lb[node] = {'value': obj, 'type': 'relaxation', 'solution': None}
            total_time += runtime
        return total_time

    def __LiftedRounding(self, sub_com, sub_node, model):
        cap = self.data.capacity
        total_nodes = len(sub_node)
        remaining_flow = total_nodes - sum(
            model.cbGetNodeRel(model._secondary[i, j]) for i in sub_node for j in sub_com
        )
        facility = [
            sum(model.cbGetNodeRel(model._combiner[j, t]) for j in sub_com)
            for t in cap
        ]

        delta = 0
        for j, capacity in enumerate(cap):
            remainder = total_nodes % capacity
            if remainder == 0:
                continue

            coef = []
            rhs = math.ceil(total_nodes / capacity) * remainder
            lhs = remaining_flow

            for i in range(j + 1):
                c = min(cap[i], remainder)
                coef.append(c)
                lhs += facility[i] * c

            for i in range(j + 1, len(cap)):
                c = (
                    math.floor(cap[i] / capacity) * remainder
                    + min(remainder, cap[i] % capacity)
                )
                coef.append(c)
                lhs += facility[i] * c

            if rhs - lhs > 1e-5:
                model.cbCut(
                    total_nodes
                    - gp.quicksum(model._secondary[i, j] for i in sub_node for j in sub_com)
                    + gp.quicksum(
                        coef[i] * gp.quicksum(model._combiner[j, cap[i]] for j in sub_com)
                        for i in range(len(cap))
                    )
                    >= rhs
                )
            delta += rhs - lhs

        return delta

    def __Callback(self, model, where):
        if where == GRB.Callback.MIPNODE:
            if model.cbGet(GRB.Callback.MIPNODE_OBJBND) > model._upper_bound:
                model.terminate()
                return

            if model.cbGet(GRB.Callback.MIPNODE_NODCNT) == 0:
                tmp_count = len(model._vub)
                relax_com, relax_node = {}, {}

                for c in model._nodes:
                    for i in model._cover_range[c]:
                        if model.cbGetNodeRel(model._secondary[i, c]) >= 1e-5:
                            if model.cbGetNodeRel(model._secondary[i, c]) > model.cbGetNodeRel(
                                model._secondary[c, c]
                            ):
                                model.cbCut(model._secondary[i, c] <= model._secondary[c, c])
                                model._vub.add((i, c))
                            relax_com.setdefault(c, []).append(i)
                            relax_node.setdefault(i, []).append(c)

                if len(model._vub) - tmp_count > 0:
                    return

                visited = {c: False for c in relax_com}
                while sum(visited.values()) < len(visited):
                    com = rd.choice([c for c in visited if not visited[c]])
                    visited[com] = True
                    sub_com, sub_node = {com}, set(relax_com[com])

                    while True:
                        candidates = {
                            c
                            for j in sub_node
                            for c in relax_node[j]
                            if not visited[c]
                        }
                        if not candidates:
                            break

                        pos, max_delta = None, -1e5
                        for candidate in candidates:
                            delta = self.__LiftedRounding(
                                sub_com | {candidate},
                                sub_node | set(relax_com[candidate]),
                                model,
                            )
                            if delta > max_delta:
                                pos, max_delta = candidate, delta

                        visited[pos] = True
                        sub_com.add(pos)
                        sub_node.update(relax_com[pos])

        elif where == GRB.Callback.MIPSOL:
            group = {}
            for c in model._nodes:
                for i in model._cover_range[c]:
                    if model.cbGetSolution(model._secondary[i, c]) >= 1e-5:
                        if c not in group:
                            group[c] = []
                        group[c].append(i)

            for com in group:
                major = [com]
                for i in major:
                    major += [j for j in model._adjacent[i]
                                if j in group[com] and j not in major]
                major_cut = []
                for i in major:
                    major_cut += [j for j in model._adjacent[i]
                                    if j not in major and j not in major_cut]
                rest = [i for i in group[com] if i not in major]
                while len(rest) > 0:
                    l = [rest.pop()]
                    for i in l:
                        for j in model._adjacent[i]:
                            if j in rest:
                                l.append(j)
                                rest.remove(j)
                    l_cut = []
                    for i in l:
                        l_cut += [j for j in model._adjacent[i]
                                    if j not in l and j not in l_cut]
                    for i in l:
                        model.cbLazy(gp.quicksum(
                            model._secondary[k, com] for k in l_cut) >= model._secondary[i, com])
                        model.cbLazy(gp.quicksum(
                            model._secondary[k, com] for k in major_cut) >= model._secondary[i, com])

    def exact_solve(self, inverter: int, upper_bound: float):
        nodes = list(self.nodes)
        nodes.remove(inverter)

        adjacent = {node: self.adjacent[node][:] for node in nodes}
        for node in nodes:
            if inverter in adjacent[node]:
                adjacent[node].remove(inverter)

        cover_range = {node: self.cover_range[node][:] for node in nodes}
        for node in nodes:
            if inverter in cover_range[node]:
                cover_range[node].remove(inverter)

        if self.lb[inverter]["type"] == "cut-off":
            model = self.lb[inverter]["model"]
        else:
            model = gp.Model()
            combiner = model.addVars(nodes, self.data.capacity, vtype=GRB.BINARY)
            secondary = model.addVars(nodes, nodes, vtype=GRB.BINARY)

            model.addConstrs(
                secondary[i, i] == gp.quicksum(combiner[i, t] for t in self.data.capacity)
                for i in nodes
            )
            model.addConstrs(
                gp.quicksum(secondary[i, j] for j in nodes) == 1 for i in nodes
            )
            model.addConstrs(
                gp.quicksum(secondary[j, i] for j in nodes)
                <= gp.quicksum(t * combiner[i, t] for t in self.data.capacity)
                for i in nodes
            )
            model.addConstrs(
                secondary[i, j] == 0 for i in nodes for j in nodes if j not in cover_range[i]
            )
            model.addConstr(
                gp.quicksum(
                    t * combiner[i, t] for t in self.data.capacity for i in nodes
                )
                >= len(nodes)
            )
            model.addConstr(
                gp.quicksum(combiner[i, t] for i in nodes for t in self.data.capacity)
                >= math.ceil(len(nodes) / max(self.data.capacity))
            )

            p_cost = self.data.primary_cost * gp.quicksum(
                self.data.distance[i, inverter]
                * gp.quicksum(combiner[i, t] for t in self.data.capacity)
                for i in nodes
            )
            s_cost = self.data.secondary_cost * gp.quicksum(
                self.data.distance[i, j] * secondary[i, j]
                for i in nodes
                for j in nodes
            )
            fix_cost = gp.quicksum(
                gp.quicksum(self.data.price[t] * combiner[i, t] for t in self.data.capacity)
                for i in nodes
            )
            model.setObjective(p_cost + s_cost + fix_cost, GRB.MINIMIZE)

            model.Params.OutputFlag = 0
            model._nodes = nodes
            model._cover_range = cover_range
            model._adjacent = adjacent
            model._combiner = combiner
            model._secondary = secondary
            model._vub = set()
            model._upper_bound = upper_bound

            for i in nodes:
                secondary[i, i].BranchPriority = 9
            combiner.BranchPriority = 5
            model.Params.lazyConstraints = 1

        start_time = time.time()
        model.optimize(lambda model, where:self.__Callback(model, where))
        elapsed_time = time.time() - start_time

        if model.status == GRB.Status.INTERRUPTED:
            self.lb[inverter] = {"value": model.ObjBound, "type": "cut-off", "model": model}
        else:
            solution = {i[0]: [j for j in cover_range[i[0]] if model._secondary[j, i[0]].x > 1 - 1e-5]
                            for i in model._combiner if model._combiner[i].x > 1 - 1e-5}
            self.lb[inverter] = {"value": model.objVal, "type": "exact", "solution": solution}

        return model.objVal, elapsed_time
