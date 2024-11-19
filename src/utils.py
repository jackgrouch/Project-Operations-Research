import numpy as np
import gurobipy as gp
from gurobipy import GRB

import data as dt
import subproblem_solver as ss
import solution as sol


def master_problem(data: dt.Data, sp_list: list[ss.DistrictSolver]):
    lb_way = [{w: None for w in sp_list[k].way} for k in range(data.district_num)]
    big_m = [0] * data.district_num

    for k in range(data.district_num):
        for r in sp_list[k].way:
            candidates = np.intersect1d(sp_list[k].nodes, data.way_set[r])
            lb_way[k][r] = min(sp_list[k].lb[i]['value'] for i in candidates)
            big_m[k] = max(big_m[k], lb_way[k][r])

    model = gp.Model()
    y = model.addVars(data.way_num, vtype=GRB.BINARY)
    z = {
        k: model.addVars(sp_list[k].way, vtype=GRB.BINARY)
        for k in range(data.district_num)
    }

    for k in range(data.district_num):
        model.addConstr(gp.quicksum(y[r] for r in sp_list[k].way) >= 1)
        model.addConstrs(z[k][r] <= y[r] for r in sp_list[k].way)
        model.addConstr(gp.quicksum(z[k][r] for r in sp_list[k].way) == 1)

    mp_obj = data.way_cost * gp.quicksum(data.way_len[i] * y[i] for i in range(data.way_num))
    sp_obj = gp.quicksum(
        lb_way[k][r] * z[k][r]
        for k in range(data.district_num)
        for r in sp_list[k].way
    )

    model.setObjective(mp_obj + sp_obj, GRB.MINIMIZE)
    model.Params.OutputFlag = 0
    model.optimize()

    if model.status == GRB.Status.OPTIMAL:
        result = sol.ResultMaster(z, model, mp_obj, data, sp_list)

    return result, model.Runtime


def worker_problem(data: dt.Data, sp_list: list[ss.DistrictSolver], iteration: int, way_select, solution: sol.Solution):
    cut_off = False
    result = {}
    total_time = 0
    order = range(data.district_num)
    dp_dict = {k: {} for k in order}

    for k in order:
        candidates = np.intersect1d(data.way_set[way_select[k]], sp_list[k].nodes)
        dp_dict[k] = {i: sp_list[k].lb[i] for i in candidates}
        dp_dict[k] = sorted(dp_dict[k].items(), key=lambda x: x[1]['value'])

    if iteration > 1:
        order = sorted(
            order,
            key=lambda x: dp_dict[x][0][1]['value'] - solution.district_solutions[x].subproblem_obj_value,
            reverse=True,
        )

    for k in order:
        gp_time = 0
        lb = 0
        ub = float('inf')

        for inverter, state in dp_dict[k]:
            lb = max(lb, state['value'])
            if ub <= lb:
                break

            runtime = 0
            if state['type'] != 'exact':
                print(f"\t \t District[{k + 1}] : solve [{inverter}],", end=' ')
                obj, runtime = sp_list[k].exact_solve(inverter, ub)
                print(f"{sp_list[k].lb[inverter]['type']} obj = {obj:.2f}, time = {runtime:.1f}s")

            if ub > sp_list[k].lb[inverter]['value']:
                ub = sp_list[k].lb[inverter]['value']
                result[k] = sol.ResultWorker(k, ub, inverter)
            gp_time += runtime

        print(f"\t \t District[{k + 1}] : optimal value = {ub:.2f}, time = {gp_time:.1f}s")
        total_time += gp_time

        if iteration > 1 and gp_time > 0:
            master, mp_time = master_problem(data, sp_list)
            gap = (solution.objective_value - master.total_obj_value) / solution.objective_value
            print(f"\t \t mp: obj = {master.total_obj_value:.2f}, time = {mp_time:.1f}s, gap = {gap * 100:.2f}%")
            total_time += mp_time
            if gap <= 1e-5:
                print(f"\t sub-problem cut off, time = {total_time:.1f}s")
                cut_off = True
                break

    return result, total_time, cut_off
