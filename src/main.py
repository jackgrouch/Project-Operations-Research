import sys
import os

from data import Data
from solution import Solution
from subproblem_solver import DistrictSolver
import utils

if __name__ == "__main__":
    file = sys.argv[1]
    file_name = os.path.basename(file)
    print(f"Instance: {file_name.split('.')[0]}")

    data = Data(file)
    sp_list = [DistrictSolver(data, i) for i in range(data.district_num)]
    print(f'Solving LP relaxation...')
    total_time = sum(sp.relax_all() for sp in sp_list)
    print(f'Done, solving time = {total_time:.1f}s')

    lower_bound = 0
    upper_bound = float('inf')
    iteration = 1
    solution = Solution(float('inf'), None, None)

    while True:
        print(f'Iteration {iteration}:')
        master, mp_time = utils.master_problem(data, sp_list)
        lower_bound = max(lower_bound, master.total_obj_value)
        gap = (upper_bound - lower_bound) / upper_bound

        print(f'\t master problem : obj = {master.total_obj_value:.2f}, service ways = {set(master.selected_ways.values())}, time = {mp_time:.1f}s, gap = {gap * 100:2.2f}% ')
        workers, sp_time, cut_off = utils.worker_problem(data, sp_list, iteration, master.selected_ways, solution)

        if cut_off:
            master, tmp = utils.master_problem(data, sp_list)
            lower_bound = max(lower_bound, master.total_obj_value)
        else:
            sp_total_obj = sum(worker.subproblem_obj_value for worker in workers.values())
            print(f'\t sub-problem : obj = {sp_total_obj}, time = {sp_time:.1f}s')
            if upper_bound > master.master_obj_value + sp_total_obj:
                upper_bound = master.master_obj_value + sp_total_obj
                solution = Solution(upper_bound, master.selected_ways, workers)


        total_time += mp_time + sp_time
        gap = (upper_bound - lower_bound) / upper_bound
        print(f'\t ub = {upper_bound:.2f}, lb = {lower_bound:.2f}, gap = {gap * 100:2.2f}% ')

        if gap <= 1e-5:
            print(f'Done! The objective value is {upper_bound:.2f}, total solving time = {total_time:.1f}s')
            break

        iteration += 1

    if len(sys.argv) > 2 :
        solution.output_to_file(sys.argv[2], sp_list)
