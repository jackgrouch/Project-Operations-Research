import data as dt
import subproblem_solver as ss


class ResultMaster:
    def __init__(self, decision_vars, model, master_problem_obj, data: dt.Data, subproblem_solvers: list[ss.DistrictSolver]):
        self.selected_ways = {
            district: way
            for district in range(data.district_num)
            for way in subproblem_solvers[district].way
            if decision_vars[district][way].x >= 1 - 1e-5
        }
        self.master_obj_value = master_problem_obj.getValue()
        self.total_obj_value = model.objVal


class ResultWorker:
    def __init__(self, district_index, subproblem_obj, inverter_index):
        self.district = district_index + 1
        self.subproblem_obj_value = subproblem_obj
        self.inverter = inverter_index


class Solution:
    def __init__(self, objective_value, selected_ways, district_solutions):
        self.objective_value = objective_value
        self.selected_ways = selected_ways
        self.district_solutions = district_solutions

    def output_to_file(self, file_path, subproblem_solvers: list[ss.DistrictSolver]):
        with open(file_path, 'w') as file:
            print(f"The optimal value is {self.objective_value:.4f}", file=file)
            print(f"Service ways: {set(self.selected_ways.values())}", file=file)
            print("Inverters:", end=" ", file=file)
            for district in self.district_solutions.values():
                print(district.inverter, end=" ", file=file)
            print("\nCable routing:", file=file)
            for district_index, district_solution in self.district_solutions.items():
                solution = subproblem_solvers[district_index].lb[district_solution.inverter]["solution"]
                print(solution, file=file)
