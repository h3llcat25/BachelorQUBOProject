from dimod import ConstrainedQuadraticModel, BinaryQuadraticModel, ExactSolver
from dwave.system import LeapHybridCQMSampler
from collections import defaultdict

from param_and_solve_automation import *

# ------- Main program -------
if __name__ == "__main__":
    file_path = "/workspaces/graph-coloring/Polyketidecpo.dot"
    qubo_matrix, graph, bin_vars_dict, marked_nodes = get_qubo_dict_for_dwave_qa(file_path)
    # "C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\smallDots\\eco_filtering_dot"
    # "\\Glycerolipid_marked.dot")
    # "C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\smallDots_w_marked_and_tests"
    # "\\eco_filtering_dot\\Polyketide.dot")

    bqm = BinaryQuadraticModel.from_qubo(qubo_matrix)
    solver = ExactSolver()

    # Solve the BQM using the ExactSolver
    result = solver.sample(bqm)

    best_sample = result.first
    print(best_sample)

    node_name_to_solulu_dict = node_to_solutionsValue_assignment_for_result(bin_vars_dict, best_sample.sample)

    modify_dot_graph(graph, node_name_to_solulu_dict, file_path,
                     marked_nodes)  # graph, solution_dictionary, dot_file, target_nodes
    # modify_dot_graph(graph, solution_dictionary, dot_file, target_nodes)

    # Print the results
    # print("Solutions:")
    # for sample, energy in result.data(['sample', 'energy']):
    #     print(f"Sample: {sample}, Energy: {energy}")

    # print(bqm)
    # cqm = build_cqm(G, num_colors)

    # sample = run_hybrid_solver(cqm)
