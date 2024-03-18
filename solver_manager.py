from dimod import ConstrainedQuadraticModel, BinaryQuadraticModel, ExactSolver
from dwave.system import LeapHybridCQMSampler
from collections import defaultdict

from param_and_solve_automation import *

# ------- Main program -------
if __name__ == "__main__":

    qubo_matrix = param_search_and_auto_solver_auto_penalty_gurobipy_model_enzymes_and_matrix(
        # "C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\smallDots\\eco_filtering_dot"
        # "\\Glycerolipid_marked.dot")
        "C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\smallDots_w_marked_and_tests"
        "\\eco_filtering_dot\\Polyketide.dot")

    bqm = BinaryQuadraticModel.from_qubo(qubo_matrix)
    solver = ExactSolver()

    # Solve the BQM using the ExactSolver
    result = solver.sample(bqm)

    # Print the results
    print("Solutions:")
    for sample, energy in result.data(['sample', 'energy']):
        print(f"Sample: {sample}, Energy: {energy}")
    #print(bqm)
    #cqm = build_cqm(G, num_colors)

    #sample = run_hybrid_solver(cqm)
