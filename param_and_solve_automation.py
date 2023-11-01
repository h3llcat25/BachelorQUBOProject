import sympy
import numpy as np
import os
import pydot
from gurobi_optimods.qubo import solve_qubo
import qubo_term_from_graph_generator as qtfgg
import qMatrixBuilder as qmtrxb

#from dwave.samplers import SimulatedAnnealingSampler
#from dwave.system import DWaveSampler, FixedEmbeddingComposite, LeapHybridSampler
#import dimod

# from dwave_qbsolv import QBSolv


# To replace k1 till k6 of the objective functions with different weights
def replace_variables(input_string, k1, k2, k3, k4, k5, k6):
    input_string = input_string.replace('k1', str(k1))
    input_string = input_string.replace('k2', str(k2))
    input_string = input_string.replace('k3', str(k3))
    input_string = input_string.replace('k4', str(k4))
    input_string = input_string.replace('k5', str(k5))
    input_string = input_string.replace('k6', str(k6))

    expr = sympy.sympify(input_string)
    expr = sympy.simplify(expr)
    expr = sympy.expand(expr)

    return expr


# Simplifying the Input Equation instead of putting different penalty values. This function can be called,
# instead of the replace_variables when the auto_meta_penalties (or smth) function has been chosen for the equation
# generation
def just_simplifying_objective_function(input_string):
    print(input_string)
    expr = sympy.sympify(input_string)
    expr = sympy.simplify(expr)
    expr = sympy.expand(expr)

    return expr


def node_to_binvalue_assignment_for_result(bin_var_dict, solution_array):
    # Filter the dictionary to only include keys before a key that starts with 'extra Var for'
    filtered_keys = []
    for k in bin_var_dict.keys():
        if str(k).startswith('extra Var for'):
            break
        filtered_keys.append(k)

    # Cut the length of the solution array and round the numbers to integers
    solution_array = solution_array[:len(filtered_keys)]
    solution_array = np.round(solution_array).astype(int)

    # Create a new dictionary where the keys are the filtered keys from the input
    # dictionary and the values are the first n entries of the np.ndarray
    output_dict = dict(zip(filtered_keys, solution_array))

    return output_dict


def modify_dot_graph(dot_file, dictionary):

    # Load the .dot graph file
    graph = pydot.graph_from_dot_file(dot_file)[0]
    attribute="color"
    color="purple"

    # Iterate over the dictionary and modify the graph
    for key, value in dictionary.items():
        if value == 1:
            try:
                node = graph.get_node(key)[0]
                if len(node.get_attributes()) == 1 and 'type' in node.get_attributes():
                    node.set(attribute, color)
            except IndexError:
                node = graph.get_node(f'"{key}"')[0]
                if len(node.get_attributes()) == 1 and 'type' in node.get_attributes():
                    node.set(attribute, color)

        # Create a new filename with 'modified_' prefix
    base_name = os.path.basename(dot_file)
    dir_name = os.path.dirname(dot_file)
    new_filename = f"modified_{base_name}"
    new_filepath = os.path.join(dir_name, new_filename)

    # Write the modified graph back to the new .dot file
    graph.write(new_filepath)


def param_search_and_auto_solver_auto_penalty(file_path):
    bin_vars_dict, optimization_term = qtfgg.generating_qubo_term_from_graph(file_path)
    modified_opt_term = just_simplifying_objective_function(optimization_term)

    print(modified_opt_term)

    qmatrix = qmtrxb.generate_final_np_matrix(bin_vars_dict, str(modified_opt_term))
    result = solve_qubo(qmatrix)

    nodeNames_to_binaryMap = node_to_binvalue_assignment_for_result(bin_vars_dict, result.solution)

    # Print the solution
    print(nodeNames_to_binaryMap)
    print(result.solution)
    print(result.objective_value)
    modify_dot_graph(file_path, nodeNames_to_binaryMap)


def param_search_and_auto_solver_with_manual_penalty_setting(is_gurobi, file_path, k1, k2, k3, k4, k5=1, k6=1):
    bin_vars_dict, optimization_term = qtfgg.generating_qubo_term_from_graph(file_path)
    modified_opt_term = replace_variables(optimization_term, k1, k2, k3, k4, k5, k6)

    print(modified_opt_term)

    qmatrix = qmtrxb.generate_final_np_matrix(bin_vars_dict, str(modified_opt_term))

    if is_gurobi:
        result = solve_qubo(qmatrix)

        nodeNames_to_binaryMap = node_to_binvalue_assignment_for_result(bin_vars_dict, result.solution)

        # Print the solution
        print(nodeNames_to_binaryMap)
        # print(result.solution)
        print(result.objective_value)
    else:
        # Populate the matrix with values from the dictionary
        matrix_dict = {index: value for index, value in np.ndenumerate(qmatrix) if value != 0}

        # Use the QBSolv module to solve the QUBO problem
        #response = QBSolv().sample_qubo(matrix_dict)

        # Print the samples and energies
        #print("samples=" + str(list(response.samples())))
        #print("energies=" + str(list(response.data_vectors['energy'])))

# MAIN
def main():
    param_search_and_auto_solver_auto_penalty(
        "C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\smallDots\\eco_filtering_dot"
        "\\Citrate_cycle_marked.dot")

        #"C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\largeDots\\hsa_filtering_dot"
        #"\\Nucleotide metabolism_marked.dot")


if __name__ == '__main__':
    main()
