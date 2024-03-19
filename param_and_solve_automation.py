import os
import sympy
from collections import defaultdict
from gurobipy import Model, GRB, QuadExpr
from sympy import symbols

from creating_qubo_dict import *
from qMatrixBuilder import *
from act_qubo_term_generation import *


# Simplifying the Input Equation instead of putting different penalty values. This function can be called,
# instead of the replace_variables when the auto_meta_penalties (or smth) function has been chosen for the equation
# generation
def just_simplifying_objective_function(input_string, locals):
    print(input_string)
    expr = sympy.sympify(input_string, locals=locals)
    expr = sympy.simplify(expr)
    expr = sympy.expand(expr)

    print(expr)
    return expr


def node_to_binvalue_assignment_for_result(bin_var_dict, solution_array):
    # Filter the dictionary to only include keys before a key that starts with 'extra Var for'
    filtered_keys = []
    for k in bin_var_dict.keys():  # bin_var_dict ist das Array
        if k.startswith('extra Var for'):
            break
        filtered_keys.append(k)

    # Cut the length of the solution array and round the numbers to integers
    solution_array = solution_array[:len(filtered_keys)]
    solution_array = np.round(solution_array).astype(int)

    # Create a new dictionary where the keys are the filtered keys from the input
    # dictionary and the values are the first n entries of the np.ndarray
    output_dict = dict(zip(filtered_keys, solution_array))

    return output_dict


# bin_var_dict ist ein Dictionary <Node : Bin_var> wo alle Variablen aus dem Graph her mit den Knoten stehen.
# solution dict ist ein dictionary <Bin_var : solutionValue (-0, 1)>
def node_to_solutionsValue_assignment_for_result(bin_var_dict, solution_dict):
    node_to_sol = {}

    for key, value in bin_var_dict.items(): # key is node-name and value is x_i values
        if key.startswith('extra Var for'):
            continue
        solution_val = solution_dict.get(value)
        node_to_sol[key] = solution_val  # node-name to solution value

    return node_to_sol


def modify_dot_graph(graph, solution_dictionary, dot_file, target_nodes):
    # Load the .dot graph file
    rimColorAttribute = "color"
    fillColorAttribute = "fillcolor"
    styleAttribute = "style"
    filledStyle = "filled"
    rimColor = "purple"
    fillColor = "skyblue"
    targetFillColor = "red"
    enzymeFillColor = "green"

    enzymeDamage = 0
    compoundDamage = 0

    # Iterate over the dictionary and modify the graph
    for key, value in solution_dictionary.items():
        if value == 1 or value == 1.0:
            node = graph.get_node(key)
            if not node:
                node = graph.get_node(f"\"{key}\"")
            node = node[0]
            if not node:
                print("well there must be something wrong with the modify graph method or the node mmmh.")
            if key in target_nodes or (f"\"{key}\"") in target_nodes: # Deckt den Fall ab, dass
                # die Nodes nicht im Dot Graph File markiert wurden, sondern randomly oder von der disease
                # Liste des Qutie Papers
                node.set(fillColorAttribute, targetFillColor)
            else:
                node.set(rimColorAttribute, rimColor)
                node.set(fillColorAttribute, fillColor)
                node_type = node.get("type")
                if node_type == "C":
                    compoundDamage += 1
                if node_type == "E":
                    enzymeDamage += 1
                    node.set(fillColorAttribute, enzymeFillColor)
            node.set(styleAttribute, filledStyle)

        # Create a new filename with 'modified_' prefix
    base_name = os.path.basename(dot_file)
    dir_name = os.path.dirname(dot_file)
    new_filename = f"xtra_modified_{base_name}"
    new_filepath = os.path.join(dir_name, new_filename)

    print(f'Enzyme Damage: {enzymeDamage}')
    print(f'Compound Damage: {compoundDamage}')

    # Write the modified graph back to the new .dot file
    graph.write(new_filepath)


def solve_qubo_with_gurobi(input_dict, objective_expr, damage_expression=None): # (bin_vars_dict, modified_opt_term,
    # modified_damage_opt_term) bin_vars_dict hat die namen der nodes und die variablen, die anderen beiden haben einfach expressions.
    m = Model()

    # Add binary variables to the model
    variables = {}
    reverse_mapping = {}

    for key, value in input_dict.items():
        variables[value] = m.addVar(vtype=GRB.BINARY, name=value)
        reverse_mapping[value] = key
    m.update()

    # Convert the sympy expression to a dictionary
    objective_dict_whole = objective_expr.as_coefficients_dict()
    if damage_expression is not None:
        objective_dict_damage = damage_expression.as_coefficients_dict()

        # Merge the two dictionaries
        for term, coeff in objective_dict_damage.items():
            if term in objective_dict_whole:
                objective_dict_whole[term] += coeff
            else:
                objective_dict_whole[term] = coeff

    # Set the objective
    objective = QuadExpr()

    for term, coeff in objective_dict_whole.items():
        if isinstance(term, sympy.Symbol):  # linear term
            objective.add(variables[str(term)], coeff)
        elif len(term.args) == 2:  # quadratic term
            var1, var2 = map(str, term.args)
            if var1 == var2:
                print(var1)
            if var2 == "2":
                objective.add(variables[var1], coeff)
            else:
                objective.add(variables[var1] * variables[var2], coeff)
        else:  # constant term
            objective.addConstant(coeff)
    m.setObjective(objective, GRB.MINIMIZE)

    # Solve the model
    m.optimize()

    # Print the solution
    solution = {}
    for var_name, var in variables.items():
        if not reverse_mapping[var_name].startswith("extra Var"):
            solution[reverse_mapping[var_name]] = var.X
    return solution


# TODO: Aktuell.
def param_search_and_auto_solver_auto_penalty_gurobipy_model_enzymes_and_matrix(file_path):
    bin_vars_dict, optimization_term, output_damage_only, graph, marked_nodes = generating_qubo_term_from_graph_two_part(file_path)
    variables = {symbol: symbols(symbol) for symbol in set(bin_vars_dict.values())}

    modified_opt_term = just_simplifying_objective_function(optimization_term, variables)
    modified_damage_opt_term = sympy.sympify(output_damage_only, locals=variables)

    print(modified_damage_opt_term)
    print(bin_vars_dict)

    qmatrix = generate_final_np_matrix(bin_vars_dict, str(modified_opt_term))
    qmatrix_damage = generate_final_np_matrix(bin_vars_dict, str(modified_damage_opt_term))

    cmplt_matrix = add_damage_matrix_to_q_matrix(qmatrix, qmatrix_damage)

    qubo_dict = defaultdict(int)
    n = cmplt_matrix.shape[0]  # Assuming a square matrix
    bin_vars_list = list(bin_vars_dict.values())

    for i in range(n):
        for j in range(i, n):  # Only need to iterate over the upper triangle due to symmetry
            if cmplt_matrix[i][j] != 0:  # Only consider non-zero entries
                qubo_dict[(bin_vars_list[i], bin_vars_list[j])] = cmplt_matrix[i][j]

    solution_dict = solve_qubo_with_gurobi(bin_vars_dict, modified_opt_term, modified_damage_opt_term)


    # nodeNames_to_binaryMap = node_to_binvalue_assignment_for_result(bin_vars_dict, result.solution)

    # Print the solution
    # print(nodeNames_to_binaryMap)
    # print(result.solution)
    # print(result.objective_value)
    modify_dot_graph(graph, solution_dict, file_path, marked_nodes)  # graph, solution_dictionary, dot_file, target_nodes
    return qubo_dict

    # print(solution_dict)


def old_get_qubo_matrix_for_dwave_qa(file_path):
    bin_vars_dict, optimization_term, output_damage_only, graph, marked_nodes = generating_qubo_term_from_graph_two_part(file_path)
    variables = {symbol: symbols(symbol) for symbol in set(bin_vars_dict.values())}

    modified_opt_term = just_simplifying_objective_function(optimization_term, variables)
    modified_damage_opt_term = sympy.sympify(output_damage_only, locals=variables)

    print(modified_damage_opt_term)
    print(bin_vars_dict)

    qmatrix = generate_final_np_matrix(bin_vars_dict, str(modified_opt_term))
    qmatrix_damage = generate_final_np_matrix(bin_vars_dict, str(modified_damage_opt_term))

    cmplt_matrix = add_damage_matrix_to_q_matrix(qmatrix, qmatrix_damage)

    qubo_dict = defaultdict(int)
    n = cmplt_matrix.shape[0]  # Assuming a square matrix
    bin_vars_list = list(bin_vars_dict.values())

    for i in range(n):
        for j in range(i, n):  # Only need to iterate over the upper triangle due to symmetry
            if cmplt_matrix[i][j] != 0:  # Only consider non-zero entries
                qubo_dict[(bin_vars_list[i], bin_vars_list[j])] = cmplt_matrix[i][j]

    return qubo_dict, graph, bin_vars_dict, marked_nodes

    # print(solution_dict)


def get_qubo_dict_for_dwave_qa(file_path):
    bin_vars_dict, optimization_term, output_damage_only, graph, marked_nodes = generating_qubo_term_from_graph_two_part(file_path)
    variables = {symbol: symbols(symbol) for symbol in set(bin_vars_dict.values())}

    modified_opt_term = just_simplifying_objective_function(optimization_term, variables)
    modified_damage_opt_term = sympy.sympify(output_damage_only, locals=variables)

    q_dict = creating_qubo_dict_sympy(modified_opt_term, modified_damage_opt_term)

    return q_dict, graph, bin_vars_dict, marked_nodes


# MAIN
def main():
    # param_search_and_auto_solver_auto_penalty(
    # param_search_and_auto_solver_auto_penalty_two_parted(
    tests_for_getting_qubo_matrix(
    # param_search_and_auto_solver_auto_penalty_gurobipy_model_enzymes_and_matrix(
        # "C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\smallDots\\eco_filtering_dot"
        # "\\Glycerolipid_marked.dot")
        "C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\smallDots_w_marked_and_tests"
        "\\eco_filtering_dot\\Polyketide.dot")
        # "C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\graphoz\\a.dot")
        #  "C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\graphoz\\a.dot")
    # "\\Biosynthesis of amino acids marked.dot")


if __name__ == '__main__':
    main()

# embedding = self.graph_embedder.get_embedding_for_bundle(Q,bundle)
#             sampler = FixedEmbeddingComposite(DWaveSampler(), embedding)
#             sampleset = sampler.sample_qubo(Q, num_reads=self.annealing_runs)
#             response = sampleset.first.sample