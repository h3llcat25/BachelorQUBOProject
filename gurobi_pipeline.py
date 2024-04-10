import itertools

import sympy
from sympy import symbols

from act_qubo_term_generation import read_dot_file_pydotplus, \
    create_binary_variables_plus_enzymes, generate_equation
from param_and_solve_automation import just_simplifying_objective_function, solve_qubo_with_gurobi, modify_dot_graph


class GurobiResult:
    def __init__(self, compound_d_weight, enzyme_d_weight, result_graph=None, solution_result_dict=None, file_path=None):
        self.enzyme_d_weight = enzyme_d_weight
        self.compound_d_weight = compound_d_weight
        self.result_graph = result_graph
        self.solution_result_dict = solution_result_dict if solution_result_dict is not None else {}
        self.file_path = file_path

        self.enzyme_damage=0
        self.compound_damage=0

    def get_string(self):
        return f'Enzyme/Compound-Weights: {self.enzyme_d_weight}/{self.compound_d_weight}: Final Enzyme/Compound-Damage: {self.enzyme_damage}/{self.compound_damage} \n'

    def set_from_modify_graph_data(self, modify_graph_return):
        enzymeDamage, compoundDamage, graph, new_filepath = modify_graph_return

        self.enzyme_damage = enzymeDamage
        self.compound_damage = compoundDamage
        self.result_graph = graph
        self.file_path = new_filepath

    def get_enz_comp_dmg(self):
        return self.enzyme_damage, self.compound_damage


def run_gurobi_on_different_biases(file_path, input_arguments):
    #combinations = list(itertools.product(input_arguments, repeat=2))
    combinations = list(itertools.combinations(input_arguments, 2))
    results_list = []

    graph = read_dot_file_pydotplus(file_path)
    if not graph:
        print("The Graph is weirdly, None...")
        return None

    # Open a file to write the results
    with (open("gurobi_results.txt", "w") as file):
        # Iterate over the combinations and apply the function
        for combo in combinations:
            for _ in range(3):
                dict_and_sorted_nodes = create_binary_variables_plus_enzymes(graph)  # TODO WICHTIG, ich habe hier  k=3
                # gesetzt, weil der Polyketide Graph keine diseased Knoten hat!!
                if not dict_and_sorted_nodes:
                    print("there is somehow not a dict_and_sorted_nodes")
                    return
                if dict_and_sorted_nodes.has_seed():
                    seed = dict_and_sorted_nodes.get_seed()

                enzy, compo = combo
                tie_qubo_struct = generate_equation(graph, dict_and_sorted_nodes, compo, enzy)
                bin_vars_dict, optimization_term, output_damage_only, graph, marked_nodes = tie_qubo_struct.get_dict_objectives_graph_targets()
                variables = {symbol: symbols(symbol) for symbol in set(bin_vars_dict.values())}

                modified_opt_term = just_simplifying_objective_function(optimization_term, variables)
                modified_damage_opt_term = sympy.sympify(output_damage_only, locals=variables)

                solution_dict = solve_qubo_with_gurobi(bin_vars_dict, modified_opt_term, modified_damage_opt_term)

                # The first value of combo is the compound, and the second the enzyme weight value
                temp_result = GurobiResult(compo, enzy, solution_result_dict=solution_dict)

                # Returns: enz_damag, cmp_damag, res_graph, end_filepath
                modify_return_data = modify_dot_graph(bin_vars_dict, solution_dict, file_path, marked_nodes)
                temp_result.set_from_modify_graph_data(modify_return_data)

                # Write the best energy and chain strength to the file
                file.write(temp_result.get_string())

                if temp_result.get_enz_comp_dmg() != (3,5):
                    results_list.append(temp_result)
                else:
                    results_list.insert(0, temp_result)


        '''directory = "best_graphs"
        if not os.path.exists(directory):
            os.makedirs(directory)

        # Step 2: Process the first 5 graph objects
        for obj in results_list[:5]:
            # Generate the filename using a, b, and c values
            filename = f"{obj.enzyme_damage}(){obj.compound_damage}(){obj.energy}.dot"
            filepath = os.path.join(directory, filename)

            obj.result_graph.del_node('\"\\r\\n\"')

            obj.result_graph.write(filepath)'''


def main():
    file_path = ("C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\largeDots_w_marked_and_tests\\eco_filtering_dot\\Biosynthesis of amino acids_m.dot")

        # ("C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\smallDots_w_marked_and_tests"
        #          "\\eco_filtering_dot\\Citrate_cycle_marked.dot")
    input_arguments = [2, 5]
    run_gurobi_on_different_biases(file_path, input_arguments)


if __name__ == '__main__':
    main()
