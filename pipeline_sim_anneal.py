from solver_manager import *
import sympy
from sympy import symbols
import time
import os

from param_and_solve_automation import just_simplifying_objective_function, solve_qubo_with_gurobi, modify_dot_graph, \
    get_qubo_dict_for_dwave_qa


class SamplerResult:
    def __init__(self, energy=None, enzyme_damage=None,
                 compound_damage=None, result_graph=None, sample_result_dict=None, file_path=None):
        self.energy = energy if energy is not None else 10000
        self.compound_damage = compound_damage if compound_damage is not None else 1000
        self.result_graph = result_graph
        self.sample_result_dict = sample_result_dict if sample_result_dict is not None else {}
        self.file_path = file_path

    def get_comp_damage(self):
        return self.compound_damage

    def set_from_modify_graph_data(self, modify_graph_return):
        enzymeDamage, compoundDamage, graph, new_filepath = modify_graph_return

        self.enzyme_damage = enzymeDamage
        self.compound_damage = compoundDamage
        self.result_graph = graph
        self.file_path = new_filepath


def run_sim_ann_different_random_targets(file_path):
    # combinations = list(itertools.product(input_arguments, repeat=2))
    combinations = [2,7,12,17,22,27]
    #combinations = [22,27]
    results_list = []

    # Open a file to write the results
    with (open(f"simm_ann_for_diff_rndm_trgts_eco_small.txt", "a") as file):
        file.write(file_path)
        file.write("\n")
        # Iterate over the combinations and apply the function
        for target_nr in combinations:
            for _ in range(5):
                temp_var = get_qubo_dict_for_dwave_qa(file_path, target_nr=target_nr)
                if temp_var == -1:
                    continue
                qubo_matrix, graph, bin_vars_dict, marked_nodes = get_qubo_dict_for_dwave_qa(file_path, target_nr=target_nr)
                start_time = time.time()
                sampler = SimulatedAnnealingSampler()
                result = sampler.sample_qubo(qubo_matrix, num_reads=2500)

                end_time = time.time()
                # Calculate the total time it took to run the loop
                elapsed_time = end_time - start_time

                best_sample_first = result.first.sample
                best_energy = result.first.energy

                temp_result = SamplerResult(sample_result_dict=best_sample_first, energy=best_energy)

                # Returns: enz_damag, cmp_damag, res_graph, end_filepath
                modify_return_data = modify_dot_graph(bin_vars_dict, best_sample_first, file_path, marked_nodes)

                temp_result.set_from_modify_graph_data(modify_return_data)

                # Write the best energy and chain strength to the file
                file.write(f'{target_nr} , {temp_result.get_comp_damage()} , {len(bin_vars_dict)}, {elapsed_time}\n')





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
    #directory = "graphStuff\\largeDots\\eco_filterin_dot"
    directory = "graphStuff\\smallDots\\eco_filtering_dot"

    abs_directory = os.path.abspath(directory)

    # List all entries in the directory
    with os.scandir(abs_directory) as entries:
        for entry in entries:
            if entry.is_file():  # Check if the entry is a file
                # Get the full path of the file
                full_path = entry.path
                # Split the path to manipulate and extract the desired part
                path_parts = full_path.split(os.sep)
                # Include up to two directory levels above
                relevant_path = os.sep.join(path_parts[max(0, len(path_parts) - 4):])
                run_sim_ann_different_random_targets(relevant_path)


    #file_path = ("C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\smallDots_w_marked_and_tests\\eco_filtering_dot\\Citrate_cycle.dot")

        # ("C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\smallDots_w_marked_and_tests"
        #          "\\eco_filtering_dot\\Citrate_cycle_marked.dot")
    #input_arguments = [2, 5]
    #run_gurobi_on_different_biases(file_path, input_arguments)
    #run_sim_ann_different_random_targets(file_path)

if __name__ == '__main__':
    #run_sim_ann_different_random_targets("C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\smallDots\\mmu_filtering_dot\\Pyruvate.dot")
    main()
