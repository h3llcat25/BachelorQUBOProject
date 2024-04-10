from solver_manager import *


class SamplerResult:
    def __init__(self, chainstrength=None, chainstrength_multiplier=None, energy=None, enzyme_damage=None,
                 compound_damage=None, result_graph=None, sample_result_dict=None, file_path=None):
        if chainstrength:
            self.chainstrength = int(round(chainstrength))
        else:
            self.chainstrength = 0.0
        self.chainstrength_multiplier = chainstrength_multiplier if chainstrength_multiplier is not None else 0
        self.energy = energy if energy is not None else 10000
        self.enzyme_damage = enzyme_damage if enzyme_damage is not None else 1000
        self.compound_damage = compound_damage if compound_damage is not None else 1000
        self.result_graph = result_graph
        self.sample_result_dict = sample_result_dict if sample_result_dict is not None else {}
        self.file_path = file_path

    def get_string(self):
        return f'Chainstrength/-multiplier: {self.chainstrength}/={self.chainstrength_multiplier}*(abs. Bias), Enrgy: {self.energy}, Enzyme/Compound-Damage: {self.enzyme_damage}/{self.compound_damage} '

    def set_from_modify_graph_data(self, modify_graph_return):
        enzymeDamage, compoundDamage, graph, new_filepath = modify_graph_return

        self.enzyme_damage = enzymeDamage
        self.compound_damage = compoundDamage
        self.result_graph = graph
        self.file_path = new_filepath

    def is_better_than(self, other_sampler_result):
        if other_sampler_result.energy == other_sampler_result.compound_damage == other_sampler_result.enzyme_damage == 1000:
            return True
        elif self.enzyme_damage <= other_sampler_result.enzyme_damage:
            return True
        else:
            return False


def run_dwave_sampler_or_solver(sampler_or_solver, file_path):
    qubo_matrix, graph, bin_vars_dict, marked_nodes = get_qubo_dict_for_dwave_qa(file_path)

    results_list = []

    chainstrength_multiplier = [5.0, 10.0, 30.0]
    # chainstrength_multiplier = [0.25, 2]

    # Wichtig für die Chainstrength ist wieviel der absolut value ist von den bias Unterschieden
    '''
    max_val = max(qubo_matrix.values())
    min_val = min(qubo_matrix.values())
    if min_val < 0:
        absolut_bias = max_val - min_val
    else:
        absolut_bias = max_val + min_val
    '''
    # Open a file to write the results
    with (open("qubo_results.txt", "w") as file):
        for strength in chainstrength_multiplier:
            chainstrength = strength
            # Run each chain strength twice
            for _ in range(2):

                # Part to decide which Sampler/Solver
                if sampler_or_solver == 'hybrid':
                    sampler = LeapHybridSampler()
                elif sampler_or_solver == 'simulated_annealing':
                    sampler = SimulatedAnnealingSampler()
                elif sampler_or_solver == "dwave_auto_emb":
                    sampler = EmbeddingComposite(DWaveSampler())
                else:
                    raise ValueError("Invalid solver type specified.")

                if sampler_or_solver in ['hybrid', 'simulated_annealing']:
                    result = sampler.sample_qubo(qubo_matrix)
                else:
                    result = sampler.sample_qubo(qubo_matrix, chain_strength=int(chainstrength),
                                                 num_reads=2500)  # Specify num_reads for quantum samplers

                # Find the sample with the lowest energy
                best_sample_first = result.first.sample
                best_energy = result.first.energy

                temp_result = SamplerResult(sample_result_dict=best_sample_first, energy=best_energy,
                                            chainstrength=chainstrength, chainstrength_multiplier=strength)

                # Returns: enz_damag, cmp_damag, res_graph, end_filepath
                modify_return_data = modify_dot_graph(bin_vars_dict, best_sample_first, file_path, marked_nodes)

                temp_result.set_from_modify_graph_data(modify_return_data)

                # Write the best energy and chain strength to the file
                file.write(temp_result.get_string())
                file.write("\n")

                if results_list:
                    if temp_result.is_better_than(results_list[0]):
                        results_list.insert(0, temp_result)
                        continue
                results_list.append(temp_result)

        # Write the best Result at the end:
        file.write(f'The best total result is: {results_list[0].get_string()}')
        file.write("\n\n")

        # Write the qubo matrix dictionary inside of the file
        i = 0
        for key, value in qubo_matrix.items():
            file.write(f'{key}: {value},  ')
            i += 1
            if i == 9:
                file.write("\n")
                i = 0

    directory = "best_graphs"
    if not os.path.exists(directory):
        os.makedirs(directory)

    # Step 2: Process the first 5 graph objects
    for obj in results_list[:5]:
        # Generate the filename using a, b, and c values
        filename = f"{obj.enzyme_damage}(){obj.compound_damage}(){obj.energy}.dot"
        filepath = os.path.join(directory, filename)

        obj.result_graph.del_node('\"\\r\\n\"')

        obj.result_graph.write(filepath)


def main():
    file_path = "/workspaces/graph-coloring/Citrate_cycle_marked.dot"
    sampler_or_solver = "dwave_auto_emb"
    run_dwave_sampler_or_solver(sampler_or_solver, file_path)


if __name__ == '__main__':
    main()
