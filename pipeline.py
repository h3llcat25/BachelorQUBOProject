from solver_manager import *


class SamplerResult:
    def __init__(self, chainstrength=None, chainstrength_multiplier=None, energy=None, enzyme_damage=None, compound_damage=None, result_graph=None, sample_result_dict=None, file_path=None):
        self.chainstrength = chainstrength if chainstrength is not None else 0
        self.chainstrength_multiplier = chainstrength_multiplier if chainstrength_multiplier is not None else 0
        self.energy = energy if energy is not None else 10000
        self.enzyme_damage = enzyme_damage if enzyme_damage is not None else 1000
        self.compound_damage = compound_damage if compound_damage is not None else 1000
        self.result_graph = result_graph
        self.sample_result_dict = sample_result_dict if sample_result_dict is not None else {}
        self.file_path = file_path

    def get_string(self):
        return f'Chainstrength/-multiplier: {self.chainstrength}/={self.chainstrength_multiplier}*(abs. Bias), Enrgy: {self.energy}, Enzyme/Compound-Damage: {self.enzyme_damage}/{self.compound_damage}'

    def set_from_modify_graph_data(self, modify_graph_return):
        enzymeDamage, compoundDamage, graph, new_filepath = modify_graph_return

        self.enzyme_damage = enzymeDamage
        self.compound_damage = compoundDamage
        self.result_graph = graph
        self.file_path = new_filepath

    def is_better_than(self, other_sampler_result):
        if other_sampler_result.energy == other_sampler_result.compound_damage == other_sampler_result.enzyme_damage == 1000:
            return True
        elif self.enzyme_damage <= other_sampler_result.enzyme_damage and self.compound_damage <= other_sampler_result.compound_damage:
            return True
        else:
            return False


def run_dwave_auto_embedding_diff_chainstrengths(file_path):
    qubo_matrix, graph, bin_vars_dict, marked_nodes = get_qubo_dict_for_dwave_qa(file_path)

    # chainstrength_multiplier = [0.25, 0.5, 0.75, 1.0, 0.25, 0.5, 0.75, 2]
    chainstrength_multiplier = [0.25, 2]

    # Wichtig für die Chainstrength ist wieviel der absolut value ist von den bias Unterschieden
    max_val = max(qubo_matrix.values())
    min_val = min(qubo_matrix.values())
    if min_val < 0:
        absolut_bias = max_val - min_val
    else:
        absolut_bias = max_val + min_val

    best_result = SamplerResult()

    # Open a file to write the results
    with (open("qubo_results.txt", "w") as file):
        for strength in chainstrength_multiplier:
            chainstrength = absolut_bias * strength
            # Run each chain strength twice
            for _ in range(2):
                # Use EmbeddingComposite for auto-embedding with specified chain strength
                sampler = DWaveSampler()
                embedding_sampler = EmbeddingComposite(sampler)
                # Sample the QUBO problem
                response = embedding_sampler.sample_qubo(qubo_matrix, chain_strength=chainstrength, num_reads=2500)

                # Find the sample with the lowest energy
                best_sample1 = min(response.data(['sample', 'energy']), key=lambda x: x.energy)
                best_sample_first = response.first
                best_energy = best_sample_first.energy

                temp_result = SamplerResult(sample_result_dict=best_sample_first, energy=best_energy, chainstrength=chainstrength, chainstrength_multiplier=strength)

                # Returns: enz_damag, cmp_damag, res_graph, end_filepath
                modify_return_data = modify_dot_graph(bin_vars_dict, graph, best_sample_first, file_path, marked_nodes)

                temp_result.set_from_modify_graph_data(modify_return_data)

                if temp_result.is_better_than(best_result):
                    best_result = temp_result

                # Write the best energy and chain strength to the file
                file.write(temp_result.get_string())
                file.write(f'best_sample1: {best_sample1}, best_sample_first: {best_sample_first}') # TODO

        # Write the best Result at the end:
        file.write(f'The best result is: {best_result.get_string()}')


def main():
    file_path = "/workspaces/graph-coloring/Citrate_cycle_marked.dot"
    run_dwave_auto_embedding_diff_chainstrengths(file_path)


if __name__ == '__main__':
    main()
