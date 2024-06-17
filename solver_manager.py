from dimod import BinaryQuadraticModel, ExactSolver
from dwave.system import LeapHybridSampler, DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
from dwave.samplers import SimulatedAnnealingSampler
from minorminer import find_embedding

import dwave.inspector

from param_and_solve_automation import *


def solve_w_exact_solver(bqm):
    solver = ExactSolver()
    return solver.sample(bqm).first.sample


def solve_hybrid(qubo):
    # Initialize the hybrid solver
    sampler = LeapHybridSampler()

    # Submit the QUBO problem to the hybrid solver
    sampleset = sampler.sample_qubo(qubo, time_limit=20)

    # dwave.inspector.show(sampleset)
    solution = sampleset.first.sample
    return solution


def solve_w_simulated_annealing_sampler(qubo):
    sampler = SimulatedAnnealingSampler()

    # Solve the QUBO
    sample_set = sampler.sample_qubo(qubo, num_reads=2500)
    return sample_set.first.sample


def solve_w_dwave_sampler_fixed_empedding(q_matrix):
    # Anzahl der Runs
    num_runs = 500  # Kann angepasst werden, um optimale Ergebnisse zu finden

    # Anzahl der Runs
    max_val = max(q_matrix.values())
    min_val = min(q_matrix.values())
    if min_val < 0:
        chainstrength = max_val - min_val
        # chainstrength = 1

    # Initialize the D-Wave sampler
    sampler = DWaveSampler()
    interactions = set(q_matrix.keys())
    embedding = find_embedding(interactions, sampler.edgelist)
    fixed_sampler = FixedEmbeddingComposite(sampler, embedding)

    response = fixed_sampler.sample_qubo(q_matrix, chain_strength=float(chainstrength) * 1.5, num_reads=num_runs)
    dwave.inspector.show(response)
    # best_sample = min(response.data(['sample', 'energy']), key=lambda x: x.energy)
    return response.first.sample


def solve_w_dwave_sampler_auto_empedding(Q):
    # Using AutoEmbeddingComposite.sample(bqm, **parameters):
    # - Flexibility: Allows for advanced preprocessing and manipulation of both Ising and QUBO models.
    # - Unified Interface: Provides a consistent approach for various problem types and samplers.
    # Using     AutoEmbeddingComposite.sample_qubo(Q, ...):
    # Convenience: Directly accepts QUBO dictionaries, ideal for straightforward QUBO problems.
    # Simplicity: Offers an easy - to - use method with potentially less overhead for QUBO-specific workflows.
    sampler = EmbeddingComposite(DWaveSampler())

    # Anzahl der Runs
    max_val = max(Q.values())
    min_val = min(Q.values())
    if min_val < 0:
        chainstrength = max_val - min_val
        # chainstrength = 1
    num_runs = 2500  # Kann angepasst werden, um optimale Ergebnisse zu finden

    # Simulated Annealing durchführen
    bqm = BinaryQuadraticModel.from_qubo(Q)
    # response = sampler.sample(bqm, chain_strength=float(chainstrength), num_reads=num_runs)
    response = sampler.sample_qubo(Q, chain_strength=float(chainstrength), num_reads=num_runs)

    # Ergebnisse ausgeben
    for datum in response.data(['sample', 'energy', 'num_occurrences']):
        print(datum.sample, "Energie:", datum.energy, "Anzahl Vorkommen:", datum.num_occurrences)

    dwave.inspector.show(response)

    # Das optimale Ergebnis ermitteln
    best_sample = response.first.sample
    best_energy = response.first.energy
    print(best_energy)

    return best_sample


def switch_case_if(value, bqm, q_matrix):
    if value == 'exact_bqm':
        return solve_w_exact_solver(bqm)
    elif value == 'simulated_anneal_qubo':
        return solve_w_simulated_annealing_sampler(q_matrix)
    elif value == 'dwave_quantum':
        return solve_w_dwave_sampler_auto_empedding(q_matrix)
    elif value == 'dwave_sampler_fixed_empedding':
        return solve_w_dwave_sampler_fixed_empedding(q_matrix)
    elif value == "hybrid":
        return solve_hybrid(q_matrix)
    elif value == 'c':
        return "Third case"
    else:
        return "Default case"


def auto_emb_solv_chainstrength(file_path, ch_str_vals):
    qubo_matrix, graph, bin_vars_dict, marked_nodes = get_qubo_dict_for_dwave_qa(file_path)
    sampler = EmbeddingComposite(DWaveSampler())

    num_runs = 2500  # Kann angepasst werden, um optimale Ergebnisse zu finden

    # Simulated Annealing durchführen
    bqm = BinaryQuadraticModel.from_qubo(qubo_matrix)
    # response = sampler.sample(bqm, chain_strength=float(chainstrength), num_reads=num_runs)
    response = sampler.sample_qubo(qubo_matrix, chain_strength=10.0, num_reads=num_runs)

    dwave.inspector.show(response)

    # Das optimale Ergebnis ermitteln
    best_sample = response.first.sample
    best_energy = response.first.energy
    print(best_energy)

    print(best_sample)
    print(bin_vars_dict)

    enz_damag, cmp_damag, res_graph, end_filepath = modify_dot_graph(bin_vars_dict, best_sample, file_path, marked_nodes)

    print(f'Enzyme Damage: {enz_damag}')
    print(f'Compound Damage: {cmp_damag}')
    # Write the modified graph back to the new .dot file
    res_graph.del_node('\"\\r\\n\"')
    res_graph.write(end_filepath)
    print(qubo_matrix)


def solver_manager(file_path):
    # file_path = "/workspaces/graph-coloring/Glycine_serine_threonine_test.dot"
    qubo_matrix, graph, bin_vars_dict, marked_nodes = get_qubo_dict_for_dwave_qa(file_path)
    # "C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\smallDots\\eco_filtering_dot"
    # "\\Glycerolipid_marked.dot")
    # "C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\smallDots_w_marked_and_tests"
    # "\\eco_filtering_dot\\Polyketide.dot")

    bqm = BinaryQuadraticModel.from_qubo(qubo_matrix)

    value = "simulated_anneal_qubo"
    # value = "dwave_sampler_fixed_empedding"
    #value = "dwave_quantum"
    # value = "hybrid"
    # value = "exact_bqm"

    best_sample = switch_case_if(value, bqm, qubo_matrix)
    #print(best_sample)
    #print(bin_vars_dict)

    enz_damag, cmp_damag, res_graph, end_filepath = modify_dot_graph(bin_vars_dict, best_sample, file_path, marked_nodes)

    print(f'Enzyme Damage: {enz_damag}')
    print(f'Compound Damage: {cmp_damag}')
    # Write the modified graph back to the new .dot file
    res_graph.del_node('\"\\r\\n\"')
    res_graph.write(end_filepath)
    #print(qubo_matrix)


# ------- Main program -------
if __name__ == "__main__":
    file_path = "C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\largeDots_w_marked_and_tests\\hsa_filtering_dot\\Nucleotide metabolism_marked.dot"
    solver_manager(file_path)
    # solver_manager(file_path)

