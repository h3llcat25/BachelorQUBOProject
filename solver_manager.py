from dimod import BinaryQuadraticModel, ExactSolver
from dwave.system import DWaveSampler, EmbeddingComposite, FixedEmbeddingComposite
from dwave.samplers import SimulatedAnnealingSampler
from minorminer import find_embedding

from param_and_solve_automation import *


def solve_w_exact_solver(bqm):
    solver = ExactSolver()
    return solver.sample(bqm).first.sample


def solve_w_simulated_annealing_sampler(qubo):
    sampler = SimulatedAnnealingSampler()

    # Solve the QUBO
    sample_set = sampler.sample_qubo(qubo)
    return sample_set.first.sample


def solve_w_dwave_sampler_fixed_empedding(q_matrix):
    # Initialize the D-Wave sampler
    sampler = DWaveSampler()
    interactions = set(qubo_matrix.keys())
    embedding = find_embedding(interactions, sampler.edgelist)
    fixed_sampler = FixedEmbeddingComposite(sampler, embedding)

    # Store best samples and their energies for each pack of 7 QUBOs
    response = fixed_sampler.sample_qubo(q_matrix)
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
    num_runs = 100  # Kann angepasst werden, um optimale Ergebnisse zu finden

    # Simulated Annealing durchführen
    response = sampler.sample_qubo(Q, num_reads=num_runs)

    # Ergebnisse ausgeben
    for datum in response.data(['sample', 'energy', 'num_occurrences']):
        print(datum.sample, "Energie:", datum.energy, "Anzahl Vorkommen:", datum.num_occurrences)

    # Das optimale Ergebnis ermitteln
    best_sample = response.first.sample
    best_energy = response.first.energy

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
    elif value == 'c':
        return "Third case"
    else:
        return "Default case"


# ------- Main program -------
if __name__ == "__main__":
    file_path = "/workspaces/graph-coloring/test.dot"
    qubo_matrix, graph, bin_vars_dict, marked_nodes = get_qubo_dict_for_dwave_qa(file_path)
    # "C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\smallDots\\eco_filtering_dot"
    # "\\Glycerolipid_marked.dot")
    # "C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\smallDots_w_marked_and_tests"
    # "\\eco_filtering_dot\\Polyketide.dot")

    bqm = BinaryQuadraticModel.from_qubo(qubo_matrix)

    # value = "simulated_anneal_qubo"
    # value = "dwave_sampler_fixed_empedding"
    value = "dwave_quantum"
    # value = "exact_bqm"

    best_sample = switch_case_if(value, bqm, qubo_matrix)
    print(best_sample)

    node_name_to_solulu_dict = node_to_solutionsValue_assignment_for_result(bin_vars_dict, best_sample)

    modify_dot_graph(graph, node_name_to_solulu_dict, file_path, marked_nodes)
    # graph, solution_dictionary, dot_file, target_nodes
    # modify_dot_graph(graph, solution_dictionary, dot_file, target_nodes)

    # Print the results
    # print("Solutions:")
    # for sample, energy in result.data(['sample', 'energy']):
    #     print(f"Sample: {sample}, Energy: {energy}")

    # print(bqm)
    # cqm = build_cqm(G, num_colors)

    # sample = run_hybrid_solver(cqm)
