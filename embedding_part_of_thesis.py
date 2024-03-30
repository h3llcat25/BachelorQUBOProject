from dwave.system import DWaveSampler, FixedEmbeddingComposite
from minorminer import find_embedding

# Assume qubos_variations is a list of lists, with each sublist containing 7 variations of a QUBO
qubos_variations = [
    # [qubo_variation_1, qubo_variation_2, ..., qubo_variation_7], for the first base QUBO
    # Other base QUBOs and their variations...
]

# Initialize the D-Wave sampler
sampler = DWaveSampler()

for idx, variations in enumerate(qubos_variations, start=1):
    all_interactions = set().union(*[set(qubo.keys()) for qubo in variations])
    embedding = find_embedding(all_interactions, sampler.edgelist)
    fixed_sampler = FixedEmbeddingComposite(sampler, embedding)

    # Store best samples and their energies for each pack of 7 QUBOs
    best_samples = []

    for qubo in variations:
        response = fixed_sampler.sample_qubo(qubo, label='QUBO Variation')
        best_sample = min(response.data(['sample', 'energy']), key=lambda x: x.energy)
        best_samples.append(best_sample)

    # Write the results to a file
    filename = f'qubo_pack_{idx}_results.txt'
    with open(filename, 'w') as file:
        for sample_num, best_sample in enumerate(best_samples, start=1):
            file.write(f"Best Sample for QUBO {sample_num}: {best_sample.sample}\n")
            file.write(f"Energy: {best_sample.energy}\n\n")

    print(f"Results for pack {idx} written to {filename}")
