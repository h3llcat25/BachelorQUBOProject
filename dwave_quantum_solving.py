import dwave_qbsolv

# Define the QUBO problem. Here, the QUBO problem is represented as a dictionary
Q = {(0, 0): 1, (1, 1): 1, (0, 1): 1}

# Use the QBSolv module to solve the QUBO problem
response = QBSolv().sample_qubo(Q)

# Print the samples and energies
print("samples=" + str(list(response.samples())))
print("energies=" + str(list(response.data_vectors['energy'])))





