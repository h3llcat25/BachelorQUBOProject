import numpy as np
from gurobi_optimods.qubo import solve_qubo

# Define the weights for your problem
# This is a 3x3 matrix as an example, you should replace this with your own weights
Q = np.array([[0, -1, -2], [0, -3, 3], [0, 0, 2]])

# Solve the QUBO problem
result = solve_qubo(Q)

# Print the solution
print(result.solution)
print(result.objective_value)
