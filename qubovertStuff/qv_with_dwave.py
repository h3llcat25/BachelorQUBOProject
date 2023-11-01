import qubovert as qv
from qubovert import boolean_var

N = 10

# create the variables
x = {i: boolean_var('x(%d)' % i) for i in range(N)}

# minimize \sum_{i=0}^{N-2} (1-2x_{i}) x_{i+1}
model = 0
for i in range(N-1):
    model += (1 - 2 * x[i]) * x[i+1]

# subject to the constraint that x_1 equals the XOR of x_3 and x_5
# enforce with a penalty factor of 3
model.add_constraint_eq_XOR(x[1], x[3], x[5], lam=3)

model_solution = model.solve_bruteforce()
print("Variable assignment:", model_solution)
print("Model value:", model.value(model_solution))
print("Constraints satisfied?", model.is_solution_valid(model_solution))

from qubovert.sim import anneal_pubo

res = anneal_pubo(model, num_anneals=10)
model_solution = res.best.state

print("Variable assignment:", model_solution)
print("Model value:", res.best.value)
print("Constraints satisfied?", model.is_solution_valid(model_solution))

from neal import SimulatedAnnealingSampler

# Get the QUBO form of the model
qubo = model.to_qubo()

# D-Wave accept QUBOs in a different format than qubovert's format
# to get the qubo in this form, use the .Q property
dwave_qubo = qubo.Q

# solve with D-Wave
res = SimulatedAnnealingSampler().sample_qubo(dwave_qubo, num_reads=10)
qubo_solution = res.first.sample

# convert the qubo solution back to the solution to the model
model_solution = model.convert_solution(qubo_solution)

print("Variable assignment:", model_solution)
print("Model value:", model.value(model_solution))
print("Constraints satisfied?", model.is_solution_valid(model_solution))

