import numpy as np

# Dictionary with data
Q = {(0, 0): 1, (1, 1): 1, (0, 1): 1}

# Determine the dimensions of the matrix
max_row = max(key[0] for key in Q.keys())
max_col = max(key[1] for key in Q.keys())
num_rows = max_row + 1
num_cols = max_col + 1

# Create a NumPy array and initialize it with zeros
matrix = np.zeros((num_rows, num_cols), np.int8)
print(matrix)

# Populate the matrix with values from the dictionary
for key, value in Q.items():
    if value != 0:
        matrix[key] = value

print(matrix)


Q_dict = {index: value for index, value in np.ndenumerate(matrix)}

print(Q_dict)