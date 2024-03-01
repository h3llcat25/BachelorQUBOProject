import re
import numpy as np


# Takes a list of (binary variables) and the whole qubo optimization term and returns a list of tuples,
# that each contain a unit of the problem formulation. These tuples are the units that are separated by + and - signs.
def build_matrix_tuples(vari_list, term_string):
    reversed_var_list = vari_list.copy()
    reversed_var_list.reverse()

    mathTerms = re.split(r'(?=[+-])', term_string)
    # Initialize the output list
    output = []

    # Iterate over the mathTerms list. Each element is a by a + or - sign isolated part of the whole objective fct
    for mathTerm in mathTerms:
        if not mathTerm:
            continue
        # Initialize a tuple for each mathTerm
        tuple_elements = []

        # Iterate over the variable list to take out the variables and leave only constants or minus signs left.
        for var in reversed_var_list:
            if len(tuple_elements) == 2:  # break out, if there are already two variables in the tuple
                break
            # Check if the mathTerm contains the current var
            if var in mathTerm:
                # If it does, find the index of var in var_list and add to the tuple
                index = vari_list.index(var)
                tuple_elements.insert(0, index)

                # Remove the string "*" + var and all "+" characters from the mathTerm
                mathTerm = term_splicer(var, mathTerm.replace(" ", ""))

        # If there is a linear term, then this part squares it to get to the quadratic form
        if len(tuple_elements) == 1:
            tuple_elements.append(tuple_elements[0])

        # "mathTerm" is what is left, when removing all variables, whitespaces and "+"-signs. If there is nothing
        # left or only a "-"-sign, puts a 1
        if not mathTerm or mathTerm == "-":
            mathTerm = mathTerm + "1"

        # Add the modified mathTerm to the tuple
        tuple_elements.append(f' + ({mathTerm})')

        # Convert the tuple_elements list to a tuple and add to the output list
        output.append(tuple(tuple_elements))

    return output


def parse_polynomial(variables, poly):
    # Regular expression to match polynomial terms, capturing the sign, coefficient, and variable parts
    pattern = r'([+-]?\s*\d*)(\*?[a-zA-Z_][a-zA-Z0-9_]*)?'

    # Find all matches of the pattern in the polynomial string
    matches = re.findall(pattern, poly)

    # Initialize an empty list to hold the formatted terms
    terms = []

    for coeff, var in matches:
        # Remove whitespace and handle the case where the coefficient is missing or is just a sign
        coeff = coeff.replace(" ", "")
        if coeff == "+" or coeff == "" or coeff == "-":
            coeff += "1"  # Append '1' to make the coefficient explicit

        # Construct the term by combining the coefficient and variable, if present
        term = coeff + ('*' + var if var else '')
        terms.append(term)

    # Sort variables by length in descending order to match longer variables first
    variables = sorted(variables, key=len, reverse=True)

    # Initialize a list to hold tuples for each term and a variable to sum constant terms
    constant_sum = 0

    # Split the polynomial into terms
    for term in terms:
        # Initialize the coefficient and variables for this term
        term_vars = []

        # Check if the term has a coefficient
        if '*' in term:
            coeff, var_part = term.split('*')

            # Match variables in the term
            for var in variables:
                if var in var_part:
                    term_vars.append(var)
                    var_part = var_part.replace(var, '', 1)  # Remove the matched variable from the term
        else:
            # The term might be a constant or a single variable
            is_constant = True
            for var in variables:
                if var in term:
                    term_vars.append(var)
                    is_constant = False
                    break  # Stop after the first match to avoid matching subparts of variables

            if is_constant:
                constant_sum += float(term)  # Add the term to the constant sum
                continue  # Skip to the next term

        # Append the tuple for this term to the list
        terms.append((coeff, *term_vars))

    # If there were constant terms, add their sum as a tuple
    if constant_sum != 0:
        terms.append((constant_sum,))

    return terms


# Takes out the "outtake" of a string ("stringer") and also other stuff to create the tuples for filling in the matrix.
def term_splicer(outtake, stringer):
    strLen = len(stringer)
    stringer = stringer.replace("*" + outtake, "")
    if strLen == len(stringer):
        stringer = stringer.replace(outtake, "")

    stringer = stringer.replace("+", "")

    return stringer


def add_damage_matrix_to_q_matrix(mat1, mat2):
    # Check if both matrices have the same dimensions
    if mat1.shape != mat2.shape:
        raise ValueError("Both matrices should have the same dimensions")

    # Extract the diagonal of the second matrix
    diag2 = np.copy(np.diag(mat2))

    # Replace zero values in the diagonal with the corresponding values from the first matrix
    diag2[diag2 == 0] = np.diag(mat1)[diag2 == 0]

    # Add the modified diagonal of the second matrix to the first matrix
    mat1 += np.diag(diag2)

    return mat1


def generate_ltx_matrix(tuples, matrix_size):
    matrix = [["0" for _ in range(matrix_size)] for _ in range(matrix_size)]
    constants = ""

    for tup in tuples:
        if len(tup) == 1:
            # Append the content of the tuple (as a string in brackets) and a " + " string to the output string if
            # the size is 1
            constants += "({}) + ".format(tup[0])
        elif len(tup) == 3:
            if matrix[tup[0]][tup[1]] == "0":
                matrix[tup[0]][tup[1]] = tup[2].lstrip(" + ")
            else:
                matrix[tup[0]][tup[1]] += tup[2]
        else:
            return f"Error: encountered tuple of unsupported size: {tup}"

    latex_matrix = "\\begin{equation}\n\\begin{bmatrix}"
    for i in range(matrix_size):
        latex_matrix += "\n" + " & ".join(str(matrix[i][j]) for j in range(matrix_size))
        if i < matrix_size - 1:
            latex_matrix += " \\\\"
    latex_matrix += "\n\\end{bmatrix}\n\\end{equation}"
    return latex_matrix


# Only possible with integer values
def np_matrix_of_tuples_and_size_generator(tuples, matrix_size):
    matrix = np.array([[0 for _ in range(matrix_size)] for _ in range(matrix_size)])
    constants = 0

    for tup in tuples:
        if len(tup) == 1:
            tempConst = tup[0].replace("+", "").replace("(", "").replace(")", "").replace(" ", "")
            constants += int(tempConst)
        elif len(tup) == 3:
            if matrix[tup[0]][tup[1]] == 0:
                tempConst = tup[2].replace("+", "").replace("(", "").replace(")", "").replace(" ", "")
                matrix[tup[0]][tup[1]] = int(tempConst)
            else:
                matrix[tup[0]][tup[1]] += int(tup[2])
        else:
            return f"Error: encountered tuple of unsupported size: {tup}"

    # Save the matrix to a text file
    # np.savetxt('matrixo.txt', matrix, fmt='%d')
    return matrix


# Print generated latex matrix
def prnt_gnrtd_ltx_mtrx(variable_list, objective_fct):
    mtrx_tuples = build_matrix_tuples(variable_list, objective_fct)
    print(generate_ltx_matrix(mtrx_tuples, len(variable_list)))


def generate_final_np_matrix(bin_var_dict, objective_fct):
    variable_list = list(bin_var_dict.values())
    mtrx_tuples = build_matrix_tuples(variable_list, objective_fct)
    return np_matrix_of_tuples_and_size_generator(mtrx_tuples, len(variable_list))


def generate_final_np_matrix_var_list(variable_list, objective_fct):
    mtrx_tuples = build_matrix_tuples(variable_list, objective_fct)
    return np_matrix_of_tuples_and_size_generator(mtrx_tuples, len(variable_list))


# This code is executed only when the script is run directly
if __name__ == '__main__':
    stringo = (
        # "-2*C1*k3*x_9 + C1*k3 - 2*C2*k3*x_10 + C2*k3 - 2*C3*k3*x_11 + C3*k3 - 2*C4*k3*x_12 + C4*k3 - 2*C5*k3*x_13 + "
        # "C5*k3 - 2*C6*k3*x_14 + C6*k3 - 2*C7*k3*x_15 + C7*k3 + k1*x_1 + k1*x_2 + k1*x_3 + k1*x_4 + k1*x_5 + k1*x_6 "
        # "+ k1*x_7 - k2*x_8 + k2 + k3*x_10 + k3*x_11 + k3*x_12 + k3*x_13 + k3*x_14 + k3*x_15 + k3*x_9 - "
        # "2*k4*x_10*x_3 + k4*x_10 - 2*k4*x_11*x_4 + k4*x_11 - 2*k4*x_12*x_5 + k4*x_12 - 2*k4*x_13*x_6 + k4*x_13 - "
        # "2*k4*x_14*x_7 + k4*x_14 - 2*k4*x_15*x_8 + k4*x_15 - 2*k4*x_2*x_9 + k4*x_2 + k4*x_3 + k4*x_4 + k4*x_5 + "
        # "k4*x_6 + k4*x_7 + k4*x_8 + k4*x_9")
        "-2*x_9*y_c1 + 1 - 2*x_10 + 1 - 2*x_11 + 2*x_12*y_r1 + - 2*x_13 + 4 - 2*x_14 - 2*x_15*y_r2 + x_1 + x_2 + x_3 + x_4 + x_5 + x_6 "
        "+ x_7 - x_8 + 4 + x_10 + x_11 + x_12 + x_13 + x_14 + x_15 + x_9 - "
        "2*x_3 + x_10 - 2x_11*x_4 + x_11 - 2*x_12*x_5 + x_12 - 2*x_13*x_6 + x_13 - "
        "2*x_14*x_7 + x_14 - 2*x_15*x_8 + 4*x_15 - 2*x_2*x_9 + x_2 + x_3 + x_4 + x_5 + "
        "x_6 + x_7 + x_8*y_r1_1 + x_9")
    var_list = ["x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "x_9", "y_r1_1", "y_r1", "y_r2", "y_c1"]

    # print(generate_final_np_matrix_var_list(var_list, stringo))
    print(parse_polynomial(var_list, stringo))
