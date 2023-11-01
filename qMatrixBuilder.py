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


# Takes out the "outtake" of a string ("stringer") and also other stuff to create the tuples for filling in the matrix.
def term_splicer(outtake, stringer):
    strLen = len(stringer)
    stringer = stringer.replace("*" + outtake, "")
    if strLen == len(stringer):
        stringer = stringer.replace(outtake, "")

    stringer = stringer.replace("+", "")

    return stringer


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
    return np_matrix_of_tuples_and_size_generator(mtrx_tuples,len(variable_list))


# This code is executed only when the script is run directly
if __name__ == '__main__':
    stringo = (
        "-2*C1*k3*x_9 + C1*k3 - 2*C2*k3*x_10 + C2*k3 - 2*C3*k3*x_11 + C3*k3 - 2*C4*k3*x_12 + C4*k3 - 2*C5*k3*x_13 + "
        "C5*k3 - 2*C6*k3*x_14 + C6*k3 - 2*C7*k3*x_15 + C7*k3 + k1*x_1 + k1*x_2 + k1*x_3 + k1*x_4 + k1*x_5 + k1*x_6 "
        "+ k1*x_7 - k2*x_8 + k2 + k3*x_10 + k3*x_11 + k3*x_12 + k3*x_13 + k3*x_14 + k3*x_15 + k3*x_9 - "
        "2*k4*x_10*x_3 + k4*x_10 - 2*k4*x_11*x_4 + k4*x_11 - 2*k4*x_12*x_5 + k4*x_12 - 2*k4*x_13*x_6 + k4*x_13 - "
        "2*k4*x_14*x_7 + k4*x_14 - 2*k4*x_15*x_8 + k4*x_15 - 2*k4*x_2*x_9 + k4*x_2 + k4*x_3 + k4*x_4 + k4*x_5 + "
        "k4*x_6 + k4*x_7 + k4*x_8 + k4*x_9")
    var_list = ["x_1", "x_2", "x_3", "x_4", "x_5", "x_6", "x_7", "x_8", "y_r1_1", "y_r1", "y_r2", "y_c1"]

    prnt_gnrtd_ltx_mtrx(var_list, stringo)
