import sympy
from collections import defaultdict


'''
def create_qubo_dict(bin_vars_dict, objective_function, damage_objective=None):
    bin_vars_list = list(bin_vars_dict.values())
    # Create a dictionary of binary variables
    binary_vars = {f'{i}': Binary(f'{i}') for i in bin_vars_list}

    # if damage_objective:

    model = objective_function.compile()

    # Step 6: Get the QUBO
    Q, offset = model.to_qubo()

    # Q is the QUBO dictionary
    print("QUBO dictionary:", Q)

    # return Q
'''


def creating_qubo_dict_sympy(modified_opt_term, modified_damage_opt_term=None):
    q_dict = defaultdict(lambda: 0)

    coeff_dict = modified_opt_term.as_coefficients_dict()
    '''
    # For Testing Purposes
    print("Coefficients dictionary:")
    for term, coeff in coeff_dict.items():
        print(f"{term}: {coeff}")
    '''

    for term, coeff in coeff_dict.items():
        # Normalize squared terms
        term = str(term)
        normalized_term = term.replace('**2', '')

        # Split the normalized term by '*' to separate variables
        variables = normalized_term.split('*')

        if len(variables) == 2:
            # Quadratic term: add to QUBO with variable names as keys
            q_dict[(variables[0], variables[1])] += coeff
        elif len(variables) == 1 and variables[0] != '1':
            # Linear or squared term: add to QUBO, using the same variable name for both keys
            q_dict[(variables[0], variables[0])] += coeff
        elif variables[0] != '1':
            print(f"There is something fishy here with {variables}")

    if modified_damage_opt_term:
        coeff_dict_damage = modified_damage_opt_term.as_coefficients_dict()

        for term, coeff in coeff_dict_damage.items():
            term = str(term)

            normalized_term = term.replace('**2', '')
            variables = normalized_term.split('*')

            if len(variables) == 2:
                q_dict[(variables[0], variables[1])] += coeff
            elif len(variables) == 1 and variables[0] != '1':
                q_dict[(variables[0], variables[0])] += coeff
            elif variables[0] != '1':
                print(f"There is something fishy here with {variables}")

    # Convert defaultdict back to a regular dict for the final output
    qubo_dict = dict(q_dict)

    return qubo_dict


def main():
    str_expr = "x_1*y_r3*24 + x_1*24*y_r3 -  x_1*y_r3_1*24 + y_r3*x_1*24 -1 + y_r3**2 + y_r3*y_r3  + y_r3"
    modified_opt_term = sympy.sympify(str_expr)
    coeff_dict = modified_opt_term.as_coefficients_dict()

    print("Coefficients dictionary:")
    for term, coeff in coeff_dict.items():
        print(f"{term}: {coeff}")


if __name__ == '__main__':
    main()