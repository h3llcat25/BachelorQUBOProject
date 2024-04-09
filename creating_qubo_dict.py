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


def creating_qubo_dict_sympy(bin_var_dict, objective_expr, damage_expression=None): # (bin_vars_dict, modified_opt_term,
    # modified_damage_opt_term) bin_vars_dict hat die namen der nodes und die variablen, die anderen beiden haben einfach expressions.

    og_node_names = True
    var_bin_dict = {value: key for key, value in bin_var_dict.items()}

    # Convert the sympy expression to a dictionary
    objective_dict_whole = objective_expr.as_coefficients_dict()
    if damage_expression is not None:
        objective_dict_damage = damage_expression.as_coefficients_dict()

        # Merge the two dictionaries
        for term, coeff in objective_dict_damage.items():
            if term in objective_dict_whole:
                objective_dict_whole[term] += coeff
            else:
                objective_dict_whole[term] = coeff

    q_dict = defaultdict(lambda: 0)

    for term, coeff in objective_dict_whole.items():
        if isinstance(term, sympy.Symbol):  # linear term
            if og_node_names:
                q_dict[var_bin_dict[(str(term))], var_bin_dict[(str(term))]] += coeff
            else:
                q_dict[(str(term), str(term))] += coeff
        elif len(term.args) == 2:  # quadratic term
            var1, var2 = map(str, term.args)
            if var2 == "2":
                if og_node_names:
                    q_dict[var_bin_dict[var1], var_bin_dict[var1]] += coeff
                else:
                    q_dict[(var1, var1)] += coeff
            else:
                if og_node_names:
                    q_dict[var_bin_dict[var1], var_bin_dict[var2]] += coeff
                else:
                    q_dict[(var1, var2)] += coeff
        else:  # constant term
            continue


    # Convert defaultdict back to a regular dict for the final output
    qubo_dict = dict(q_dict)

    return qubo_dict


def test_check_q_creation(objective_expr, damage_expression=None):  # (bin_vars_dict, modified_opt_term,
    # modified_damage_opt_term) bin_vars_dict hat die namen der nodes und die variablen, die anderen beiden haben einfach expressions.

    # Convert the sympy expression to a dictionary
    objective_dict_whole = objective_expr.as_coefficients_dict()
    if damage_expression is not None:
        objective_dict_damage = damage_expression.as_coefficients_dict()

        # Merge the two dictionaries
        for term, coeff in objective_dict_damage.items():
            if term in objective_dict_whole:
                objective_dict_whole[term] += coeff
            else:
                objective_dict_whole[term] = coeff

    q_dict = defaultdict(lambda: 0)

    for term, coeff in objective_dict_whole.items():
        if isinstance(term, sympy.Symbol):  # linear term
            q_dict[(str(term), str(term))] += coeff
        elif len(term.args) == 2:  # quadratic term
            var1, var2 = map(str, term.args)
            if var2 == "2":
                q_dict[(var1, var1)] += coeff
            else:
                q_dict[(var1, var2)] += coeff
        else:  # constant term
            continue

    # Convert defaultdict back to a regular dict for the final output
    qubo_dict = dict(q_dict)

    test_output = True

    if test_output:
        return qubo_dict, objective_dict_whole

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