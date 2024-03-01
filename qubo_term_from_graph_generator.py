import math

import networkx as nx


# Step 1: Read the .dot file and extract the graph.
def read_dot_file(file_path):
    try:
        graph = nx.drawing.nx_agraph.read_dot(file_path)
        return graph
    except Exception as e:
        print(f'Error reading the .dot file: {str(e)}')
        return None


# Step 2: Check the input format.
def check_input_format(graph):
    try:
        if graph:
            # Check if the graph is directed
            if not graph.is_directed():
                print("Error: The graph is not directed.")
                return False
            return True
        else:
            return False
    except Exception as e:
        print(f'Error checking the input format: {str(e)}')
        return False


# Step 3: Create binary variables for each node and save nodes depending on type and color
def create_binary_variables(graph):
    binary_var_dict = {}
    marked_c_nodes = []
    c_nodes = []
    r_nodes = []
    node_index = 1

    try:
        for node, data in graph.nodes(data=True):
            # Error if the node has more or less than 2 attributes (type and color)
            if len(data) < 1 or len(data) > 2:
                print(f"Node {node} has invalid attribute count. Found attributes: {data}")
                return None

            node_type = data.get('type')
            color = data.get('color')

            # Error if the node has no type
            if not node_type:
                print(f"Node {node} has no type attribute")
                return None

            if node_type not in ["E", "C", "R"]:
                print(f"Node {node} has invalid type {node_type}")
                return None

            # Error if a node has color attribute, and it's not type C
            if color and node_type != "C":
                print(f"Node {node} has color attribute but is not type C")
                return None

            if color == "red" and node_type == "C":
                marked_c_nodes.append(node)

            if node_type == "C":
                c_nodes.append(node)
            if node_type == "R":
                r_nodes.append(node)

            binary_var_dict[node] = f'x_{node_index}'
            node_index += 1

    except Exception as e:
        print(f'Error creating binary variables: {str(e)}')
        return None
    return binary_var_dict, marked_c_nodes, c_nodes, r_nodes


# Step 3(a): Create binary variables for each node and save nodes depending on type and color, BUT also Enzyme List
# output.
def create_binary_variables_plus_enzymes(graph):
    binary_var_dict = {}
    marked_c_nodes = []
    c_nodes = []
    r_nodes = []
    e_nodes = []
    node_index = 1

    try:
        for node, data in graph.nodes(data=True):
            # Error if the node has more or less than 2 attributes (type and color)
            if len(data) < 1 or len(data) > 2:
                print(f"Node {node} has invalid attribute count. Found attributes: {data}")
                return None

            node_type = data.get('type')
            color = data.get('color')

            # Error if the node has no type
            if not node_type:
                print(f"Node {node} has no type attribute")
                return None

            if node_type not in ["E", "C", "R"]:
                print(f"Node {node} has invalid type {node_type}")
                return None

            # Error if a node has color attribute, and it's not type C
            if color and node_type != "C":
                print(f"Node {node} has color attribute but is not type C")
                return None

            if color == "red" and node_type == "C":
                marked_c_nodes.append(node)

            if node_type == "E":
                e_nodes.append(node)
            if node_type == "C":
                c_nodes.append(node)
            if node_type == "R":
                r_nodes.append(node)

            binary_var_dict[node] = f'x_{node_index}'
            node_index += 1

    except Exception as e:
        print(f'Error creating binary variables: {str(e)}')
        return None
    return binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes


# Step 4: Create the equation, where you can later decide penalty values manually.
def generate_equation_manual_penalties(graph, binary_var_dict, marked_c_nodes, c_nodes, r_nodes):
    try:
        # start from where to count for the extra variables for reactions r and compounds c.
        extra_var_index_r = 1
        extra_var_index_c = 1
        equation = ''
        react_pena_term = ''
        compo_pena_term = ''

        # Case 1: Damage Term counting how many un-targeted compounds got inhibited
        if len(c_nodes) - len(marked_c_nodes) > 0:
            equation += "k1*("

        for c_node in c_nodes:
            if c_node not in marked_c_nodes:
                equation += f'{binary_var_dict[c_node]} + '

        if equation:
            equation = equation.rstrip(' + ')
            equation += ") + "

        # Case 2  Target Term, how many targeted compounds didn't get inhibited.

        if len(marked_c_nodes) > 0:
            equation += "k2*("
            for marked_node in marked_c_nodes:
                equation += f'(1 - {binary_var_dict[marked_node]}) + '
            equation = equation.rstrip(' + ')
            equation += ") + "

        # Case 3 Reactions and Reduction Penalty:
        if len(r_nodes) > 0:
            equation += "k3*("
            for r_node in r_nodes:  # For all Reaction nodes
                if len(list(graph.predecessors(r_node))) > 0:  # If the number of predecessors is bigger than zero.
                    currBinVar_r = binary_var_dict[r_node]  # currBinVar contains r_node
                    if len(list(graph.predecessors(r_node))) == 1:  # If r_node has just one predecessor
                        equation += f'{binary_var_dict[list(graph.predecessors(r_node))[0]]} + {currBinVar_r}'
                        equation += f' - 2 * {binary_var_dict[list(graph.predecessors(r_node))[0]]} * {currBinVar_r} + '
                    else:  # If there are more Predecessors than 1:
                        equation += f'(1 - {currBinVar_r} - y_r{extra_var_index_r}'
                        equation += f' + 2 * {currBinVar_r} * y_r{extra_var_index_r}) + '
                        binary_var_dict[f'extra Var for r{extra_var_index_r}'] = f'y_r{extra_var_index_r}'
                        react_pena_term += reaction_reduction_penalty(list(graph.predecessors(r_node)), binary_var_dict,
                                                                      extra_var_index_r)
                        extra_var_index_r += 1
            equation = equation.rstrip(' + ')
            equation += ") + k4*("

        # Case 4 Compound and Reduction Penalty
        for c_node in c_nodes:
            if len(list(graph.predecessors(c_node))) != 0:
                currBinVar_c = binary_var_dict[c_node]
                if len(list(graph.predecessors(c_node))) == 1:
                    fill4Var = binary_var_dict[list(graph.predecessors(c_node))[0]]
                    equation += f'({currBinVar_c} + {fill4Var} - 2 * {currBinVar_c} * {fill4Var}) + '
                else:
                    equation += f'({currBinVar_c} + y_c{extra_var_index_c} - 2 * {currBinVar_c} * y_c{extra_var_index_c}) + '
                    binary_var_dict[f'extra Var for c{extra_var_index_c}'] = f'y_c{extra_var_index_c}'
                    compo_pena_term += compound_reduction_penalty(list(graph.predecessors(c_node)), binary_var_dict,
                                                                  extra_var_index_c)
                    extra_var_index_c += 1
        equation = equation.rstrip(' + ')
        equation += ") + "

        if react_pena_term:
            react_pena_term = react_pena_term.rstrip(' + ')
            equation += f'k5*({react_pena_term}) + '
        if compo_pena_term:
            compo_pena_term = compo_pena_term.rstrip(' + ')
            equation += f'k6*({compo_pena_term}) + '
        equation = equation.rstrip(' + ')
        # print(f'Generated equation: {equation}')
        return binary_var_dict, equation

    except Exception as e:
        print(f'Error generating equation: {str(e)}')


# Step 4: Generate the equation with auto penalty
def generate_equation_auto_penalty(graph, binary_var_dict, marked_c_nodes, c_nodes, r_nodes):
    try:
        # start from where to count for the extra variables for reactions r and compounds c.
        extra_var_index_r = 1
        extra_var_index_c = 1
        equation = ''
        react_pena_term = ''
        compo_pena_term = ''

        # Automatic Penalty calculation, based on if all Reactions, Enzymes and marked Compounds of the metabolic
        # network are inhibited + 10%
        meta_rulebreak_penalty = math.floor(len(c_nodes) + len(marked_c_nodes))
        print(f'Meta Rulebreaker Penalty is: {meta_rulebreak_penalty}')

        # Damage Term has penalty factor 1 (not explicitly written out)
        for c_node in c_nodes:
            if c_node not in marked_c_nodes:
                equation += f'{binary_var_dict[c_node]} + '

        # Case 2  Target Term, how many targeted compounds didn't get inhibited.
        if len(marked_c_nodes) > 0:
            equation += f'{meta_rulebreak_penalty}*('
            for marked_node in marked_c_nodes:
                equation += f'(1 - {binary_var_dict[marked_node]}) + '

        # Case 3 Reactions and Reduction Penalty:
        if len(r_nodes) > 0:
            for r_node in r_nodes:  # For all Reaction nodes
                if len(list(graph.predecessors(r_node))) > 0:  # If the number of predecessors is bigger than zero.
                    currBinVar_r = binary_var_dict[r_node]  # currBinVar contains r_node
                    if len(list(graph.predecessors(r_node))) == 1:  # If r_node has just one predecessor
                        equation += f'{binary_var_dict[list(graph.predecessors(r_node))[0]]} + {currBinVar_r}'
                        equation += f' - 2 * {binary_var_dict[list(graph.predecessors(r_node))[0]]} * {currBinVar_r} + '
                    else:  # If there are more Predecessors than 1:
                        equation += f'(1 - {currBinVar_r} - y_r{extra_var_index_r}'
                        equation += f' + 2 * {currBinVar_r} * y_r{extra_var_index_r}) + '
                        binary_var_dict[f'extra Var for r{extra_var_index_r}'] = f'y_r{extra_var_index_r}'
                        react_pena_term += reaction_reduction_penalty(list(graph.predecessors(r_node)), binary_var_dict,
                                                                      extra_var_index_r)
                        extra_var_index_r += 1

        # Case 4 Compound and Reduction Penalty
        for c_node in c_nodes:
            if len(list(graph.predecessors(c_node))) != 0:
                currBinVar_c = binary_var_dict[c_node]
                if len(list(graph.predecessors(c_node))) == 1:
                    fill4Var = binary_var_dict[list(graph.predecessors(c_node))[0]]
                    equation += f'({currBinVar_c} + {fill4Var} - 2 * {currBinVar_c} * {fill4Var}) + '
                else:
                    equation += f'({currBinVar_c} + y_c{extra_var_index_c} - 2 * {currBinVar_c} * y_c{extra_var_index_c}) + '
                    binary_var_dict[f'extra Var for c{extra_var_index_c}'] = f'y_c{extra_var_index_c}'
                    compo_pena_term += compound_reduction_penalty(list(graph.predecessors(c_node)), binary_var_dict,
                                                                  extra_var_index_c)
                    extra_var_index_c += 1

        if react_pena_term:
            react_pena_term = react_pena_term.rstrip(' + ')
            equation += f'({react_pena_term}) + '
        if compo_pena_term:
            compo_pena_term = compo_pena_term.rstrip(' + ')
            equation += f'({compo_pena_term}) + '
        equation = equation.rstrip(' + ')
        equation += ")"
        # print(f'Generated equation: {equation}')
        return binary_var_dict, equation

    except Exception as e:
        print(f'Error generating equation: {str(e)}')


# Step 4: Generate the equation with auto penalty and Enzyme damage calculation 
def generate_equation_auto_penalty_and_enzyme_damage(graph, binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes):
    try:
        # start from where to count for the extra variables for reactions r and compounds c.
        extra_var_index_r = 1
        extra_var_index_c = 1
        equation = ''
        damage_only_equation = ''
        react_pena_term = ''
        compo_pena_term = ''

        double_edged_reactions = set()
        circular_marked_c_nodes = []

        # Automatic Penalty calculation, based on if all Reactions, Enzymes and marked Compounds of the metabolic
        # network are inhibited + 10%
        meta_rulebreak_penalty = int(math.floor(len(binary_var_dict)))
        print(f'Meta Rulebreaker Penalty is: {meta_rulebreak_penalty}')

        # Damage Term has penalty factor 1 (not explicitly written out)
        damage_only_equation += '2 * ('
        for c_node in c_nodes:
            if c_node not in marked_c_nodes:
                damage_only_equation += f'{binary_var_dict[c_node]} + '

        damage_only_equation = damage_only_equation.rstrip(' + ')
        damage_only_equation += ") + "

        # Damage Term has penalty factor 1 (not explicitly written out)
        for e_node in e_nodes:
            damage_only_equation += f'{binary_var_dict[e_node]} + '
        damage_only_equation = damage_only_equation.rstrip(' + ')

        # Case 2  Target Term, how many targeted compounds didn't get inhibited.
        if len(marked_c_nodes) > 0:
            equation += f'{meta_rulebreak_penalty}*('
            for marked_node in marked_c_nodes:
                equation += f'(1 - {binary_var_dict[marked_node]}) + '
                # check for circular subgraphs (reaction nodes connected with 2 edges in both directions to the same compound.
                predecessors_set = set(graph.predecessors(marked_node))
                successors_set = set(graph.successors(marked_node))
                if predecessors_set & successors_set:
                    double_edged_reactions.update(predecessors_set & successors_set)
                    circular_marked_c_nodes.append(marked_node)

            if double_edged_reactions:
                for r_node in double_edged_reactions:  # For all Reaction nodes
                    r_nodes.remove(r_node)
                    equation += f'(1 - {binary_var_dict[r_node]}) + '
                    non_marked_reaction_predecessors = [e_node for e_node in graph.predecessors(r_node) if graph.nodes[e_node]['type'] == "E"]
                    if non_marked_reaction_predecessors:
                        if len(non_marked_reaction_predecessors) == 1:  # If r_node has just one predecessor
                            equation += f'(1 - {binary_var_dict[non_marked_reaction_predecessors[0]]}) + '
                        elif len(non_marked_reaction_predecessors) == 2:
                            equation += f'(1 - {binary_var_dict[non_marked_reaction_predecessors[0]]}) * (1 - {binary_var_dict[non_marked_reaction_predecessors[1]]}) +'
                        else:  # If there are more Predecessors left than 1:
                            equation += f'(1 - {binary_var_dict[non_marked_reaction_predecessors[0]]} * y_r{extra_var_index_r}) +'
                            binary_var_dict[f'extra Var for r{extra_var_index_r}'] = f'y_r{extra_var_index_r}'
                            react_pena_term += reaction_reduction_penalty(non_marked_reaction_predecessors[1:],
                                                                          binary_var_dict,
                                                                          extra_var_index_r)
                            extra_var_index_r += 1

        # Case 3 Reactions and Reduction Penalty:
        if len(r_nodes) > 0:
            for r_node in r_nodes:  # For all Reaction nodes
                if len(list(graph.predecessors(r_node))) > 0:  # If the number of predecessors is bigger than zero.
                    currBinVar_r = binary_var_dict[r_node]  # currBinVar contains r_node
                    if len(list(graph.predecessors(r_node))) == 1:  # If r_node has just one predecessor
                        equation += f'{binary_var_dict[list(graph.predecessors(r_node))[0]]} + {currBinVar_r}'
                        equation += f' - 2 * {binary_var_dict[list(graph.predecessors(r_node))[0]]} * {currBinVar_r} + '
                    else:  # If there are more Predecessors than 1:
                        equation += f'(1 - {currBinVar_r} - y_r{extra_var_index_r}'
                        equation += f' + 2 * {currBinVar_r} * y_r{extra_var_index_r}) + '
                        binary_var_dict[f'extra Var for r{extra_var_index_r}'] = f'y_r{extra_var_index_r}'
                        react_pena_term += reaction_reduction_penalty(list(graph.predecessors(r_node)), binary_var_dict,
                                                                      extra_var_index_r)
                        extra_var_index_r += 1

        # Case 4 Compound and Reduction Penalty
        for c_node in c_nodes:
            if len(list(graph.predecessors(c_node))) != 0:
                currBinVar_c = binary_var_dict[c_node]
                if len(list(graph.predecessors(c_node))) == 1:
                    fill4Var = binary_var_dict[list(graph.predecessors(c_node))[0]]
                    equation += f'({currBinVar_c} + {fill4Var} - 2 * {currBinVar_c} * {fill4Var}) + '
                else:
                    equation += f'({currBinVar_c} + y_c{extra_var_index_c} - 2 * {currBinVar_c} * y_c{extra_var_index_c}) + '
                    binary_var_dict[f'extra Var for c{extra_var_index_c}'] = f'y_c{extra_var_index_c}'
                    compo_pena_term += compound_reduction_penalty(list(graph.predecessors(c_node)), binary_var_dict,
                                                                  extra_var_index_c)
                    extra_var_index_c += 1

        if react_pena_term:
            react_pena_term = react_pena_term.rstrip(' + ')
            equation += f'({react_pena_term}) + '
        if compo_pena_term:
            compo_pena_term = compo_pena_term.rstrip(' + ')
            equation += f'({compo_pena_term}) + '
        equation = equation.rstrip(' + ')
        equation += ")"
        # print(f'Generated equation: {equation}')
        return binary_var_dict, equation, damage_only_equation

    except Exception as e:
        print(f'Error generating equation: {str(e)}')


# Version that takes into account taking out as little enzymes as possible, and special casee for circles.
def generate_equation_auto_penalty_enzyme_reaction_circle(graph, binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes):
    try:
        # start from where to count for the extra variables for reactions r and compounds c.
        extra_var_index_r = 1
        extra_var_index_c = 1
        equation = ''
        damage_and_circle_react_equation = ''
        react_pena_term = ''
        compo_pena_term = ''

        # Automatic Penalty calculation, based on if all Reactions, Enzymes and marked Compounds of the metabolic
        meta_rulebreak_penalty = int(math.floor(len(binary_var_dict) - len(r_nodes)))
        print(f'Meta Rulebreaker Penalty is: {meta_rulebreak_penalty}')

        # Damage Term has penalty factor 2
        damage_and_circle_react_equation += f'2*('
        for c_node in c_nodes:
            if c_node not in marked_c_nodes:
                damage_and_circle_react_equation += f'{binary_var_dict[c_node]} + '
        damage_and_circle_react_equation = damage_and_circle_react_equation.rstrip(' + ')
        damage_and_circle_react_equation += ") + "

        # Damage Term has penalty factor 1 (not explicitly written out)
        for e_node in e_nodes:
            damage_and_circle_react_equation += f'{binary_var_dict[e_node]} + '

        # Case 2  Target Term, how many targeted compounds didn't get inhibited, and handle Circle_compound problem
        if len(marked_c_nodes) > 0:
            equation += f'{meta_rulebreak_penalty}*('
            for marked_node in marked_c_nodes:
                equation += f'(1 - {binary_var_dict[marked_node]}) + '

        # Case 3 Reactions and Reduction Penalty:
        if len(r_nodes) > 0:
            for r_node in r_nodes:  # For all Reaction nodes
                if len(list(graph.predecessors(r_node))) > 0:  # If the number of predecessors is bigger than zero.
                    currBinVar_r = binary_var_dict[r_node]  # currBinVar_r contains r_node
                    if len(list(graph.predecessors(r_node))) == 1:  # If r_node has just one predecessor
                        equation += f'{binary_var_dict[list(graph.predecessors(r_node))[0]]} + {currBinVar_r}'
                        equation += f' - 2 * {binary_var_dict[list(graph.predecessors(r_node))[0]]} * {currBinVar_r} + '
                    else:  # If there are more Predecessors than 1:
                        equation += f'(1 - {currBinVar_r} - y_r{extra_var_index_r}'
                        equation += f' + 2 * {currBinVar_r} * y_r{extra_var_index_r}) + '
                        binary_var_dict[f'extra Var for r{extra_var_index_r}'] = f'y_r{extra_var_index_r}'
                        react_pena_term += reaction_reduction_penalty(list(graph.predecessors(r_node)), binary_var_dict,
                                                                      extra_var_index_r)
                        extra_var_index_r += 1

                        # # Check for Circles:
                        # if is_circular(marked_node, graph):
                        #     set1 = set(graph.predecessors(marked_node))
                        #     set2 = set(graph.successors(marked_node))
                        #     circular_react_list = list(set1 & set2)
                        #     if circular_react_list:
                        #         for react_node in circular_react_list:
                        #             predec_list = list(graph.predecessors(react_node))
                        #             enz_predec_list = [node for node in predec_list if node.get('type') == "E"]
                        #             if len(enz_predec_list) > 1:
                        #         # Wenn ein circular compound. mind. ein knoten

        # Case 4 Compound and Reduction Penalty
        for c_node in c_nodes:
            if len(list(graph.predecessors(c_node))) != 0:
                currBinVar_c = binary_var_dict[c_node]
                if len(list(graph.predecessors(c_node))) == 1:
                    fill4Var = binary_var_dict[list(graph.predecessors(c_node))[0]]
                    equation += f'({currBinVar_c} + {fill4Var} - 2 * {currBinVar_c} * {fill4Var}) + '
                else:
                    equation += f'({currBinVar_c} + y_c{extra_var_index_c} - 2 * {currBinVar_c} * y_c{extra_var_index_c}) + '
                    binary_var_dict[f'extra Var for c{extra_var_index_c}'] = f'y_c{extra_var_index_c}'
                    compo_pena_term += compound_reduction_penalty(list(graph.predecessors(c_node)), binary_var_dict,
                                                                  extra_var_index_c)
                    extra_var_index_c += 1

        if react_pena_term:
            react_pena_term = react_pena_term.rstrip(' + ')
            equation += f'({react_pena_term}) + '
        if compo_pena_term:
            compo_pena_term = compo_pena_term.rstrip(' + ')
            equation += f'({compo_pena_term}) + '
        equation = equation.rstrip(' + ')
        equation += ")"

        damage_and_circle_react_equation = damage_and_circle_react_equation.rstrip(' + ')

        return binary_var_dict, equation, damage_and_circle_react_equation

    except Exception as e:
        print(f'Error generating equation: {str(e)}')


# Help functions for Reaction and Compound Order Reduction Penalty
# Version2: fixed the error with each time the Penalty term gets negated in further recursions.
def reaction_reduction_penalty(predecessor_list, binary_var_dict, extra_var_index_r):
    pen_r_equation = ''

    try:
        y_filler = binary_var_dict[predecessor_list[0]]
        for i in range(1, (len(predecessor_list))):
            secVal = binary_var_dict[predecessor_list[i]]
            if i == 1:
                y_filler = f'(1 - {y_filler})'
            if i == len(predecessor_list) - 1:
                pen_r_equation += f'({y_filler} * (1 - {secVal}) - 2 * {y_filler} * y_r{extra_var_index_r}'
                pen_r_equation += f' - 2 * (1 - {secVal}) * y_r{extra_var_index_r} + 3 * y_r{extra_var_index_r}) + '
            else:
                pen_r_equation += f'({y_filler} * (1 - {secVal}) - 2 * {y_filler} * y_r{extra_var_index_r}_{i}'
                pen_r_equation += f' - 2 * (1 - {secVal}) * y_r{extra_var_index_r}_{i} + 3 * y_r{extra_var_index_r}_{i}) +'
                binary_var_dict[f'extra Var for r{extra_var_index_r}, Nr. {i}'] = f'y_r{extra_var_index_r}_{i}'
                y_filler = f'y_r{extra_var_index_r}_{i}'

    except Exception as e:
        print(f'Error creating binary variables for reaction order reduction penalty: {str(e)}')
        return None
    return pen_r_equation


def compound_reduction_penalty(predecessor_list, binary_var_dict, extra_var_index_c):
    pen_c_equation = ''

    try:
        y_filler = binary_var_dict[predecessor_list[0]]
        for i in range(1, (len(predecessor_list))):
            secVal = binary_var_dict[predecessor_list[i]]
            if i == len(predecessor_list) - 1:
                pen_c_equation += f'({y_filler} * {secVal} - 2 * {y_filler} * y_c{extra_var_index_c}'
                pen_c_equation += f' - 2 * {secVal} * y_c{extra_var_index_c} + 3 * y_c{extra_var_index_c}) + '
            else:
                pen_c_equation += f'({y_filler} * {secVal} - 2 * {y_filler} * y_c{extra_var_index_c}_{i}'
                pen_c_equation += f' - 2 * {secVal} * y_c{extra_var_index_c}_{i} + 3 * y_c{extra_var_index_c}_{i}) + '
                binary_var_dict[f'extra Var for c{extra_var_index_c}, Nr. {i}'] = f'y_c{extra_var_index_c}_{i}'
                y_filler = f'y_c{extra_var_index_c}_{i}'

    except Exception as e:
        print(f'Error creating binary variables for reaction order reduction penalty: {str(e)}')
        return None
    return pen_c_equation


def generating_qubo_term_from_graph(file_path):
    graph = read_dot_file(file_path)

    if check_input_format(graph):
        dict_and_sorted_nodes = create_binary_variables(graph)
        if dict_and_sorted_nodes:
            binary_vars, output_term = generate_equation_auto_penalty(graph, dict_and_sorted_nodes[0],
                                                                      dict_and_sorted_nodes[1],
                                                                      dict_and_sorted_nodes[2],
                                                                      dict_and_sorted_nodes[3])
            return binary_vars, output_term

    print(f'generating the QUBO Term failed')


# Generate qubo term from graph, but generates two Objective functions, to lower the complexity. -> AKTUELL!
def generating_qubo_term_from_graph_two_part(file_path):
    graph = read_dot_file(file_path)

    if check_input_format(graph):
        print("hier")
        dict_and_sorted_nodes = create_binary_variables_plus_enzymes(graph)
        if dict_and_sorted_nodes:
            # binary_vars, output_term = generate_equation_auto_penalty(graph, dict_and_sorted_nodes[0],
            binary_vars, output_term, output_damage_only = generate_equation_auto_penalty_and_enzyme_damage(graph,
                                                                                                            dict_and_sorted_nodes[
                                                                                                                0],
                                                                                                            dict_and_sorted_nodes[
                                                                                                                1],
                                                                                                            dict_and_sorted_nodes[
                                                                                                                2],
                                                                                                            dict_and_sorted_nodes[
                                                                                                                3],
                                                                                                            dict_and_sorted_nodes[
                                                                                                                4])
            return binary_vars, output_term, output_damage_only

    print(f'generating the QUBO Term failed')



if __name__ == '__main__':
    print(generating_qubo_term_from_graph_two_part(
        "C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\smallDots\\eco_filtering_dot"
        "\\Glycerolipid_marked.dot")[0])
