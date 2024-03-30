import math
import pydotplus
import random
import time

from tie_graph_problem import TieGraphProblem
from tie_qubo_struct import TieQuboStruct

# Helpfunction for getting predecessors of a specific Node in the .Dot file, Wichtig, gibt ne list of strings zurück
def predecessors(graph, target_node_name):
    # Specify the node for which you want to find the predecessors
    # Initialize an empty list to hold the predecessors
    predecessors = []

    # Iterate over the edges in the graph
    for edge in graph.get_edges():
        # Check if the edge points to the target node
        if edge.get_destination() == target_node_name or edge.get_destination() == f'\"{target_node_name}\"':
            # If so, add the source node of the edge to the predecessors list
            predecessors.append(edge.get_source())

    # Now 'predecessors' contains the IDs of the nodes that are predecessors of 'target_node'
    return predecessors


# Helpfunction for getting successors of a specific Node in the .Dot file
def successors(graph, target_node_name):
    # Specify the node for which you want to find the successors

    # Initialize an empty list to hold the successors
    successors = []

    # Iterate over the edges in the graph
    for edge in graph.get_edges():
        # Check if the edge comes from the target node
        if edge.get_source() == target_node_name or edge.get_destination() == f'\"{target_node_name}\"':
            # If so, add the get_destination node of the edge to the successors list
            successors.append(edge.get_destination())

    # Now 'successors' contains the IDs of the nodes that are successors of 'target_node'
    return successors


# Step 1: Read the .dot file and extract the graph.
def read_dot_file_pydotplus(file_path):
    # Load the .dot file
    try:
        graph = pydotplus.graph_from_dot_file(file_path)
        return graph
    except Exception as e:
        print(f'Error reading the .dot file: {str(e)}')
        return None


# Step 2: Create binary variables for each node and save nodes depending on type and color
# If k is an integer Value, choose k random C Nodes as Target nodes
# If seed not
def create_binary_variables_plus_enzymes(graph, k=None, seed=None):
    binary_var_dict = {}
    marked_c_nodes = []
    c_nodes = []
    r_nodes = []
    e_nodes = []
    node_index = 1

    graph.del_node("\"\\r\\n\"")  # Muss gelöscht werden, da als Artefakt immer zusätzlich als Knoten generiert wird
    graph.del_node("\"\\n\"")

    try:  # Check ob die Datei auch nur mit Knoten eines bestimmten Formats ausgestattet ist.
        for node in graph.get_nodes():
            if len(node.get_attributes()) < 1 or len(node.get_attributes()) > 2:
                raise ValueError(f"The Node \"{node.get_name()}\" has 0 or more than 3 Attributes! ")

            node_type = node.get('type')
            color = node.get_color()

            # Error if the node has no type
            if not node_type:
                raise ValueError(f"The Node \"{node.get_name()}\" has no type attribute")

            if node_type not in ["E", "C", "R"]:
                raise ValueError(f"The Node \"{node.get_name()}\" has invalid type {node_type}")

            # Error if a node has color attribute, and it's not type C
            if color:
                if node_type != "C":
                    raise ValueError(f"The Node \"{node.get_name()}\" has color attribute but is not type C")

            if color == "red" and node_type == "C":
                marked_c_nodes.append(node.get_name())

            if node_type == "E":
                e_nodes.append(node.get_name())
            if node_type == "C":
                c_nodes.append(node.get_name())
            if node_type == "R":
                r_nodes.append(node.get_name())

            binary_var_dict[node.get_name()] = f'x_{node_index}'
            node_index += 1

        if marked_c_nodes:
            tie_graph_problem = TieGraphProblem(binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes)
            return tie_graph_problem

        if k:
            if seed is None:
                # Use the current time to generate a random seed
                seed = int(time.time())

            # Set the seed for the random number generator for reproducibility
            random.seed(seed)

            # Select k random items from the input list
            marked_c_nodes = random.sample(c_nodes, k)

            # Return the selected items and the seed used
            tie_graph_problem = TieGraphProblem(binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes, seed)
            return tie_graph_problem

        # Check Auto Case, where if no target compound is given, it will automatically look, which on ist there from
        # the QUTIE Paper Disease List
        disease_comp_set = {"C00036", "C01165", "C00183", "C00148", "C00074", "C00082", "C05382", "C00169", "C00049",
                            "C00236", "C03287", "C00327", "C01005",
                            "C00199", "C01005" "C00141", "C00407", "C00233", "C00119"}

        marked_c_nodes = list(set(c_nodes).intersection(disease_comp_set))
        if not marked_c_nodes:
            raise ValueError(
                "The given file has no marked Nodes or no node with the name listed in the OG Paper Disease "
                "list")

        tie_graph_problem = TieGraphProblem(binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes)
        return tie_graph_problem


    except Exception as e:
        print(f'Error creating binary variables: {str(e)}')
        return None


def generate_equation_auto_penalty_and_enzyme_damage(graph, dict_and_sorted_nodes):
    try:
        binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes = dict_and_sorted_nodes.get_dict_and_lists()
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

        # Damage Term has penalty factor 2 (not explicitly written out)
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
                # check for circular subgraphs (reaction nodes connected with 2 edges in both directions to the same
                # compound.
                predecessors_set = set(predecessors(graph, marked_node))
                successors_set = set(successors(graph, marked_node))
                if predecessors_set & successors_set:
                    double_edged_reactions.update(predecessors_set & successors_set)
                    circular_marked_c_nodes.append(marked_node)

            if double_edged_reactions:
                for r_node in double_edged_reactions:  # For all Reaction nodes
                    r_nodes.remove(r_node)
                    equation += f'(1 - {binary_var_dict[r_node]}) + '

                    non_marked_reaction_predecessors = []

                    for e_node in predecessors(graph, r_node):
                        curr_node = graph.get_node(e_node)  # curr_node here is a Node List Object, not a node_Name
                        # nor a Node!
                        if not curr_node:
                            curr_node = graph.get_node(f"\"{e_node}\"")  # curr_node here is a Node List Object, not
                            # a node_Name nor a Node!
                        if curr_node[0].get("type") == "E":
                            non_marked_reaction_predecessors.append(e_node)

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
                if len(list(predecessors(graph, r_node))) > 0:  # If the number of predecessors is bigger than zero.
                    currBinVar_r = binary_var_dict[r_node]  # currBinVar contains r_node
                    if len(list(predecessors(graph, r_node))) == 1:  # If r_node has just one predecessor
                        equation += f'{binary_var_dict[list(predecessors(graph, r_node))[0]]} + {currBinVar_r}'
                        equation += f' - 2 * {binary_var_dict[list(predecessors(graph, r_node))[0]]} * {currBinVar_r} + '
                    else:  # If there are more Predecessors than 1:
                        equation += f'(1 - {currBinVar_r} - y_r{extra_var_index_r}'
                        equation += f' + 2 * {currBinVar_r} * y_r{extra_var_index_r}) + '
                        binary_var_dict[f'extra Var for r{extra_var_index_r}'] = f'y_r{extra_var_index_r}'
                        react_pena_term += reaction_reduction_penalty(list(predecessors(graph, r_node)),
                                                                      binary_var_dict,
                                                                      extra_var_index_r)
                        extra_var_index_r += 1

        # Case 4 Compound and Reduction Penalty
        for c_node in c_nodes:
            if len(list(predecessors(graph, c_node))) != 0:
                currBinVar_c = binary_var_dict[c_node]
                if len(list(predecessors(graph, c_node))) == 1:
                    fill4Var = binary_var_dict[list(predecessors(graph, c_node))[0]]
                    equation += f'({currBinVar_c} + {fill4Var} - 2 * {currBinVar_c} * {fill4Var}) + '
                else:
                    equation += f'({currBinVar_c} + y_c{extra_var_index_c} - 2 * {currBinVar_c} * y_c{extra_var_index_c}) + '
                    binary_var_dict[f'extra Var for c{extra_var_index_c}'] = f'y_c{extra_var_index_c}'
                    compo_pena_term += compound_reduction_penalty(list(predecessors(graph, c_node)), binary_var_dict,
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


# Generate qubo term from graph, but generates two Objective functions, to lower the complexity. -> AKTUELL!
def generating_qubo_term_from_graph_two_part(file_path):
    graph = read_dot_file_pydotplus(file_path)
    if not graph:
        print("The Graph is weirdly, None...")
        return None
    dict_and_sorted_nodes = create_binary_variables_plus_enzymes(graph, 3)  # TODO WICHTIG, ich habe hier  k=3
    # gesetzt, weil der Polyketide Graph keine diseased Knoten hat!!
    if dict_and_sorted_nodes:
        if dict_and_sorted_nodes.has_seed():
            seed = dict_and_sorted_nodes.get_seed()
        # binary_vars, output_term = generate_equation_auto_penalty(graph, dict_and_sorted_nodes[0],
        binary_vars, output_term, output_damage_only = generate_equation_auto_penalty_and_enzyme_damage(graph,

                                                                                                        dict_and_sorted_nodes)
        return binary_vars, output_term, output_damage_only, graph, dict_and_sorted_nodes[1]  # Added the
        # information of the target nodes, to color the result node, if no node was marked from the start or the
        # Target(s) was/were randomly chosen


# Example usage
# file_path = '/workspaces/graph-coloring/Citrate_cycle_marked.dot'
# graph = read_dot_file_pydotplus(file_path)
# binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes = create_binary_variables_plus_enzymes(graph)

if __name__ == '__main__':
    graph = read_dot_file_pydotplus(
        "C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\graphoz\\c.dot")
    # print(graph.get_node("H"))
    nodii = "1.013.2"
    print(graph.get_node(f"\"{nodii}\""))
    # print(generating_qubo_term_from_graph_two_part("/workspaces/graph-coloring/Citrate_cycle_marked.dot")[2])

# print(graph.get_node('"\r\n"'))
# print(node.get_name())
# print(node.get('type'))
# print(node.get_attributes())
# print(node.get_color())

# my_function("This is required", optional_arg1="This is optional 1")  # Providing the first optional argument
# my_function("This is required", optional_arg2="This is optional 2")  # Providing the second optional argument
# my_function("This is required", "This is optional 1", "This is optional 2")  # Providing both optional arguments
