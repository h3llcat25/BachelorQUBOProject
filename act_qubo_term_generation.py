import math
import pydotplus
import random
import time

from tie_graph_problem import TieGraphProblem
from tie_qubo_struct import TieQuboStruct

# Helpfunction for getting predecessors of a specific Node in the .Dot file, Wichtig, gibt ne list of strings zurück
# Ist Smart, weil wenn es durch die
def smart_predecessors(graph, target_node_name):
    target_node = graph.get_node(target_node_name)[0]
    trgt_node_type = target_node.get('type')

    if trgt_node_type == "C":
        predecessor_set = {edge.get_source() for edge in graph.get_edges() if edge.get_destination() == target_node_name}
        successor_set = {edge.get_destination() for edge in graph.get_edges() if edge.get_source() == target_node_name}

        for r_node_name in (predecessor_set & successor_set):
            r_node = graph.get_node(r_node_name)[0]
            if r_node.get('color') is None:
                r_node.set('color', 'lime')

        return list(predecessor_set)

    if trgt_node_type == "R":
        predecessor_set = {edge.get_source() for edge in graph.get_edges() if edge.get_destination() == target_node_name}
        if target_node.get('color') is None:
            return list(predecessor_set)

        # If it is lime, it means, it is a revertable Reaktion and it needs different treatment
        if target_node.get('color') == "lime":
            successor_set = {edge.get_destination() for edge in graph.get_edges() if edge.get_source() == target_node_name}
            non_rev_predecessors = predecessor_set - successor_set
            return list(non_rev_predecessors)

        else:
            raise ValueError("sometin wrong with predecessors function")


# Helpfunction for getting successors of a specific Node in the .Dot file
def successors(graph, target_node_name):
    # Specify the node for which you want to find the successors

    # Iterate over the edges in the graph, Check if the edge comes from the target node, If so, add the get_destination node of the edge to the successors list
    successors = [edge.get_destination() for edge in graph.get_edges() if edge.get_source() == target_node_name]

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

    abstrct_var_dict = {}
    r_var_count = 1
    e_var_count = 1


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

            if node_type == "C":
                if color == "red":
                    marked_c_nodes.append(node.get_name())
                else:
                    c_nodes.append(node.get_name())
            if node_type == "E":
                e_nodes.append(node.get_name())
                abstrct_var_dict[node.get_name()] = f'E{e_var_count}'
                e_var_count += 1
            if node_type == "R":
                r_nodes.append(node.get_name())
                abstrct_var_dict[node.get_name()] = f'R{r_var_count}'
                r_var_count += 1

            binary_var_dict[node.get_name()] = f'x_{node_index}'
            node_index += 1

        tie_graph_problem = TieGraphProblem(bin_var_dict=binary_var_dict, list3=e_nodes, list4=r_nodes)
        tie_graph_problem.set_abstrct_dict(abstrct_var_dict)

        if marked_c_nodes: # 1. Case: There are marked C nodes in the dot file, so we don't have to generate them
            abstrct_var_dict.update({node: f'C{i + 1}' for i, node in enumerate(c_nodes)})
            abstrct_var_dict.update({node: f'TC{i + 1}' for i, node in enumerate(marked_c_nodes)})
            tie_graph_problem.set_c_nodes(c_nodes)
            tie_graph_problem.set_marked_c_nodes(marked_c_nodes)
            return tie_graph_problem

        if k: # 2. Case to randomly take marked c_nodes
            if seed is None:
                # Use the current time to generate a random seed
                seed = int(time.time())

            # Set the seed for the random number generator for reproducibility
            random.seed(seed)

            # Select k random items from the input list
            marked_c_nodes = random.sample(c_nodes, k)

            # And take them out of the c_nodes list
            c_nodes_set = set(c_nodes)
            marked_c_nodes_set = set(marked_c_nodes)
            # Use set difference to remove marked_c_nodes from c_nodes
            c_nodes_set -= marked_c_nodes_set
            c_nodes = list(c_nodes_set)

            # Return the selected items and the seed used
            abstrct_var_dict.update({node: f'C{i + 1}' for i, node in enumerate(c_nodes)})
            abstrct_var_dict.update({node: f'TC{i + 1}' for i, node in enumerate(marked_c_nodes)})

            tie_graph_problem.set_abstrct_dict(abstrct_var_dict)
            tie_graph_problem.set_marked_c_nodes(marked_c_nodes)
            tie_graph_problem.set_c_nodes(c_nodes)
            tie_graph_problem.set_seed(seed)
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

        abstrct_var_dict.update({node: f'C{i + 1}' for i, node in enumerate(c_nodes)})
        abstrct_var_dict.update({node: f'TC{i + 1}' for i, node in enumerate(marked_c_nodes)})

        tie_graph_problem.set_abstrct_dict(abstrct_var_dict)
        tie_graph_problem.set_marked_c_nodes(marked_c_nodes)
        tie_graph_problem.set_c_nodes(c_nodes)
        return tie_graph_problem


    except Exception as e:
        print(f'Error creating binary variables: {str(e)}')
        return None


# AKTUELL
def generate_equation(graph, dict_and_sorted_nodes, cmp_dmg_weight, enzym_dmg_weight):
    try:
        binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes = dict_and_sorted_nodes.get_dict_and_lists()
        abstract_var_dict = dict_and_sorted_nodes.get_abstrct_dict()
        abstrct_terms = ''

        # start from where to count for the extra variables for reactions r and compounds c.
        # Setting up Abstract Dictionary to better check in the results
        extra_var_index_r = 1
        extra_var_index_c = 1

        damage_only_equation = ''
        equation = ''
        react_pena_term = ''
        compo_pena_term = ''
        # Start of the set-up for inevitable inhibitions:
        # 1. Marked compound nodes and reactions
        all_successors_set_up_reacts = set()
        all_predecessors_set_up_reacts = set()

################ Start Set up ##########

        # 1.1 Generate the terms for each marked_node x, that it has to get inhibited (1-x)
        # And create the predecessor and successor sets for the next step
        for marked_node in marked_c_nodes:
            equation += f'(1 - {binary_var_dict[marked_node]}) + '

            marked_node_successors = set(successors(graph, marked_node))
            marked_node_predecessors = set(smart_predecessors(graph, marked_node))

            all_successors_set_up_reacts.update(marked_node_successors)
            all_predecessors_set_up_reacts.update(marked_node_predecessors)

        # 2. Generate the terms for each predecessor and successor node y, that it has to get inhibited (1-y)
        all_set_up_reacts = all_predecessors_set_up_reacts | all_successors_set_up_reacts
        equation += ''.join([f'(1 - {binary_var_dict[node]}) + ' for node in all_set_up_reacts])

        # 3. take out every r_node that is only a successor and not a reversable or a predecessor of a target
        for r_node in (all_successors_set_up_reacts - all_predecessors_set_up_reacts):
            r_nodes.remove(r_node)

        ####### Ende Set UP ##############

        # Automatic Penalty calculation, based on if all Reactions, Enzymes and marked Compounds of the metabolic
        # network are inhibited
        meta_rulebreak_penalty = int(math.floor(0.5 * len(binary_var_dict)))
        print(f'Meta Rulebreaker Penalty is: {meta_rulebreak_penalty}')

        equation = f'{meta_rulebreak_penalty}*({equation}'

        # Case 2 Compound and Reduction Penalty
        for c_node in c_nodes:
            c_predecessors = smart_predecessors(graph, c_node)
            c_pred_len = len(c_predecessors) # length of the predecessor list of c
            if c_pred_len != 0:
                currBinVar_c = binary_var_dict[c_node]
                if c_pred_len == 1:
                    fill4Var = binary_var_dict[c_predecessors[0]]
                    equation += f'({currBinVar_c} + {fill4Var} - 2 * {currBinVar_c} * {fill4Var}) + '
                else:
                    equation += f'({currBinVar_c} + y_c{extra_var_index_c} - 2 * {currBinVar_c} * y_c{extra_var_index_c}) + '
                    binary_var_dict[f'extra Var for c{extra_var_index_c}'] = f'y_c{extra_var_index_c}'
                    compo_pena_term += compound_reduction_penalty(c_predecessors, binary_var_dict, extra_var_index_c)
                    extra_var_index_c += 1

        # 3. Reactions and Reduction Penalty:
        # With smart predecessors, the list will automatically handle revertable reactions
        if len(r_nodes) > 0:
            for r_node in r_nodes:  # For all Reaction nodes
                r_smrt_prdcssrs = smart_predecessors(graph, r_node)
                r_smrt_pred_len = len(r_smrt_prdcssrs)
                if r_smrt_pred_len > 0:  # If the number of predecessors is bigger than zero.
                    currBinVar_r = binary_var_dict[r_node]  # currBinVar contains r_node
                    if r_smrt_pred_len == 1:  # If r_node has just one predecessor
                        equation += f'{binary_var_dict[r_smrt_prdcssrs[0]]} + {currBinVar_r}'
                        equation += f' - 2 * {binary_var_dict[r_smrt_prdcssrs[0]]} * {currBinVar_r} + '
                    else:  # If there are more Predecessors than 1:
                        equation += f'(1 - {currBinVar_r} - y_r{extra_var_index_r}'
                        equation += f' + 2 * {currBinVar_r} * y_r{extra_var_index_r}) + '
                        binary_var_dict[f'extra Var for r{extra_var_index_r}'] = f'y_r{extra_var_index_r}'
                        react_pena_term += reaction_reduction_penalty(r_smrt_prdcssrs, binary_var_dict, extra_var_index_r)
                        extra_var_index_r += 1

        # FIND ME
        # Case 3 Damage Term has different penalty weights: cmp_dmg_weight ++++++++++++++++++++++++++++++++++++
        damage_only_equation += f'{cmp_dmg_weight}* ('
        for c_node in c_nodes:
            if c_node not in marked_c_nodes:
                damage_only_equation += f'{binary_var_dict[c_node]} + '
        damage_only_equation = damage_only_equation.rstrip(' + ')

        # Enzyme Damage Term has penalty factor enzym_dmg_weight
        damage_only_equation += f") + {enzym_dmg_weight}*("
        for e_node in e_nodes:
            damage_only_equation += f'{binary_var_dict[e_node]} + '
        damage_only_equation = damage_only_equation.rstrip(' + ')
        #damage_only_equation += ")"

        damage_only_equation += f") + 3*("
        for r_node in r_nodes:
            damage_only_equation += f'{binary_var_dict[r_node]} + '
        damage_only_equation = damage_only_equation.rstrip(' + ')
        damage_only_equation += ")"

        if react_pena_term:
            react_pena_term = react_pena_term.rstrip(' + ')
            equation += f'({react_pena_term}) + '
        if compo_pena_term:
            compo_pena_term = compo_pena_term.rstrip(' + ')
            equation += f'({compo_pena_term}) + '
        equation = equation.rstrip(' + ')
        equation += ")"
        # print(f'Generated equation: {equation}')
        tie_qubo_struct = TieQuboStruct(binary_var_dict, equation, damage_only_equation, graph,
                                        dict_and_sorted_nodes.get_marked_c_nodes())
        return tie_qubo_struct

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
                pen_r_equation += f' - 2 * (1 - {secVal}) * y_r{extra_var_index_r}_{i} + 3 * y_r{extra_var_index_r}_{i}) + '
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
def generating_qubo_term_from_graph_two_part(file_path, cmp_dmg_weight=5, enzym_dmg_weight=2):
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
        tie_qubo_struct = generate_equation(graph, dict_and_sorted_nodes, cmp_dmg_weight=cmp_dmg_weight, enzym_dmg_weight=enzym_dmg_weight)
        return tie_qubo_struct  # Added the
        # information of the target nodes, to color the result node, if no node was marked from the start or the
        # Target(s) was/were randomly chosen
    else:
        print("there is somehow not a dict_and_sorted_nodes")

# Example usage
# file_path = '/workspaces/graph-coloring/Citrate_cycle_marked.dot'
# graph = read_dot_file_pydotplus(file_path)
# binary_var_dict, marked_c_nodes, c_nodes, r_nodes, e_nodes = create_binary_variables_plus_enzymes(graph)

if __name__ == '__main__':

    file_path = ("C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\largeDots_w_marked_and_tests\\eco_filtering_dot\\Biosynthesis of amino acids_m.dot")


    # Load the .dot file
    graph = pydotplus.graph_from_dot_file(file_path)

    # Write the modified graph to a new .dot file
    graph.write('ZHELLO_graph.dot')

