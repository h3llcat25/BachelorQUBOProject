#9.4.24
# not aktuell
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


        abstrct_dict = dict_and_sorted_nodes.get_abstrct_dict()
        abstrct_terms = ''

        double_edged_reactions = set()
        circular_marked_c_nodes = []

        # Automatic Penalty calculation, based on if all Reactions, Enzymes and marked Compounds of the metabolic
        # network are inhibited + 10%
        meta_rulebreak_penalty = 2*int(math.floor(0.5*len(binary_var_dict)))
        print(f'Meta Rulebreaker Penalty is: {meta_rulebreak_penalty}')

        # Damage Term has penalty factor 2 (not explicitly written out)
        damage_only_equation += '8 * ('
        for c_node in c_nodes:
            if c_node not in marked_c_nodes:
                damage_only_equation += f'{binary_var_dict[c_node]} + '

        damage_only_equation = damage_only_equation.rstrip(' + ')
        damage_only_equation += ") + 3*("

        # Enzyme Damage Term has penalty factor 3 (not explicitly written out)
        for e_node in e_nodes:
            damage_only_equation += f'{binary_var_dict[e_node]} + '
        damage_only_equation = damage_only_equation.rstrip(' + ')
        damage_only_equation += ") + 5*("

        # Reaction Damage Term has penalty factor 1 (not explicitly written out)
        for r_node in r_nodes:
            damage_only_equation += f'{binary_var_dict[r_node]} + '
        damage_only_equation = damage_only_equation.rstrip(' + ')
        damage_only_equation += ")"

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
        tie_qubo_struct = TieQuboStruct(binary_var_dict, equation, damage_only_equation, graph, dict_and_sorted_nodes.get_marked_c_nodes())
        return tie_qubo_struct

    except Exception as e:
        print(f'Error generating equation: {str(e)}')
