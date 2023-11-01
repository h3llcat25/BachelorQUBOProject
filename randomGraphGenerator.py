import networkx as nx
import random as rnd
import math


class React:

    def __init__(self, id_nr):
        self.id_nr = id_nr
        self.input_C = []
        self.output_C = []
        self.enzyme_list = []

    def write_edges(self, graph):
        for e in self.enzyme_list:
            graph.add_edge(f"E{e}", f"R{self.id_nr}")
        for c_in in self.input_C:
            graph.add_edge(f"C{c_in}", f"R{self.id_nr}")
        for c_out in self.output_C:
            graph.add_edge(f"R{self.id_nr}", f"C{c_out}")


def generate_graph(enzy_nodes_nr, reac_nodes_nr, comp_nodes_nr, red_comp_nodes_nr=None):
    if enzy_nodes_nr < 1:
        return "There must be at least 1 enzyme node"
    if reac_nodes_nr < 1:
        return "There must be at least 1 reaction node"
    if comp_nodes_nr < 2:
        return "There must be at least 2 compound nodes"

    if red_comp_nodes_nr is None:
        red_comp_nodes_nr = math.ceil(rnd.uniform(1, comp_nodes_nr / 4))
    elif red_comp_nodes_nr < 1 or red_comp_nodes_nr > comp_nodes_nr:
        return "Red compound nodes must be between 1 and the total number of compound nodes"

    G = nx.MultiDiGraph()
    # Determine randomly, which Compound nodes to use as target Compounds
    targetList = rnd.sample(range(comp_nodes_nr), red_comp_nodes_nr)

    # Add Enzyme, Reaction and Compound Nodes to the Graph G
    for i in range(enzy_nodes_nr):
        G.add_node(f"E{i}", type="E")
    for i in range(reac_nodes_nr):
        G.add_node(f"R{i}", type="R")
    # The randomly assigned Target Compounds get marked with the color red, while the not-targeted ones get no color
    # attribute
    for i in range(comp_nodes_nr):
        if i in targetList:
            G.add_node(f"C{i}", type="C", color="red")
        else:
            G.add_node(f"C{i}", type="C")

    # Create random Edges following the rules:
    # 1. A reaction has at least one Enzyme and one Compound connected with an incoming edge and one different Compound
    #   connected with an outgoing edge
    # 2. Enzymes and Compounds can be part and output of multiple reactions.
    all_reacts_list = []
    all_enzymes_list = list(range(enzy_nodes_nr))
    all_compounds_list = list(range(comp_nodes_nr))
    choosable_comp_list = list(range(comp_nodes_nr))
    all_used_enzymes_set = set()
    all_used_compounds_set = set()

    # Create all the Reaction Objects as sceleton for the network and save them in "all_reacts_list".
    all_reacts_list = [React(i) for i in range(reac_nodes_nr)]

    # Give every reaction object a random Enzyme(-id), a random "Input Compound" and a different random "Output
    # Compound". Save random enzyme in a set for later, and take the chosen Input Compound out, for choosing a random
    # Output Compound
    for r in all_reacts_list:
        temp_enzy = rnd.choice(all_enzymes_list)
        all_used_enzymes_set.add(temp_enzy)
        r.enzyme_list.append(temp_enzy)
        # Input Compound Part:
        tmp_input_comp = rnd.choice(choosable_comp_list)  # choose a random Compound of all Compounds and store it
        r.input_C.append(tmp_input_comp)  # put it as input in the Reaction object
        all_used_compounds_set.add(tmp_input_comp)  # Store what Compounds were already used at least once
        choosable_comp_list.remove(tmp_input_comp)  # For all coming React Objects, take that Compound out the pool
        # Output compound Part:
        all_comps_copy = all_compounds_list[:]  # make a copy, to choose a random Compound as "output"
        all_comps_copy.remove(tmp_input_comp)  # remove the input compound to not have it as input compound
        temp_comp = rnd.choice(all_comps_copy)  # chose a random comp as output from remaining compounds
        all_used_compounds_set.add(temp_comp)  # put it in the set of already used compounds
        r.output_C.append(temp_comp)  # Put it in the output list of the R object.

    if len(all_enzymes_list) != len(all_used_enzymes_set):
        remain_enzymes = [x for x in all_enzymes_list if x not in all_used_enzymes_set]
        for renz in remain_enzymes:
            rnd.choice(all_reacts_list).enzyme_list.append(renz)

    if len(all_compounds_list) != len(all_used_compounds_set):
        remain_compounds = [x for x in all_compounds_list if x not in all_used_compounds_set]
        for recomp in remain_compounds:
            tmp_random = rnd.choice([0, 1])
            if tmp_random == 0:
                rnd.choice(all_reacts_list).input_C.append(recomp)
            else:
                rnd.choice(all_reacts_list).output_C.append(recomp)

    for react in all_reacts_list:
        react.write_edges(G)

    # Convert to PyGraphviz graph, write to .dot file, and create a PNG image
    A = nx.nx_agraph.to_agraph(G)
    A.write('graph.dot')


generate_graph(3, 3, 3)
