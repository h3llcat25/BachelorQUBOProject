import os
from contextlib import redirect_stdout


def run_filo_reader(directory):
    # Iterate through each subdirectory in the given directory
    for subdir, dirs, files in os.walk(directory):
        for dirct in dirs:
            # Define the path to the current subdirectory
            subdirectory_path = os.path.join(subdir, dirct)

            # Create a corresponding "_solutions" directory
            output_directory_path = os.path.join(subdir, f"{dirct}_solutions")
            os.makedirs(output_directory_path, exist_ok=True)

            # Iterate through each file in the subdirectory
            for file in os.listdir(subdirectory_path):
                file_path = os.path.join(subdirectory_path, file)

                # Define the output file path
                output_file_path = os.path.join(output_directory_path, f"{file}_output.txt")

                # Redirect stdout to the output file
                with open(output_file_path, 'w') as output_file, redirect_stdout(output_file):
                    print(file)
                    # psa.param_search_and_auto_solver_auto_penalty_gurobipy_model_enzymes(file_path)

run_filo_reader("C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\largeDots")

'''import pydotplus
    import random
    import copy

    file_path = ("C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\largeDots_w_marked_and_tests\\eco_filtering_dot\\Biosynthesis of amino acids_m.dot")


    # Load the .dot file
    graph = pydotplus.graph_from_dot_file(file_path)

    # Copy the list of edges to avoid modifying the list while iterating
    edges = copy.deepcopy(graph.get_edges())
    # Delete every third edge
    for i in range(2, len(edges), 3):
        graph.del_edge(edges[i].get_source(), edges[i].get_destination())

    # Get all node names
    node_names = [node.get_name() for node in graph.get_nodes()]

    # Select 3 random nodes, ensuring we have enough nodes
    if len(node_names) >= 3:
        selected_nodes = random.sample(node_names, 3)
    else:
        print("The graph does not have enough nodes.")
        selected_nodes = node_names

    # For each selected node, print out its specific predecessors and successors
    for node_name in selected_nodes:
        predecessors = [edge.get_source() for edge in graph.get_edges() if edge.get_destination() == node_name]
        successors = [edge.get_destination() for edge in graph.get_edges() if edge.get_source() == node_name]
        print(f"Node: {node_name}")
        print(f"Predecessors: {predecessors}")
        print(f"Successors: {successors}")
        print()

    # Write the modified graph to a new .dot file
    graph.write('ZHELLO_graph.dot')
'''
