import param_and_solve_automation as psa
import os
import sys
import pydotplus
import pandas as pd


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


def run_filo_reader(directory):
    disease_comp_list = ["C00036", "C01165", "C00183", "C00148", "C00074", "C00082", "C05382", "C00169", "C00049",
                         "C00236", "C03287", "C00327", "C01005",
                         "C00199", "C01005" "C00141", "C00407", "C00233", "C00119"]

    # Iterate through each subdirectory in the given directory
    for subdir, dirs, files in os.walk(directory):
        for dirct in dirs:
            # Define the path to the current subdirectory
            subdirectory_path = os.path.join(subdir, dirct)

            # Create a corresponding "texts_" directory
            output_directory_path = os.path.join(subdir, f"texts_{dirct}")
            os.makedirs(output_directory_path, exist_ok=True)

            # Iterate through each file in the subdirectory
            for file in os.listdir(subdirectory_path):
                file_path = os.path.join(subdirectory_path, file)

                # Define the output file path
                output_file_path = os.path.join(output_directory_path, f"{file}_output.txt")

                # Redirect stdout to the output file
                with open(output_file_path, 'w') as output_file, redirect_stdout(output_file):
                    psa.param_search_and_auto_solver_auto_penalty_gurobipy_model_enzymes(file_path)


def write_matching_nodes_from_files(start_directory, output_file='matchings.txt'):
    disease_comp_list = ["C00036", "C01165", "C00183", "C00148", "C00074", "C00082", "C05382", "C00169", "C00049",
                         "C00236", "C03287", "C00327", "C01005",
                         "C00199", "C01005" "C00141", "C00407", "C00233", "C00119"]

    with open(output_file, 'w') as outfile:  # Open the output file for writing
        for subdir, dirs, files in os.walk(start_directory):  # Traverse the directory
            for file in files:
                file_path = os.path.join(subdir, file)  # Get the full path of the file
                outfile.write(f"File: {file_path}\nMatching Nodes: {search_node_names(file_path, disease_comp_list)}\n\n")


def count_difference_of_bin_vars_for_each_method(start_directory):
    # Define the DataFrame columns
    columns = ["File Name", "Nr. of Nodes", "Nr. of total Bin Vars (Qutie)", "Nr. of total Bin Vars (new Method)",
               "Aux. Variable (Qutie)", "Aux. Variables (new Method)", "% difference of total Bin Vars",
               "% difference Nr. of aux. Variables"]
    data = []

    # Walk through the directory and its subdirectories
    for subdir, dirs, files in os.walk(start_directory):  # Traverse the directory
        for file in files:
            file_path = os.path.join(subdir, file)  # Get the full path of the fil

            relative_subdir = os.path.relpath(subdir, start_directory)
            # Join the relative subdirectory path with the file name
            rel_file_path = os.path.join(relative_subdir, file)

            # Apply the get_tuples function to each file
            tuple_result = calculate_total_bin_vars_needed(file_path)
            # Calculate the percentage differences
            percent_diff_bin_vars = ((tuple_result[2] - tuple_result[1]) / tuple_result[1]) * 100 if tuple_result[
                1] else 0
            percent_diff_aux_vars = ((tuple_result[4] - tuple_result[3]) / tuple_result[3]) * 100 if tuple_result[
                3] else 0
            # Append the results including the file name and the calculated percentages to the data list
            data.append([rel_file_path] + list(tuple_result) + [round(percent_diff_bin_vars, 2), round(percent_diff_aux_vars, 2)])

        # Create a DataFrame from the data
        df = pd.DataFrame(data, columns=columns)

        # Export the DataFrame to an Excel file
        df.to_excel("large_dots_stats.xlsx", index=False)


# Context manager to redirect stdout
class redirect_stdout:
    def __init__(self, new_target):
        self.new_target = new_target

    def __enter__(self):
        self.old_target = sys.stdout
        sys.stdout = self.new_target
        return self.new_target

    def __exit__(self, exc_type, exc_val, exc_tb):
        sys.stdout = self.old_target


def search_node_names(dot_file_path, names_to_search):
    # Load the .dot file and create a graph object
    graph = pydotplus.graph_from_dot_file(dot_file_path)

    # Get all node names from the graph and remove quotes
    node_names = {node.get_name().strip('"') for node in graph.get_nodes()}

    # Convert the list of names to search into a set for faster lookup
    search_set = set(names_to_search)

    # Find the intersection of node names and names to search
    matching_names = list(node_names.intersection(search_set))

    return matching_names


def calculate_total_bin_vars_needed(dot_file_path):
    graph = pydotplus.graph_from_dot_file(dot_file_path)
    graph.del_node("\"\\r\\n\"")  # Muss gelöscht werden, da als Artefakt immer zusätzlich als Knoten generiert wird

    nr_of_nodes = 0
    # nodes_with_type_E = 0
    sum_predecessors_minus_one = 0
    sum_predecessors_plus_one = 0

    # Iterate over all nodes in the graph
    for node in graph.get_nodes():
        nr_of_nodes += 1
        # Get the list of predecessors for the current node
        predecessor_lst = predecessors(graph, node.get_name())

        # Check if the node has at least one predecessor
        if predecessor_lst:
            # attributes = node.get_attributes()
            # if attributes.get('type') == '"E"':
            #     nodes_with_type_E += 1
            num_predecessors = len(predecessor_lst)
            sum_predecessors_minus_one += (num_predecessors - 1)
            sum_predecessors_plus_one += (num_predecessors + 1)
    return  nr_of_nodes, nr_of_nodes + sum_predecessors_plus_one, nr_of_nodes + sum_predecessors_minus_one, sum_predecessors_plus_one, sum_predecessors_minus_one

# Example usage
# run_filo_reader("largeDots")
disease_comp_list = ["C00036", "C01165", "C00183", "C00148", "C00074", "C00082", "C05382", "C00169", "C00049",
                     "C00236", "C03287", "C00327", "C00199", "C01005", "C00141", "C00407", "C00233", "C00119"]

# print(disease_comp_list)
# print(set(disease_comp_list))

# Example usage
start_directory = 'C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\largeDots'  # Replace with your start directory
count_difference_of_bin_vars_for_each_method(start_directory)
