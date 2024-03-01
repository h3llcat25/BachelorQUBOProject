import param_and_solve_automation as psa
import os
import sys
import pydotplus


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


def write_first_lines_to_output(start_directory, output_file='matchings.txt'):
    disease_comp_list = ["C00036", "C01165", "C00183", "C00148", "C00074", "C00082", "C05382", "C00169", "C00049",
                         "C00236", "C03287", "C00327", "C01005",
                         "C00199", "C01005" "C00141", "C00407", "C00233", "C00119"]

    with open(output_file, 'w') as outfile:  # Open the output file for writing
        for subdir, dirs, files in os.walk(start_directory):  # Traverse the directory
            for file in files:
                file_path = os.path.join(subdir, file)  # Get the full path of the file
                outfile.write(f"File: {file_path}\nMatching Nodes: {search_node_names(file_path, disease_comp_list)}\n\n")


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


# Example usage
# run_filo_reader("largeDots")
disease_comp_list = ["C00036", "C01165", "C00183", "C00148", "C00074", "C00082", "C05382", "C00169", "C00049",
                     "C00236", "C03287", "C00327", "C00199", "C01005", "C00141", "C00407", "C00233", "C00119"]

# print(disease_comp_list)
# print(set(disease_comp_list))

# Example usage
start_directory = 'C:\\Users\\marsh\\Documents\\Python Bachelor\\QUBO_Project_BA\\graphStuff\\smallDots'  # Replace with your start directory
write_first_lines_to_output(start_directory)
