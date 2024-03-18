import pydotplus


def check_node_attribute(dot_file_path, node_name):
    # Load the .dot graph
    graph = pydotplus.graph_from_dot_file(dot_file_path)

    # Iterate over the nodes in the graph
    for node in graph.get_nodes():
        # Node names are quoted strings in pydotplus, so we strip the quotes
        if node.get_name().strip('"') == node_name:
            # Check if the node has the attribute 'type' with value 'E'
            if node.get("type") == "E":
                print(f"The node '{node_name}' has the type 'E'.")
                return True
            else:
                print(f"The node '{node_name}' does not have the type 'E'.")
                return False

    # If the node was not found in the graph
    print(f"The node '{node_name}' was not found in the graph.")
    return False

# Example usage
dot_file_path = 'C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\graphStuff\\graphoz\\a.dot'  # Path to your .dot file
node_name = 'YourNodeName'  # The name of the node you want to check

print(check_node_attribute(dot_file_path, "5.1.3.13"))
