import pydot


# Function that, if you give it lines, will return a .dot graph.
# DOESN'T READ IN A FILE!! ONLY THE FUNCTION TO READ LINES AND CREATE .DOT GRAPH
def input2dotfile(lines):
    # Open the text file and read the data

    # Create a new directed graph
    graph = pydot.Dot(graph_type='digraph')

    # Parse the data and add it to the graph
    line_index = 0

    # Add enzyme nodes
    n = int(lines[line_index].strip())
    line_index += 1
    for _ in range(n):
        node_name = lines[line_index].strip()
        line_index += 1
        graph.add_node(pydot.Node(node_name, type='E'))

    # Add reaction nodes
    k = int(lines[line_index].strip())
    line_index += 1
    for _ in range(k):
        node_name = lines[line_index].strip()
        line_index += 1
        graph.add_node(pydot.Node(node_name, type='R'))

    # Add compound nodes
    i = int(lines[line_index].strip())
    line_index += 1
    for _ in range(i):
        node_name = lines[line_index].strip()
        line_index += 1
        graph.add_node(pydot.Node(node_name, type='C'))

    # Add edges
    j = int(lines[line_index].strip())
    line_index += 1
    for _ in range(j):
        node_names = lines[line_index].strip().split()
        line_index += 1
        graph.add_edge(pydot.Edge(node_names[0], node_names[1]))

    return graph

