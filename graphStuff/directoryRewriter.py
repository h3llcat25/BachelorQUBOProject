import os
import inputToDotFile

# Script to change all .input graph files of a directory to .dot Graph files
# Define the source directory and the destination directory
src_dir = "paperData/large_data/mmu_filtering"
dst_dir = "largeDots/mmu_filtering_dot"

# List all the .txt files in the source directory
input_files = [f for f in os.listdir(src_dir) if f.endswith('.input')]

# Create the destination directory if it doesn't exist
if not os.path.exists(dst_dir):
    os.makedirs(dst_dir)

# For each .input file in the source directory, count the number of characters
# and create a new .txt file in the destination directory that contains this count
for file in input_files:
    # Open the file in read mode
    with open(os.path.join(src_dir, file), 'r') as f:
        lines = f.readlines()

    graph = inputToDotFile.input2dotfile(lines)
    new_file = os.path.join(dst_dir, os.path.splitext(file)[0] + '.dot')
    graph.write_raw(new_file)



