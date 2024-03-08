import matplotlib.pyplot as plt
import pandas as pd
import numpy as np  # Ensure numpy is imported


# Step 1: Read the Excel file
file_path = 'C:\\Users\\marsh\\Documents\\GitHub\\BachelorQUBOProject\\small_dots_stats.xlsx'
df = pd.read_excel(file_path, usecols=range(8), nrows=23)  # Read the first 8 columns and 23 rows

# Step 2: Sort the data
df = df.iloc[1:].sort_values(by=df.columns[1])
# Assuming df is your DataFrame and is already sorted by the 'Base Value' column

# Generate a numerical range for x-axis positions
x_positions = np.arange(len(df))

plt.figure(figsize=(14, 8))


# Step 3: Overlay line plots for Method A and Method B
# plt.plot(x_positions, df['Nr. of total Bin Vars (Qutie)'], color='red', marker='o', linestyle='-', linewidth=2, markersize=5, label='Nr. of total Bin Vars (Qutie)')
plt.bar(x_positions, df['Nr. of total Bin Vars (Qutie)'], color='red',label='Nr. of total Bin Vars (Qutie)', width=0.4)
plt.bar(x_positions, df['Nr. of total Bin Vars (new Method)'], color='green', label='Nr. of total Bin Vars (new Method)', width=0.4)
plt.bar(x_positions, df['Nr. of Nodes'], color='skyblue', tick_label=df['File Name'] label='Nr. of Nodes', width=0.4)

# Adding some chart elements for clarity
plt.title('Comparison of Base Values and Total Values by Qutie and the new Method')
plt.xlabel('Metabolic Networks (small)')
plt.ylabel('Nr of Binary Variables')
plt.xticks(x_positions, df['Nr. of Nodes'], rotation=45)  # Set x-ticks to be the sorted row identifiers; rotate for readability
plt.legend()

# Show the plot
plt.tight_layout()
plt.show()

