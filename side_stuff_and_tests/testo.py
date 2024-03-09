import matplotlib.pyplot as plt

# Example 1: Basic bar chart to show the use of 'x'
plt.figure(figsize=(14, 4))
plt.subplot(1, 4, 1)
x = [1, 2, 3, 4]  # x coordinates of the bars
height = [10, 15, 7, 10]
aigw = [3, 9,4,7]
plt.bar(x, height, color='skyblue')
plt.title("Basic Bar Chart\n(x coordinates)")

# Example 2: Demonstration of 'align' parameter
plt.subplot(1, 4, 2)
plt.bar(x, height, width=0.5, align='center', color='lightgreen')
plt.title("Align='center'")

plt.bar([p + 0.5 for p in x], height, width=0.5, align='edge', color='orange', alpha=0.7)
plt.title("Align='center' & 'edge'")

plt.subplot(1, 4, 3)

plt.bar([p - 0.2 for p in x], [p - 2 for p in height], width=0.4, align='center', color='red',edgecolor="black")

plt.bar([p + 0.2 for p in x], [p + 1 for p in height], width=0.4, align='center', color='orange',edgecolor="black")

plt.bar([p - 0.2 for p in x], aigw, width=0.4, align='center', color='red', edgecolor="black")
# alpha=0.2
plt.bar([p + 0.2 for p in x], aigw, width=0.4, align='center', color='orange', edgecolor="black")
plt.title("Align='center'")

# Example 3: Demonstration of 'tick_label'
plt.subplot(1, 4, 4)
tick_labels = ['A', 'B', 'C', 'D']  # Labels for each bar
plt.bar(x, height, tick_label=tick_labels, color='violet')
plt.title("Tick Labels")

plt.show()