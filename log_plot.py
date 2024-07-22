"""
Opens log_runs.txt
Expects 10 coloumns
Assigns each column an identity
Plots each coloumn wrt the first coloumn
Save
"""
import numpy as np
import matplotlib.pyplot as plt

col_details = [
    [int, "Nx"],
    [float, "Error"],
    [float, "T0"],
    [float, "Triangluation time"],
    [float, "Plotting time of the exact solution"],
    [float, "Dictionary forming time"],
    [float, "Known-unknown separation time"],
    [float, "np matrix forming time"],
    [float, "Time taken by linear solver"],
    [float, "Beg time junk"],
    ]

for col in col_details:
    col.append([])

with open("log_runs.txt", "r") as log_file:
    for line in log_file:
        for i,entry in enumerate(line.split()):
            print(i, entry, col_details[i][1])
            col_details[i][2].append( col_details[i][0](entry) )

for i_col, col in enumerate(col_details):
    title = col[1]
    vals = col[2]
    plt.title(title)
    plt.plot(col_details[0][2], vals)
    plt.xlabel("Nx")
    plt.grid()
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(f"log_plot{i_col}_NX.png")
    plt.close()
