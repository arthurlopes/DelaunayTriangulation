import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import pylab as pl
from matplotlib import collections  as mc
from tqdm import tqdm
import ast

with open('points.txt', 'r') as file:
    # Initialize an empty list to store the triangle points
    points = []

    # Loop through each line in the file
    for line in file:
        line = line[:-2]
        # Convert the string representation of the tuple to a Python tuple object
        point = ast.literal_eval(line)

        # Add the triangle to the list of triangles
        points.append(point)

with open('triangles.txt', 'r') as file:
    # Initialize an empty list to store the triangle points
    triangles = []

    # Loop through each line in the file
    for line in file:
        line = line[:-2]
        # Convert the string representation of the tuple to a Python tuple object
        triangle = ast.literal_eval(line)

        # Add the triangle to the list of triangles
        triangles.append(triangle)

lines = []
for t in tqdm(triangles):
    lines.append([(t[0][0], t[0][1]), (t[1][0], t[1][1])])
    lines.append([(t[1][0], t[1][1]), (t[2][0], t[2][1])])
    lines.append([(t[2][0], t[2][1]), (t[0][0], t[0][1])])

fig, ax = pl.subplots(figsize=(15,15))

lc = mc.LineCollection(lines, linewidths=1)
ax.add_collection(lc)
ax.autoscale()
plt.savefig('delaunay.png')
