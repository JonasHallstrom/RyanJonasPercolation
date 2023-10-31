# Jonas Hallstrom, 10/08/2023

# I want to implement a way of drawing bond percolation systems that is generic
#  for different types of 2D lattices. So we'll use a unit-cell based system.
# The order of the bonds in a unit cell should start from the top and go
#  clockwise, just like a clock.
# Updated 10/31/2023 to only draw two clusters

# Imports
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.collections import LineCollection

# Definitions
def draw_square_bonds(bond_bools, groups, graph_params=(4, 20, 4)):
    """Drawing square lattice bond percolation system

    Note that this function will not enforce any boundary conditions, the code
      that generates the bond_bools should be in charge of that, i.e. making
      sure that "bonds" that go outside the BCs are never open, so always False.
    The square lattice unit cell is a site in the top left and then a full bond
        down and a full bond to the right. So 2 bonds per unit cell, and the unit
        cells combine as a simple/primitive square lattice.
    Additionally, the (0,0) unit cell is the top left one and has a graph
        coordinate of (0,0) for the center of the unit cell

    Parameters
    ----------
    bond_bools : NxWxH array of booleans
        Bools for each whether a bond is open (True) or not (False). N is number
        of bonds per unit cell, W the width of the system in unit cells, and H
        the height of the system in unit cells
    graph_params : 1x3 array of floats
        Figure width, site circle diameter, and bond thickness respectively

    Returns
    -------
    None
    """
    n, width, height = bond_bools.shape
    lat_param = 1  # lattice parameter, length of unit cell
    # list of positions of sites within a unit cell
    site_cellpos = [lat_param*np.array([-1/2, 1/2])]  # 1 site per unit cell
    bond_cellseg = []  # 2 bonds/line segments per unit cell
    bond_cellseg.append(lat_param*np.array([[-1/2, 1/2], [1/2, 1/2]]))
    bond_cellseg.append(lat_param*np.array([[-1/2, 1/2], [-1/2, -1/2]]))

    graph_width = graph_params[0]
    plt.rcParams['figure.figsize'] = [graph_width, height*graph_width/width]
    fig, ax = plt.subplots()
    ax.set_axis_off()
    ax.set_xlim(-lat_param, (width-1)*lat_param)
    ax.set_ylim(-(height-1)*lat_param, lat_param)

    A = lat_param*np.arange(width)
    B = -lat_param*np.arange(height)
    # list of positions of cell centers in graph/plot
    cell_graphpos = lat_param*np.array([(x, y) for y in B for x in A])

    # Draw bonds
    bond_graphseg = np.array([cell+seg for cell in cell_graphpos for seg in bond_cellseg])
    flat_bool = [square_bond_bools[i, j, k] for k in range(height) for j in range(width) for i in range(2)]
    bond_graphseg = bond_graphseg[np.where(flat_bool)[0]]

    groups = -np.ones(len(bond_graphseg))
    available_groups = [num for num in range(1000)]

    for i in range(len(groups)):
        if groups[i] == -1:
            groups[i] = available_groups.pop(0)
        for j in range(i+1, len(groups)):
            if groups[i]==groups[j]:
                pass
            elif j-i>2*(width+1):
                pass
            elif (np.array_equal(bond_graphseg[j, 0, :], bond_graphseg[i, 0, :]) or
            np.array_equal(bond_graphseg[j, 0, :], bond_graphseg[i, 1, :]) or
            np.array_equal(bond_graphseg[j, 1, :], bond_graphseg[i, 0, :]) or
            np.array_equal(bond_graphseg[j, 1, :], bond_graphseg[i, 1, :])):
                if groups[j] == -1:
                    groups[j] = groups[i]
                else:
                    available_groups = [groups[j]] + available_groups
                    groups[np.where(groups == groups[j])[0]] = groups[i]

    num_groups = len(np.unique(groups))
    available_groups = list(np.sort(available_groups))
    while (m:=np.max(groups)) > available_groups[0]:
        available_groups = available_groups + [m]
        groups[np.where(groups == m)[0]] = available_groups.pop(0)
        available_groups = list(np.sort(available_groups))

    group_sizes = [len(np.where(groups == i)[0]) for i in range(num_groups)]
    groups_by_size = np.flip(np.argsort(group_sizes))
    reordered_groups = -np.ones(len(groups))
    for i in range(num_groups):
        reordered_groups[np.where(groups == groups_by_size[i])[0]] = i
    bond_color = ["tab:blue" if num==0 else "tab:grey" for num in reordered_groups]

    line_segments = LineCollection(bond_graphseg, linewidths=graph_params[2],
                                    colors=bond_color, linestyle="solid")
    ax.add_collection(line_segments)

    # Draw sites
    site_color = ["k"]*(height*width*len(site_cellpos))
    # list of positions of sites in graph/plot
    site_graphpos = np.array([cell+site for cell in cell_graphpos for site in site_cellpos])
    plt.scatter(site_graphpos[:, 0], site_graphpos[:, 1], s=graph_params[1], c=site_color, zorder=2.5)

    plt.savefig("SquareBondPercolation_{}x{}_p{:.0f}.png".format(width, height, p*100), dpi=300)
    plt.show()


# Some examples
if __name__ == "__main__":
    print("Beginning examples")

    # Square Example
    width, height = (80, 80)
    # Generating bonds with a certain open probability
    p=0.49
    square_bond_bools = np.random.uniform(0.0, 1.0, (2, width, height)) < p
    # Enforce boundary conditions, where bond 0 is right and bond 1 down
    square_bond_bools[0, -1, :] = False
    square_bond_bools[1, :, -1] = False

    draw_square_bonds(square_bond_bools, False, graph_params=(4, 0.8, 2.5))
