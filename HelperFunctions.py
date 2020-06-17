import numpy as np
from math import exp


def neighbours(t, i, j, grid_size):
    """
    Generates a Moore neighbourhood around a cell (i,j)
    :param t: time coordinate
    :param i: x-axis coordinate
    :param j: y-axis coordinate
    :param grid_size: total size of the spatial grid  
    :return: Moore neighbourhood of cell (t,i,j)
    """
    if i >= grid_size or j >= grid_size:
        raise Exception('I or J is not a valid grid entry')

    rowbound = grid_size
    colbound = grid_size

    # Define Moore Neighbourhood
    neighbours = [
        (t, i - 1, j - 1), (t, i - 1, j), (t, i - 1, j + 1),
        (t, i, j - 1), (t, i, j + 1),
        (t, i + 1, j - 1), (t, i + 1, j), (t, i + 1, j + 1)
    ]
    # Filter the valid neighbours.
    valid_tuples = []
    for (t, a, b) in neighbours:
        if -1 < a < rowbound and -1 < b < colbound:
            valid_tuples.append((t, a, b))

    return valid_tuples


def neighbour_summation(t, i, j, susceptibles, neighbour_dict):
    """
    Computes the susceptible ratio from Eq. 8 of Aadland et al. (2015)
    :param t: time coordinate
    :param i: x-axis coordinate
    :param j: y-axis coordinate
    :param neighbour_dict: Dictionary containing all of the neighbours of cell (t,i,j)
    :return: List of susceptible ratios for each neighbour
    """
    S = susceptibles
    coordinate = (t, i, j)

    # Sum all of the susceptible trees for the valid neighbours of cell (i,j).
    neighbour_trees = 0
    for (a, b, c) in neighbour_dict.get(coordinate):
        neighbour_trees += S[t, b, c]

    solution_list = [S[t, o, p] / neighbour_trees if neighbour_trees != 0 else 0 for (n, o, p) in neighbour_dict.get(coordinate)]
    return solution_list

def neighbours_ldd(i, j, grid_size):
    """
    Similar to neighbours(t,i,j,grid_size) just without the time variable
    """
    # Dimension is the length of a row in square matrix
    rowbound = grid_size
    colbound = grid_size

    # Define Moore Neighbourhood
    neighbours = [
        (i - 1, j - 1), (i - 1, j), (i - 1, j + 1),
        (i, j - 1), (i, j + 1),
        (i + 1, j - 1), (i + 1, j), (i + 1, j + 1)
    ]
    # Filter the valid neighbours.
    valid_tuples = []
    for (a, b) in neighbours:
        if -1 < a < rowbound and -1 < b < colbound:
            valid_tuples.append((a, b))

    return valid_tuples

def neighbour_check(valid_tuples_of_neighbours_ldd, grid3):
    """
    Generates a list containing appropriate values from grid3
    :param valid_tuples_of_neighbours_ldd: The return value from neighbours_ldd(i,j,grid_size)
    :param grid3: See MPB_Climate_Sim.py
    :return: A list of values that correspond to the valid_tuples_of_neighbours_ldd
    """
    valid_tuples = valid_tuples_of_neighbours_ldd
    local_grid = grid3

    soln_list = [local_grid[i][j]  for i,j in valid_tuples]

    return soln_list

def ldd_probability(grid3, grid_size):
    """

    :param grid3: See MPB_Climate_Sim.py
    :param grid_size:
    :return: The probability that cell (i,j) will experience long distance dispersal in the grid.
    """
    counter = 0
    local_grid = grid3.copy()

    for i in range(grid_size):
        for j in range(grid_size):

            if local_grid[i][j] == 1:
                valid_neighbours = neighbours_ldd(i, j, grid_size)
                valid_neighbour_vals = neighbour_check(valid_neighbours,local_grid)

                if 2 not in valid_neighbour_vals:
                    local_grid[i][j] = 3
                    counter += 1

    ldd_prob = counter / local_grid.size
    return ldd_prob, local_grid


def poisson_lambda(grid1, grid2, dimension):
    """
    Returns the Poisson_Lambda for the Poisson Distribution
    :param grid1:
    :param grid2:
    :param dimension:
    :return:
    """
    cell_count, cell_values = 0, 0
    tuple_list = []
    for i in range(dimension):
        for j in range(dimension):
            if grid1[i][j] == 3:
                cell_values += grid2[i][j]
                cell_count += 1
                tuple_list.append((i, j))

    poisson_lambda = cell_values / cell_count
    return poisson_lambda


def func3(t, i, j, neighbour_dict, target_index_dict, local_dict):
    coord = (t, i, j)
    neighbours = neighbour_dict.get(coord)
    target_indicies = target_index_dict.get(coord)

    soln_list = [local_dict[neighbours[i]][target_indicies[i]] for i, _ in enumerate(neighbours)]
    return soln_list