import numpy as np
from conversion_functions import *

def map2grid(map):
    """
    This function takes a map, defined as a list of lists where the outermost
    rows and columns are filled with NaNs to denote boundaries or areas
    outside of interest, and converts it into a grid by removing these NaN
    borders. Within the resulting grid, any list containing NaN values is
    emptied, signifying that it does not contain useful data or represents an
    area outside the region of interest.

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022

    Usage:        grid = map2grid(map)

    Inputs:
      > map:          list of lists representing the map. It is the input grid
                      with an added border of NaN values. Inside the grid,
                      empty cells are replaced with [NaN NaN]

    Returns:
      > grid:         list of lists that contains a 2-element vector for grid
                      points within the ROI, or is empty for points outside
                      the ROI or excluded from coverage.
    """
    grid = [row[1:-1] for row in map[1:-1]]
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            if np.size(grid[i][j]) == 1 and grid[i][j] == None:
                grid[i][j] = grid[i][j]
            elif (np.isnan(grid[i][j])).all():
                    grid[i][j] = None


    return grid