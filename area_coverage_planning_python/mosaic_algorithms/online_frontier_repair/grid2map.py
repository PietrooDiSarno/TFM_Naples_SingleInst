import numpy as np

def grid2map(grid):
    """
    This function creates a map from a given grid. It adds a border of NaN
    values around the entire grid. Within the grid, any lists that are empty
    or have been excluded from observation (because their coverage is deemed
    too small) are also filled with NaN values.

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022

    Usage:        map = grid2map(grid)

    Input:
      > grid:          list of lists that contains a 2-element vector for grid
                       points within the ROI, or None for points outside
                       the ROI or excluded from coverage.

    Returns:
      > map:           list of lists representing the map. It is the input grid
                       with an added border of NaN values. Inside the grid,
                       lists with None are replaced with [NaN NaN]

    """

    # Initialize map with additional rows and columns to place the NaN boundaries
    map_height = len(grid) + 2
    map_width = len(grid[0]) + 2
    map = [[[] for _ in range(map_width)] for _ in range(map_height)]


    # Populate the map with data from the grid and NaN borders
    for i in range(map_height):
        for j in range(map_width):
            if i == 0 or j == 0 or i == map_height - 1 or j == map_width - 1 or grid[i - 1][j - 1] is None:
                map[i][j] = np.array([np.nan,np.nan])
            else:
                map[i][j] = grid[i - 1][j - 1]

    return map