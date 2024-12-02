import numpy as np
from mosaic_algorithms.online_frontier_repair.getMapNeighbours import getMapNeighbours
def getFrontierTiles(map):
    """
    Given a grid of points (list of lists), this function outputs the set of points
    that have less than 8 neighbors in the grid.

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022

    Usage:        frontier, indel = getFrontierTiles(map)

    Inputs:
      > map:          cell matrix of grid points. In order to avoid
                      mapping boundaries, map is bounded by NaN rows and
                      columns (first and last)

    Outputs:
      > frontier:     cell array that contains the frontier tiles in the map
      > indel:        cell array that contains the indices where the
                      frontier tiles are located in 'map'
    """

    # Pre-allocate variables
    frontier = []
    indel = []

    for j in range(len(map[0])):
        for i in range(len(map)):
            tile = map[i][j]

            if not np.isnan(tile).any():
                # Get the observation neighbors of the current observation points
                n = getMapNeighbours(i, j, map)

                # If the observation point has less than 8 planned neighbors, then
                # it is regarded as a frontier tile
                if len(n) < 8:
                    frontier.append(tile)
                    indel.append([i, j])

    return frontier, indel