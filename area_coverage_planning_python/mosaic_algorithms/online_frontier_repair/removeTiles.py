import numpy as np


def removeTiles(map, tiles):
    """
    This function removes disposable observation points within the grid.

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         06/2023

    Usage:        map = removeTiles(map, tiles)

    Inputs:
    > map_grid:   list of lists representing grid points. In order to avoid
                  mapping boundaries, map is bounded by NaN rows and
                  columns (first and last)
    > tiles:      list of disposable observation points to be
                  removed from 'tour' and 'grid'

    Outputs:
    > map:   updated list of lists representing grid points
    """

    for i in range(len(tiles)):
        # For each observation point in the removal list...

        for ii in range(len(map)):
            for jj in range(len(map[0])):
                if np.linalg.norm(map[ii][jj] - tiles[i]) < 1e-5:
                    
                    # Remove elements from the grid ([NaN, NaN])
                    map[ii][jj] = [np.nan, np.nan]

    return map