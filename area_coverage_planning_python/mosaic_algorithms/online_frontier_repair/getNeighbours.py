import numpy as np
def getNeighbours(gamma, ind, w, h, olapx, olapy, dx, dy):
    """
    Given a point, this function outputs the 8 adjacent neighbours and their
    index location in the grid.

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022
    Last Rev.:    06/2023

    Usage:        n, nind = getNeighbours(gamma, ind, w, h, olapx, olapy, dx, dy)

    Inputs:
      > gamma:        2D point to which the 8 adjacent points are needed
      > ind:          index position of gamma in the grid
      > w:            width (x-direction) of the grid spacing
      > h:            height (y-direction) of the grid spacing
      > olapx:        grid footprint overlap in the x direction,
                      in percentage (width)
      > olapy:        grid footprint overlap in the y direction,
                      in percentage (height)
      > dx:           vector that expresses the x-direction in the grid
      > dy:           vector that expresses the y-direction in the grid

    Outputs:
      > n:            list with the 8 neighbouring points
      > nind:         list with the index position of the neighbours
                      in the grid
    """

    # Pre-allocate variables
    n = [None] * 8
    ovlapx = olapx * w / 100
    ovlapy = olapy * h / 100

    # Cardinal and diagonal neighbouring points
    n[0] = gamma + (w - ovlapx) * dx + (h - ovlapy) * dy  # northeast
    n[1] = gamma + (w - ovlapx) * dx                      # east
    n[2] = gamma + (w - ovlapx) * dx + (-h + ovlapy) * dy # southeast
    n[3] = gamma + (h - ovlapy) * dy                      # north
    n[4] = gamma + (-h + ovlapy) * dy                     # south
    n[5] = gamma + (-w + ovlapx) * dx + (h - ovlapy) * dy # northwest
    n[6] = gamma + (-w + ovlapx) * dx                      # west
    n[7] = gamma + (-w + ovlapx) * dx + (-h + ovlapy) * dy # southwest

    # Cardinal and diagonal neighbouring points grid indices
    nind = [None] * 8
    ind=np.array(ind)
    nind[0] = ind + np.array([-1, 1])
    nind[1] = ind + np.array([0, 1])
    nind[2] = ind + np.array([1, 1])
    nind[3] = ind + np.array([-1, 0])
    nind[4] = ind + np.array([1, 0])
    nind[5] = ind + np.array([-1, -1])
    nind[6] = ind + np.array([0, -1])
    nind[7] = ind + np.array([1, -1])

    return n, nind
