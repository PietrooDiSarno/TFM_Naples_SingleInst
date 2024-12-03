import numpy as np

def boustrophedon(grid, dir1, dir2):
    """
    This function plans an observation tour over a specified grid, creating a
    path that covers the area in alternating rows/columns, according to a
    specified input direction.

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022

    Usage:        tour = boustrophedon(grid, dir1, dir2)

    Inputs:
      > grid:        2D list where each cell contains the coordinates of
                     an observation point or is empty if there is no point.
      > dir1:        primary direction of the sweep ('north', 'south', 'east', 'west').
      > dir2:        secondary direction of the sweep. This defines if the
                     traversal is going to be performed either in
                     alternating rows or columns.

    Outputs:
      > tour:        ordered list of points representing the planned
                     tour. Each element is a 2-element numpy.array indicating
                     a point on the grid to be observed.
    """

    # Previous check...
    if dir1 in ['north', 'south']:
        if dir2 not in ['east', 'west']:
            raise ValueError("Sweeping direction is not well defined")
    elif dir1 in ['east', 'west']:
        if dir2 not in ['north', 'south']:
            raise ValueError("Sweeping direction is not well defined")

    # Pre-allocate variables
    tour = []

    # Particular case: len(grid) == 1

    if len(grid) == 1 and len(grid[0]) == 1:
        tour=grid[0]

    sweep = dir2 in ['east', 'south']

    # Plan tour over the grid discretization
    # The origin of the coverage path depends on the spacecraft ground track position
    if dir1 in ['north', 'south']:  # Horizontal sweep
        if sweep:
            bearing = True #left -> right
        else:
            bearing = False #right -> left

        tour = [[] for _ in range(np.count_nonzero([item is not None for row in grid for item in row]))] #list of planned observations
        ii = 0
        for i in range(len(grid)):
            #Sweep across latitude
            if dir1 == 'south':
                irow = i
            else:
                irow = len(grid) - i - 1

            for j in range(len(grid[0])):
                if not bearing:
                    icol = len(grid[0]) - j - 1
                else:
                    icol = j

                if grid[irow][icol] is not None:
                    x, y = grid[irow][icol][0], grid[irow][icol][1]
                    tour[ii] = np.array([x, y])#Save it in the coverage tour
                    ii += 1

            bearing = not bearing  # Switch coverage direction after each row sweeping, i.e. left (highest lon) to right
            # (lowest lon) or vice versa

    elif dir1 in ['east', 'west']:  # Vertical sweep
        if sweep:
            bearing = True # top -> down
        else:
            bearing = False # down -> top

        tour = [[] for _ in range(np.count_nonzero([item is not None for row in grid for item in row]))] #list of planned observations
        ii = 0
        if grid != []:
            for i in range(len(grid[0])):
                if dir1 == 'west':
                    icol = len(grid[0]) - i - 1
                else:
                    icol = i

                for j in range(len(grid)):
                    if bearing:
                        irow = j
                    else:
                        irow = len(grid) - j - 1

                    if grid[irow][icol] is not None:
                        x, y = grid[irow][icol][0], grid[irow][icol][1]
                        tour[ii] = np.array([x, y])
                        ii += 1

                bearing = not bearing  # Switch coverage direction after each column sweeping, i.e. up (highest lat) to down
                # (lowest lat) or vice versa

    return tour

