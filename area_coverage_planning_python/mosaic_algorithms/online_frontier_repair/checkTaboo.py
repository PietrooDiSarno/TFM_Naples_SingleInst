import numpy as np
from mosaic_algorithms.online_frontier_repair.map2grid import map2grid

def checkTaboo(N, Nind, map, ind_row, ind_col, indir1, indir2):
    """
    This function evaluates each tile in a list of potential observation
    points (N, Nind) and determines whether it should be considered taboo,
    i.e., unsuitable for selection based on its location relative to a specified
    direction of movement (indir1, indir2) and starting point (ind_row, ind_col).
    Taboo tiles are those that do not conform to the expected movement
    pattern across the grid.
    [Future work]: get rid of boustrophedonMod function

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022

    Usage:        N, Nind = check_taboo(N, Nind, map, ind_row, ind_col, indir1, indir2)

    Inputs:
        > N:            list containing the values of potential
                        observation pointsF
        > Nind:         list containing the indices of potential
                        observation points within the map
        > map:          list of list representing grid points, where first and last rows
                        and columns are NaN to denote boundaries
        > ind_row:      starting row index for evaluating taboo conditions
        > ind_col:      starting column index for evaluating taboo conditions
        > indir1:       primary movement direction in the path plan
        > indir2:       secondary movement direction, orthogonal to indir1

    Outputs:
        > N, Nind:      updated lists
    """

    # Persistent variables in a way suitable for Python
    if 'pdir1' not in checkTaboo.__dict__:
        checkTaboo.pdir1 = indir1
        checkTaboo.pdir2 = indir2

    pdir1 = checkTaboo.pdir1
    pdir2 = checkTaboo.pdir2
    ## Previous checks...
    #if not N or not Nind:
        #return

    # Pre-allocate variables
    nindel = np.array([])
    grid = map2grid(map)
    dir1, dir2 = boustrophedonMod(grid, indir1, indir2)
    dir_change = False
    if pdir2 != indir2:
        dir_change = True

    # Previous checks...
    if not N or not Nind:
        return N, Nind

    # Define the grid boundaries (ending rows and columns in the map)
    for i in range(len(map) - 1, -1, -1):
        el = next((j for j, val in enumerate(map[i]) if not (np.isnan(val)).any()), None) #get non-NaN elements in the map
        if el is not None:
            Nrow = i
            break

    for i in range(len(map[0]) - 1, -1, -1):
        aux = []
        for k in range(len(map)):
            aux.append(map[k][i])
        el = next((j for j, val in enumerate (aux) if not (np.isnan(val)).any()), None)
        if el is not None:
            Ncol = i
            break

    # Define the grid boundaries (starting rows and columns in the map)
    for i in range(len(map)):
        el = next((j for j, val in enumerate(map[i]) if not (np.isnan(val)).any()), None)
        if el is not None:
            Orow = i
            break


    for i in range(len(map[0])):
        aux = []
        for k in range(len(map)):
            aux.append(map[k][i])
        el = next((j for j, val in enumerate(aux) if not (np.isnan(val)).any()), None)
        if el is not None:
            Ocol = i
            break

    # Taboo tiles search
    for i in range(len(Nind)):
        taboo = False
        indel = Nind[i]

        if dir1 in ['north', 'south']:  # horizontal sweep

            if dir1 == 'south':  # spacecraft is towards roi's bottom
                if indel[0] < ind_row:
                    taboo = True
            else:  # spacecraft is towards roi's top
                if indel[0] > ind_row:
                    taboo = True

            if dir2 == 'west':  # tour is moving to the right (left -> right dir.)
                if indel[0] == ind_row and indel[1] > ind_col:
                    if dir_change and ind_col < Ncol:
                        taboo = True
                    elif not dir_change:
                        taboo = True
            else:  # tour is moving to the left (right -> left dir.)
                if indel[0] == ind_row and indel[1] < ind_col:
                    if dir_change and ind_col > Ocol:
                        taboo = True
                    elif not dir_change:
                        taboo = True

        elif dir1 in ['east', 'west']:  # vertical sweep

            if dir1 == 'east':
                if indel[1] < ind_col:
                    taboo = True
            else:
                if indel[1] > ind_col:
                    taboo = True

            if dir2 == 'south':  # downsweep = true. tour is moving to the bottom
                if indel[1] == ind_col and indel[0] < ind_row:
                    if dir_change and ind_row > Orow:
                        taboo = True
                    elif not dir_change:
                        taboo = True
            else:  # tour is moving to the top
                if indel[1] == ind_col and indel[0] > ind_row:
                    if dir_change and ind_row < Nrow:
                        taboo = True
                    elif not dir_change:
                        taboo = True

        if taboo:
            nindel = np.append(nindel, i)

    # Delete taboo tiles
    N = [N[j] for j in range(len(N)) if j not in nindel]
    Nind = [Nind[j] for j in range(len(Nind)) if j not in nindel]

    checkTaboo.pdir1 = pdir1
    checkTaboo.pdir2 = pdir2

    return N, Nind


def boustrophedonMod(grid, dir1, dir2):

    #Previous check...
    if dir1 in ['north', 'south']:
        if dir2 not in ['east', 'west']:
            raise ValueError("Sweeping direction is not well defined")
    elif dir1 in ['east', 'west']:
        if dir2 not in ['north', 'south']:
            raise ValueError("Sweeping direction is not well defined")

    # Pre-allocate variables
    sweep = dir2 in ['east', 'south']
    currdir1, currdir2 = dir1, dir2

    # Plan tour over the grid discretization
    # The origin of the coverage path depends on the spacecraft ground track position
    start = False
    if dir1 in ['north', 'south']:  # Horizontal sweep

        if sweep:
            bearing = True # left -> right
        else:
            bearing = False # right -> left

        for i in range(len(grid)):
            # Sweep across latitude
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
                    start = True
                    break
            if start:
                break
            bearing = not bearing  # Switch coverage direction after each row sweeping, i.e. left (highest lon) to right
            # (lowest lon) or vice versa


    elif dir1 in ['east', 'west']:  # Vertical sweep

        if sweep:
            bearing = True  # top -> down
        else:
            bearing = False  # down -> top

        for i in range(len(grid[0])):
            # Sweep across longitude
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
                    start = True
                    break
            if start:
                break
            bearing = not bearing  # Switch coverage direction after each row sweeping, i.e. left (highest lon) to right
            # (lowest lon) or vice versa

    # Adjust direction after sweeping
    if bearing:
        if dir2 == 'west':
            currdir2 = 'east'
        elif dir2 == 'north':
            currdir2 = 'south'
    else:
        if dir2 == 'east':
            currdir2 = 'west'
        elif dir2 == 'south':
            currdir2 = 'north'

    return currdir1, currdir2



# [Future work]: get rid of boustrophedonMod
# Shift directions
# if dir2 == 'west':
#     dir2 = 'east'
# elif dir2 == 'east':
#     dir2 = 'west'
# elif dir2 == 'north':
#     dir2 = 'south'
# elif dir2 == 'south':
#     dir2 = 'north'

# Identify in which direction are we moving
# if dir1 in ['north', 'south']:
#     if dir2 == 'west':
#         if len(map) % 2:
#             if (ind_row - 1) % 2 == 0:
#                 dir2 = 'east'
#         else:
#             if (ind_row - 1) % 2 != 0:
#                 dir2 = 'east'
#     else:
#         if len(map) % 2:
#             if (ind_row - 1) % 2 == 0:
#                 dir2 = 'west'
#         else:  n
#             if (ind_row - 1) % 2 != 0:
#                 dir2 = 'west'
# else:
#     if dir2 == 'north':
#         if len(map[0]) % 2:
#             if (ind_col - 1) % 2 == 0:
#                 dir2 = 'south'
#         else:
#             if (ind_col - 1) % 2 != 0:
#                 dir2 = 'south'
#     else:
#         if len(map[0]) % 2:
#             if (ind_col - 1) % 2 == 0:
#                 dir2 = 'north'
#         else:  # The number of columns is even
#             if (ind_col - 1) % 2 != 0:
#                 dir2 = 'north'



