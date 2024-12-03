import copy

import numpy as np
from shapely.geometry import MultiPolygon, Polygon
import matplotlib.pyplot as plt
from conversion_functions import *
from mosaic_algorithms.online_frontier_repair.checkTaboo import checkTaboo
from mosaic_algorithms.online_frontier_repair.removeTiles import removeTiles
from mosaic_algorithms.online_frontier_repair.insertTiles import insertTiles
from mosaic_algorithms.online_frontier_repair.grid2map import grid2map
from mosaic_algorithms.online_frontier_repair.map2grid import map2grid
from mosaic_algorithms.online_frontier_repair.getFrontierTiles import getFrontierTiles
from mosaic_algorithms.online_frontier_repair.getNeighbours import getNeighbours
from mosaic_algorithms.auxiliar_functions.grid_functions.inst2topo import inst2topo
from mosaic_algorithms.auxiliar_functions.grid_functions.topo2inst import topo2inst
from mosaic_algorithms.auxiliar_functions.grid_functions.boustrophedon import boustrophedon

fpref = None
pointing0 = None
sweepDir1 = None
sweepDir2 = None

def updateGrid(roi, inst_tour, inst_grid, grid_dirx, grid_diry, cx, cy, olapx, olapy,
               insweepDir1, insweepDir2, seed, old_seed, gamma, et, inst, sc, target):
    """
    This function dynamically updates the grid of observations by
    incorporating new observation points, adjusting for changes in the
    observation geometry, and removing points that no longer contribute
    to covering the region of interest (ROI). It uses a boustrophedon pattern
    for traversal. Adapted from [1].

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2023

    Usage:        [seed, inst_grid, inst_tour, topo_tour] = updateGrid(roi,
                   inst_tour, inst_grid, grid_dirx, grid_diry, cx, cy, olapx, olapy,
                   insweepDir1, insweepDir2, seed, old_seed, gamma, et, inst, sc, target)

    Inputs:
      > roi:          matrix containing the vertices of the uncovered area
                      of the ROI polygon. The vertex points are expressed in
                      2D, in latitudinal coordinates [º]
      > inst_tour:    tour path in instrument frame coordinates
      > inst_grid:    grid of potential observation points in instrument
                      frame coordinates
      > grid_dirx:    direction of the grid along the x-axis in the
                      instrument frame
      > grid_diry:    direction of the grid along the y-axis in the
                      instrument frame
      > cx, cy:       centroid coordinates of the roi
      > olapx:        grid footprint overlap in the x direction (longitude),
                      in percentage (width)
      > olapy:        grid footprint overlap in the y direction (latitude),
                      in percentage (height)
      > insweepDir1, insweepDir2: initial directions defining the sweep of the
                              Boustrophedon decomposition, derived according
                              to the spacecraft's position with respect to
                              the ROI
      > seed:         current starting point for the grid update
      > old_seed:     starting point from the previous iteration
      > gamma:        current observation point
      > et:           current time in ephemeris seconds past J2000 epoch
      > inst:         string name of the instrument
      > sc:           string name of the spacecraft
      > target:       string name of the target body

    Outputs:
      > seed, inst_grid, inst_tour: updated variables
      > topo_tour:    tour path in topographical coordinates (lat/lon on the
                      target body), in [deg]

    [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G.,
    Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis
    steerable, 2D framing sensors.
    """

    # Pre-allocate variables...
    global fpref, pointing0, sweepDir1, sweepDir2

    if sweepDir1 is None:
        sweepDir1 = insweepDir1
        sweepDir2 = insweepDir2

    if pointing0 is None:
        pointing0 = np.array([cx, cy])

    cin = []  # set of new observations (inside or outside from 'tour')
    cind = []  # map indices of the new observation (potential placement in the map)
    cout = []  # set of disposable observations (inside or outside from 'tour')
    N = []  # set of new tiles (outside from 'tour')
    Nind = []  # map indices of the new tiles
    X = []  # set of disposable tiles (inside of 'tour')
    epsilon = 0.02

    # Build reference tile (it's always going to be the same in subsequent calls)
    if fpref is None:
        # We need to craft the tile reference that we are going to use throughout the heuristic operations. To avoid
        # repetitions, we initialize fpref to None and calculate it only on the first call of the function, ensuring
        # that it is set only once.

        _,_,_,bounds = mat2py_getfov(mat2py_bodn2c(inst)[0], 4)  # get fovbounds in the instrument's reference frame
        maxx, minx = max(bounds[0]), min(bounds[0])
        maxy, miny = max(bounds[1]), min(bounds[1])
        width = maxx - minx
        height = maxy - miny
        xlimit = [-width / 2, width / 2]
        ylimit = [-height / 2, height / 2]

        xbox = np.array([xlimit[0], xlimit[0], xlimit[1], xlimit[1], xlimit[0]])
        ybox = np.array([ylimit[0], ylimit[1], ylimit[1], ylimit[0], ylimit[0]])

        fpref = {
            'width': width,
            'height': height,
            'bvertices': np.column_stack((xbox, ybox))
        }

    # Project ROI topographical coordinates to instrument's focal plane
    targetArea = topo2inst(roi, cx, cy, target, sc, inst, et)
    if (np.isnan(targetArea[:, 0])).any():
        nanindex = np.where(np.isnan(targetArea[:, 0]))[0]
        polygon_list = []
        for i in range(len(nanindex)):
            if i == 0:
                polygon_list.append(Polygon(list(zip(targetArea[:nanindex[0], 0], targetArea[:nanindex[0], 1]))))
            else:
                polygon_list.append(Polygon(
                    list(zip(targetArea[nanindex[i - 1] + 1:nanindex[i], 0], targetArea[nanindex[i - 1] + 1:nanindex[i], 1]))))
        if ~ np.isnan(targetArea[-1, 0]):
            polygon_list.append(Polygon(list(zip(targetArea[nanindex[-1] + 1:, 0], targetArea[nanindex[-1] + 1:, 1]))))
        targetpshape = MultiPolygon(polygon_list)
    else:
        targetpshape = Polygon(zip(targetArea[:, 0], targetArea[:, 1]))
    targetpshape = targetpshape.buffer(0)
    # [Future work]: orientation angle may change over the course of the mosaic
    # targetArea = topo2inst(roi, gamma_topo[0], gamma_topo[1], target, sc, inst, et) #current roi coordinates in the
    # instrument reference frame, when the instrument is pointing at the current grid origin point (next observation)
    # Since we are not calculating the grid again, but rather updating a reference, we need to have the same origin,
    # i.e., we need a 0 point of the ROI's projection that is invariant across iterations, so there is a unique
    # correspondence of the grid points over time.
    # cent = topo2inst(pointing0, gamma_topo[0], gamma_topo[1], target, sc, inst, et)
    # targetArea[:, 0] -= cent[0]
    # targetArea[:, 1] -= cent[1]
    ## Oriented area
    # anglerot = -(angle - angle0)
    # rotmat = np.array([[np.cos(np.radians(anglerot)), -np.sin(np.radians(anglerot))],
    #                   [np.sin(np.radians(anglerot)), np.cos(np.radians(anglerot))]])

    ## matrixGrid directions x and y
    # polygon = Polygon(targetArea)
    # centroid = polygon.centroid
    # cxt, cyt = centroid.x, centroid.y
    # orientedArea = np.zeros_like(len(targetArea),2)
    # for j in range(len(targetArea)):
    #     orientedArea[j, :] = np.array([cxt, cyt]) + rotmat @ (targetArea[j, :] - np.array([cxt, cyt]))

    # target_polygon = Polygon(orientedArea) # Create a polygon shape

    # Get grid shifting due to observation geometry update
    updated_seed = (topo2inst(np.array([np.array(gamma)]), cx, cy, target, sc, inst, et))

    if not np.isnan(updated_seed).all() and np.size(updated_seed)!=0:
        updated_seed = updated_seed[0]
        shift = updated_seed - np.array(seed)
    else:
        shift = 0
        updated_seed = np.array(seed)
        print("Non existent gamma")
    seed = updated_seed
    old_seed += shift

    # Shift grid and tour
    for i in range(len(inst_grid)):
        for j in range(len(inst_grid[i])):
            if inst_grid[i][j] is not None:
                inst_grid[i][j] += shift

    for i in range(len(inst_tour)):
        inst_tour[i] += shift

    # UPDATE GRID
    # Update map by removing the previous element in the tour (next observation)
    # Find which position does gamma occupy in this grid
    ind_row, ind_col = None, None
    for j in range(len(inst_grid[0])):
        for i in range(len(inst_grid)):
            if inst_grid[i][j] is not None and np.linalg.norm(inst_grid[i][j] - old_seed) < 1e-3:
                ind_row, ind_col = i, j
                break
    if old_seed is not None:
        inst_grid[ind_row][ind_col] = None

    # Enclose grid in a bigger matrix (with first and last rows and columns
    # with NaN values, so we can explore neighbours adequately)
    map = grid2map(inst_grid)

    # Obtain the frontier tiles in the map: points that have less than 8
    # neighbours in the grid
    frontier, indel = getFrontierTiles(map)

    # Update grid
    openList = copy.deepcopy(frontier) # open list starts as frontier set F
    s = copy.deepcopy(openList)  # seeds: initial open list (frontier tiles)

    while openList:

        # For each element in the openList array of grid points...
        o = openList[0]  # lon-lat coordinates of the observation point
        currind = indel[0]  # row-column indices of the observation point in 'map'

        # Visited elements are deleted
        indel.pop(0)
        openList.pop(0)

        # Pre-allocate variables
        membershipChanged = False  # boolean variable that indicates if the observation point has changed its membership
        # in 'tour'
        insideTour = False  # boolean variable that indicates if the observation point belongs to 'tour'
        inSeed = False  # boolean variable that determines if the observation is in the seeding list

        # Previous check: is the current observation point inside the seeding list? in case it is, let's analyze the
        # neighbouring points (after checking their membership in tour, next point)
        for i in range(len(s)):
            if np.array_equal(s[i], o):
                inSeed = True
                break

        # Analyze the current element's membership in tour
        # Compute the current footprint's covered area
        aux = np.array(o) + np.array(fpref['bvertices'])
        if (np.isnan(aux[:, 0])).any():
            nanindex = np.where(np.isnan(aux[:, 0]))[0]
            polygon_list = []
            for i in range(len(nanindex)):
                if i == 0:
                    polygon_list.append(Polygon(list(zip(aux[:nanindex[0], 0], aux[:nanindex[0], 1]))))
                else:
                    polygon_list.append(Polygon(
                        list(zip(aux[nanindex[i - 1] + 1:nanindex[i], 0], aux[nanindex[i - 1] + 1:nanindex[i], 1]))))
            if ~ np.isnan(aux[-1, 0]):
                polygon_list.append(Polygon(list(zip(aux[nanindex[-1] + 1:, 0], aux[nanindex[-1] + 1:, 1]))))
            fpshape = MultiPolygon(polygon_list)
        else:
            fpshape = Polygon((list(zip(aux[:, 0], aux[:, 1]))))

        fpshape = fpshape.buffer(0)

        inter = (targetpshape.difference(fpshape)).buffer(0)
        areaI = inter.area
        areaT = targetpshape.area
        fpArea = fpshape.area

        if (areaT - areaI) / fpArea >= epsilon:  # if the observation covers at least a minimum ROI area
            cin.append(o)  # Add it to the list of covering tiles
            cind.append(currind)

            # Identify if the observation was already included in the planned tour
            for i in range(len(inst_tour)):
                if np.array_equal(o, inst_tour[i].reshape(1,2)):
                    insideTour = True
                    break
            if not insideTour:
                # If it wasn't included, then its membership changed
                membershipChanged = True
        else:  # Otherwise (the observation's footprint falls outside the ROI)
            cout.append(o)  # Add it to the list of disposal tiles

            # Identify if the observation was already included in the planned tour
            for i in range(len(inst_tour)):
                if np.array_equal(o, inst_tour[i].reshape(1,2)):
                    insideTour = True
                    break
            if insideTour:
                # If it was included, then its membership changed
                membershipChanged = True

        if inSeed or membershipChanged:
            # Get the observation neighbor elements (diagonal and cardinal)
            n, nind = getNeighbours(o, currind, fpref['width'], fpref['height'], olapx, olapy, grid_dirx,
                                        grid_diry)

            # Check if the neighbors are already inside tour (in that case it is not necessary to include them in
            # the openlist for re-evaluation)
            nindel = []
            for i in range(len(n)):
                inTour=False
                for j in range(len(inst_tour)):
                    if np.linalg.norm(n[i]-inst_tour[j]) < 1e-5:
                        inTour=True
                        break
                if inTour:
                    nindel.append(i)

            n = [n[i] for i in range(len(n)) if i not in nindel]
            nind = [nind[i] for i in range(len(nind)) if i not in nindel]

            for i in range(len(n)):

                # Check if the neighbors are included in the cin or cout sets
                in1=False
                in2=False
                in1 = any(np.linalg.norm(cin_point - n[i]) < 1e-5 for cin_point in cin)
                in2 = any(np.linalg.norm(cout_point - n[i]) < 1e-5 for cout_point in cout)

                # If the neighbour node is not in the cin list nor in the cout list... then add it to the openList
                # for evaluation (if not already included)
                if not in1 and not in2:
                    if not any(np.linalg.norm(open_point - n[i]) < 1e-5 for open_point in openList): # if it's not
                    # in openList, then add it
                        openList.append(n[i])
                        indel.append(nind[i])

    # Identify new tiles N = Cin - Tour
    for i, c in enumerate(cin):
        if not any(np.linalg.norm(c - tour_point) < 1e-5 for tour_point in inst_tour): # if c is checked to be outside, include it in the new tiles set
            N.append(c)
            Nind.append(cind[i])

    # Check that new identified tiles are not taboo  (moving backwards in the coverage path)
    ind_row, ind_col = None, None

    for j in range(len(map[0])):
        for i in range(len(map)):
            if not ((j == 0 and i == 0) or (j == len(map[0]) -1 and i == len(map) -1)):
                if not np.isnan(map[i][j]).all() and np.linalg.norm(map[i][j] - seed) < 1e-3:
                    ind_row, ind_col = i , j
                    break

    # Check that N is not coincident with old_seed...
    for i in range(len(N)):
        if np.linalg.norm(N[i] - old_seed) < 1e-5:
            N.pop(i)
            Nind.pop(i)
            break
    N, Nind = checkTaboo(N, Nind, map, ind_row, ind_col, sweepDir1, sweepDir2)

    # Identify tiles to remove: X = Cout - Tour
    X = []
    for c in cout:
        if any(np.linalg.norm(c - tour_point) < 1e-5 for tour_point in inst_tour): #if c is checked to be
            # in 'tour', include it in the  disposable tiles set
            X.append(c)

    # Remove disposable tiles
    map = removeTiles(map, X)

    # Insert new tiles
    map = insertTiles(map, N, Nind)

    # # Plot grid
    # plt.figure()
    # x, y = targetpshape.exterior.xy
    # plt.plot(x, y)
    # plt.axis('equal')
    # # Loop through each grid point and plot if it exists
    # for i in range(len(grid)):
    #     for j in range(len(grid[0])):
    #         point = grid[i][j]
    #         if point is not None:
    #             plt.plot(point[0], point[1], 'b^')
    # plt.show()

    # Boustrophedon decomposition
    inst_grid = map2grid(map)
    inst_tour = boustrophedon(inst_grid, sweepDir1, sweepDir2)

    if inst_tour:
        topo_tour = inst2topo([inst_tour], cx, cy, target, sc, inst, et)[0]
        # Remove empty elements from the tour, which may result from unobservable
        # regions within the planned path
        emptyCells = [x is None for x in topo_tour]
        indEmpty = [i for i, x in enumerate(emptyCells) if x]  # find indices of empty cells
        for k in range(len(indEmpty)):
            emptyEl = inst_tour[indEmpty[k]]
            for i in range(len(map)):
                for j in range(len(map[0])):
                    if map[i][j] is not None:
                        if np.linalg.norm(map[i][j] - emptyEl) < 1e-5:
                            map[i][j] = None
        topo_tour = [x for i, x in enumerate(topo_tour) if i not in indEmpty]  # remove empty cells
        # Boustrophedon decomposition
        inst_grid = map2grid(map)
        inst_tour = boustrophedon(inst_grid, sweepDir1, sweepDir2)
        if inst_tour:
            seed = inst_tour[0]
        else:
            seed = None

        # inst_tour = [x for i, x in enumerate(inst_tour) if i not in indEmpty]  # remove empty cells
        # seed = inst_tour[0]
        # for i in range(len(emptyCells)):
        #     if not emptyCells[i]:
        #         seed = inst_tour[i]
        #         break

    else:
        seed = None
        topo_tour = []

    return seed, inst_grid, inst_tour, topo_tour


