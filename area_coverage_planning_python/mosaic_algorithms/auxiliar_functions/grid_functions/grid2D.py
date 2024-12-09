import copy

import numpy as np
from scipy.spatial import ConvexHull
from shapely.geometry import MultiPolygon, Polygon, Point
from mosaic_algorithms.auxiliar_functions.grid_functions.floodFillAlgorithm import floodFillAlgorithm


def grid2D(fpref, olapx, olapy, gamma_, targetArea):
    """
    Grid discretization (using flood-fill algorithm) of a region of interest
    given a reference footprint (unit measure to create the allocatable cells)

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         10/2022

    Usage:        matrixGrid = grid2D(fpref, ovlapx, ovlapy, gamma, targetArea)

    Inputs:
    > fpref:        dict containing the parameters that define the footprint.
                    In this function, the following are needed:
        # sizex:     footprint size in the x direction (longitude), in deg
        # sizey:     footprint size in the y direction (latitude), in deg
        # angle:     footprint's orientation angle, in deg
                    See function 'footprint' for further information.
    > olapx:        grid footprint overlap in the x direction, in percentage
    > olapy:        grid footprint overlap in the y direction, in percentage
    > gamma_:        seed ([lon, lat]) that initiates the grid flood-fill, in deg
    > targetArea:   matrix containing the vertices of the ROI polygon.
                    The vertex points are expressed in 2D.
        # targetArea[:,0] correspond to the x values of the vertices
        # targetArea[:,1] correspond to the y values of the vertices

    Outputs:
    > matrixGrid:   list of lists containing the grid discretization of the
                    region-of-interest (ROI).
                    Each point is defined by the instrument boresight
                    projection onto the body surface, in latitudinal
                    coordinates [lon lat], in deg
    The matrix sorts the discretized points (flood-fill) by latitude and
    longitude according to the following structure:

                    longitude
                    (-) --------> (+)
      latitude (+) [a11]  [a12] ⋯
                ¦  [a21]
                ¦    ⋮
                ∨
               (-)
    > dirx:         unit vector representing the direction of the x-axis in the grid
    > diry:         unit vector representing the direction of the y-axis in the grid
    """
    gamma = copy.deepcopy(gamma_)
    # Pre-allocate variables
    matrixGrid = []

    # Get the footprint angle, i.e., the angle that the 2D footprint forms with
    # respect to the meridian-equator axes
    angle = np.deg2rad(-fpref['angle'])

    # Filling the region-of-interest (roi) with a footprint that is not aligned
    # with the meridian-equator axes is equivalent to filling the oriented
    # target area with an aligned footprint (angle = 0). Therefore, we rotate
    # the region-of-interest to orient it according to the footprint
    rotmat = np.array([[np.cos(angle), -np.sin(angle)],
                       [np.sin(angle), np.cos(angle)]])

    # matrixGrid directions x and y
    dirx = rotmat[0, :]
    diry = rotmat[1, :]

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
        polygon = MultiPolygon(polygon_list)
    else:
        polygon = Polygon(targetArea)


    cx, cy = polygon.centroid.x, polygon.centroid.y

    orientedArea = np.zeros([max(np.shape(targetArea)),2])
    for j in range(max(np.shape(targetArea))):
        orientedArea[j, :] = np.array([cx, cy]) + rotmat @ (targetArea[j, :] - np.array([cx, cy]))

    gamma = np.array([cx, cy]) + rotmat @ (np.array(gamma) - np.array([cx, cy]))

    # If the area is divided in smaller regions, then we get the convex polygon
    # that encloses all of them (flood-fill)
    aux=np.column_stack((orientedArea[~np.isnan(orientedArea[:, 0]), 0] ,orientedArea[~np.isnan(orientedArea[:, 1]), 1]))
    # convhull does not accept NaN nor Inf

    k = (ConvexHull(aux)).vertices # boundary vertices that constitute the convex polygon
    periArea = aux[k]

    # Auxiliary figure
    # plt.figure()
    # plt.plot(*Polygon(orientedArea).exterior.xy)
    # plt.plot(gamma[0], gamma[1], 'r*')
    # plt.plot(*Polygon(periArea).exterior.xy)
    # plt.fill(*Polygon(periArea).exterior.xy, 'none')
    # plt.show()

    # Flood-fill algorithm to get the grid points of the oriented roi
    # gridPoints = floodFillAlgorithm(fpref['sizex'], fpref['sizey'], ovlapx, ovlapy, gamma, orientedArea, gridPoints,np.array([]),
    #                               np.array([]),'8fill')
    gridPoints,_ = floodFillAlgorithm(fpref['width'], fpref['height'], olapx, olapy, gamma, orientedArea, periArea,
                                    np.array([]), np.array([]),'4fill')


    if gridPoints.size != 0:

        # Auxiliary figure
        # plt.figure()
        # plt.plot(*Polygon(targetArea).exterior.xy)
        # plt.hold(True)
        # plt.plot(*Polygon(orientedArea).exterior.xy)
        # plt.plot(gridPoints[:, 0], gridPoints[:, 1], 'b*')
        # orientedGridPoints = np.zeros((len(gridPoints), 2))
        # for j in range(len(gridPoints)):
        #     orientedGridPoints[j, :] = np.array([cx, cy]) + np.dot(rotmat.T, gridPoints[j, :] - np.array([cx, cy]))
        # plt.plot(orientedGridPoints[:, 0], orientedGridPoints[:, 1], 'r*')
        # plt.hold(True)
        # plt.plot(gamma[0], gamma[1], 'g*')
        # plt.show()

        # Sort grid points
        sortedGrid = np.array(sorted(gridPoints, key=lambda x: -x[1])) # the elements of gridPoints are sorted by latitude (+ to -)
        uniqueLat = np.unique([pt[1] for pt in sortedGrid]) # get the different latitude values
        ind = np.abs(np.diff(uniqueLat)) < 1e-5 # double check that there are no "similar" latitude values (it may happen)
        ind = np.append(ind, False)
        uniqueLat = uniqueLat[~ind]
        uniqueLon = np.unique([pt[0] for pt in sortedGrid]) # get the different longitude values unique check
        ind = np.abs(np.diff(uniqueLon)) < 1e-5  # double check that there are no "similar" longitude values (it may happen)
        ind = np.append(ind, False)
        uniqueLon = uniqueLon[~ind]

        # Sort and rotate the grid points and insert them in the grid matrix
        matrixGrid = [[None for _ in range(max(np.shape(uniqueLon)))] for _ in range(max(np.shape(uniqueLat)))]
        for i in range(max(np.shape(uniqueLat))):
            # We will sweep across the grid by, first, latitude and, second, longitude
            lat = uniqueLat[max(np.shape(uniqueLat)) - 1 - i]
            indlat = np.abs(sortedGrid[:, 1] - lat) < 1e-5
            mrow = sortedGrid[indlat]
            mrow = np.sort(mrow[:,0], axis=0)
            for j in range(max(np.shape(mrow))):
                indlon=np.abs(uniqueLon-mrow[j]) < 1e-5
                lon = mrow[j]
                for k in range(max(np.shape(uniqueLon))):
                    if indlon[k]:
                        matrixGrid[i][k] = np.array([cx, cy]) + rotmat.T @ (np.array([lon, lat]) - np.array([cx, cy]))

    return matrixGrid, dirx, diry

