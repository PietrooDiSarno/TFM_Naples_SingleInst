import copy

import numpy as np
from shapely.geometry import MultiPolygon, Polygon, Point



def floodFillAlgorithm(w, h, olapx, olapy, gamma, targetArea, perimeterArea, gridPoints_, vPoints_, method):
    """
    Flood-fill recursive algorithm that discretizes the target area by
    "flooding" the region with 2D rectangular elements. The grid is determined
    by the input width, height and overlaps in both directions.

    Programmers: Paula Betriu (UPC/ESEIAAT)
    Date:        10/2022

    Usage:        gridPoints,vPoints = floodFillAlgorithm(w,h,ovlapx,ovlapy,gamma,targetArea,
                                                          gridPoints,method)

    Inputs:
        - w:            horizontal resolution. Units are irrelevant as long as they are consistent.
        - h:            vertical resolution. Units are irrelevant as long as they are consistent.
        - olapx:        grid footprint overlap in the horizontal direction. Units are in percentage of width.
        - olapy:        grid footprint overlap in the vertical direction. Units are in percentage of the height.
        - gamma:        grid origin point (seed)
        - targetArea:   matrix containing the vertices of the ROI polygon. The vertex points are expressed in 2D.
            # targetArea(:,1) correspond to the x values of the vertices
            # targetArea(:,2) correspond to the y values of the vertices
        - perimeterArea:matrix containing the vertices of the polygon that encloses all of the uncovered
                        area. At the beginning: perimeterArea = targetArea, but as the observations
                        advance, this is going to change. Recommended: use the convex hull function of the uncovered
                        area.
        - gridPoints_:   matrix containing the discretized grid points of the region-of-interest. The recursive calls of the
                        algorithm will fill this matrix. These grid points represent the center of the rectangular elements
                        used to fill the region.
            # When calling this function: gridPoints = np.array([])
        - vPoints_:      matrix containing the visited points (to prevent gridlock)
            # When calling this function: vPoints_ = np.array([])
        - method:       string name of the method. '4fill' fills the roi by
                        searching the cardinal directions.'8fill' considers
                        also the diagonal neighbors.


    Returns:
        - gridPoints:   matrix containing the discretized gridPoints of the
                        region-of-interest. The recursive calls of the
                        algorithm will fill this matrix. These grid points
                        represent the center of the rectangular elements used
                        to fill the region
        - vPoints:      matrix containing the visited points (to prevent
                        gridlock)

    Note: this function creates a convex polygon that encloses all the
    uncovered area (even when this is divided in portions) and tours the
    whole area. It is less computationally efficient than the classic
    flood-fill algorithm, but it is convenient to prevent sub-optimal
    fillings of the uncovered area (isolated points).
    """
    gridPoints = copy.deepcopy(gridPoints_)
    vPoints = copy.deepcopy(vPoints_)

    if isinstance(gridPoints,np.ndarray):
        gridPoints = list(gridPoints)
    if isinstance(vPoints,np.ndarray):
        vPoints = list(vPoints)

    # Variables pre-allocation
    inside= False
    ovlapx = olapx*w/100; ovlapy = olapy*h/100 # convert overlaps from
    # percentage to degrees of latitude and longitude, respectively
    epsilon = 0.05

    # Check if the cell has been previously visited
    for vp in vPoints:
        if np.linalg.norm(np.array(vp) - np.array(gamma)) < 1e-5:
            return np.array(gridPoints), np.array(vPoints)

    # Otherwise, mark this point as visited
    vPoints.append(np.array(gamma))

    # Rectangular element definition
    fpx = [gamma[0] - w / 2, gamma[0] - w / 2, gamma[0] + w / 2, gamma[0] + w / 2]
    fpy = [gamma[1] + h / 2, gamma[1] - h / 2, gamma[1] - h / 2, gamma[1] + h / 2]


    # Subtract the allocated cell (footprint) from the perimeterArea
    if (np.isnan(perimeterArea[:, 0])).any():
        nanindex = np.where(np.isnan(perimeterArea[:, 0]))[0]
        polygon_list = []
        for i in range(len(nanindex)):
            if i == 0:
                polygon_list.append(Polygon(list(zip(perimeterArea[:nanindex[0], 0], perimeterArea[:nanindex[0], 1]))))
            else:
                polygon_list.append(Polygon(
                    list(zip(perimeterArea[nanindex[i - 1] + 1:nanindex[i], 0], perimeterArea[nanindex[i - 1] + 1:nanindex[i], 1]))))
        if ~ np.isnan(perimeterArea[-1, 0]):
            polygon_list.append(Polygon(list(zip(perimeterArea[nanindex[-1] + 1:, 0], perimeterArea[nanindex[-1] + 1:, 1]))))
        peripshape = MultiPolygon(polygon_list)
    else:
        peripshape = Polygon(perimeterArea)

    fpshape = Polygon(zip(fpx, fpy))
    inter = (peripshape.difference(fpshape)).buffer(0)
    areaI = inter.area
    areaP = peripshape.area

    # Check: the footprint is larger than the region of interest...
    if areaI == 0:
        gridPoints.append(np.array(gamma))
        return np.array(gridPoints), np.array(vPoints)

    # Check if the rectangle at gamma and size [w,h] is contained in
    # the perimeter area (either partially or totally)
    if (np.isnan(targetArea[:, 0])).any():
        nanindex = np.where(np.isnan(targetArea[:, 0]))[0]
        polygon_list = []
        for i in range(len(nanindex)):
            if i == 0:
                polygon_list.append(Polygon(list(zip(targetArea[:nanindex[0], 0], targetArea[:nanindex[0], 1]))))
            else:
                polygon_list.append(Polygon(
                    list(zip(targetArea[nanindex[i - 1] + 1:nanindex[i], 0],
                             targetArea[nanindex[i - 1] + 1:nanindex[i], 1]))))
        if ~ np.isnan(targetArea[-1, 0]):
            polygon_list.append(Polygon(list(zip(targetArea[nanindex[-1] + 1:, 0], targetArea[nanindex[-1] + 1:, 1]))))
        target_polygon = MultiPolygon(polygon_list)
    else:
        target_polygon = Polygon(targetArea)
    target_polygon = target_polygon.buffer(0)
    if target_polygon.intersects(Point(gamma)) or abs(areaI - areaP) / fpshape.area > epsilon:
        inside = True


    if inside:
        # Disregard those cases where the footprint does not cover a certain
        # minimum of the roi (this also avoids sub-optimality in the
        # optimization algorithms)
        areaT = target_polygon.area
        inter = (target_polygon.difference(fpshape)).buffer(0)
        areaI = inter.area
        areaInter = areaT - areaI
        fpArea = fpshape.area


        if areaInter / fpArea > epsilon:
            gridPoints.append(np.array(gamma))
            # coordinates = [(gamma(0)-w/2, gamma(1)+ h/2),
            #                (gamma(0)-w/2, gamma(1)- h/2),
            #                (gamma(0)+w/2, gamma(1)- h/2),
            #                (gamma(0)+w/2, gamma(1)+ h/2)]
            # polygon = Polygon(coordinates)
            # plt.fill(*polygon.exterior.xy, facecolor='orange', alpha=0.2)
            # plt.plot(gamma(1), gamma(2), 'r^')
            # plt.show()

        else:
            if not gridPoints:
                return np.array(gridPoints), np.array(vPoints)
                # plot(gamma(0),gamma(1),'b^')
                # plt.show()
        # Check the cardinal (and diagonal neighbors in case the method is set
        # to 8fill) neighbors recursively
        # West
        gridPoints, vPoints = floodFillAlgorithm(w, h, olapx, olapy,
                                                     np.array([gamma[0] - w + ovlapx, gamma[1]]),
                                                     targetArea, perimeterArea,
                                                     gridPoints, vPoints, method)
        # South
        gridPoints, vPoints = floodFillAlgorithm(w, h, olapx, olapy,
                                                     np.array([gamma[0], gamma[1] - h + ovlapy]),
                                                     targetArea, perimeterArea,
                                                     gridPoints, vPoints, method)
        # North
        gridPoints, vPoints = floodFillAlgorithm(w, h, olapx, olapy,
                                                     np.array([gamma[0], gamma[1] + h - ovlapy]),
                                                     targetArea, perimeterArea,
                                                     gridPoints, vPoints, method)
        # East
        gridPoints, vPoints = floodFillAlgorithm(w, h, olapx, olapy,
                                                     np.array([gamma[0] + w - ovlapx, gamma[1]]),
                                                     targetArea, perimeterArea,
                                                     gridPoints, vPoints, method)

        # If method is '8fill', check diagonal neighbors
        if method == '8fill':
            # Northwest
            gridPoints, vPoints = floodFillAlgorithm(w, h, olapx, olapy,
                                                         np.array([gamma[0] - w + ovlapx, gamma[1] + h - ovlapy]),
                                                         targetArea, perimeterArea,
                                                         gridPoints, vPoints, method)
            # Southwest
            gridPoints, vPoints = floodFillAlgorithm(w, h, olapx, olapy,
                                                         np.array([gamma[0] - w + ovlapx, gamma[1] - h + ovlapy]),
                                                         targetArea, perimeterArea,
                                                         gridPoints, vPoints, method)
            # Northeast
            gridPoints, vPoints = floodFillAlgorithm(w, h, olapx, olapy,
                                                         np.array([gamma[0] + w - ovlapx, gamma[1] + h - ovlapy]),
                                                         targetArea, perimeterArea,
                                                         gridPoints, vPoints, method)
            # Southeast
            gridPoints, vPoints = floodFillAlgorithm(w, h, olapx, olapy,
                                                         np.array([gamma[0] + w - ovlapx, gamma[1] - h + ovlapy]),
                                                         targetArea, perimeterArea,
                                                         gridPoints, vPoints, method)

    return np.array(gridPoints), np.array(vPoints)