import numpy as np
from shapely.geometry import MultiPolygon, Polygon
from conversion_functions import *
from conversion_functions import mat2py_getfov
from mosaic_algorithms.auxiliar_functions.grid_functions.boustrophedon import boustrophedon
from mosaic_algorithms.auxiliar_functions.grid_functions.grid2D import grid2D
from mosaic_algorithms.auxiliar_functions.grid_functions.inst2topo import inst2topo
from mosaic_algorithms.auxiliar_functions.grid_functions.topo2inst import topo2inst
from mosaic_algorithms.auxiliar_functions.polygon_functions.closestSide import closestSide
from mosaic_algorithms.auxiliar_functions.polygon_functions.minimumWidthDirection import minimumWidthDirection
from mosaic_algorithms.auxiliar_functions.plot.groundtrack import groundtrack


def planSidewinderTour(target, roi, sc, inst, inittime, olapx, olapy):
    """
    This function plans an observation tour using a modified Boustrophedon
    decomposition method. It calculates an optimal path for observing a ROI
    on a target body, considering the spacecraft's starting position.
    We project the ROI's polygon onto the instrument's focal plane.
    Here, we design the grid based on the image plane's reference. The image
    plane is built according to the FOV's parameters, retrieved from the
    instrument kernel. We also devise the traversal in the image plane (tour),
    and then transform it back to the topographical coordinates.

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022

    Usage:        topo_tour, inst_grid, inst_tour, grid_dirx, grid_diry, sweepDir1, sweepDir2 = ...
                    planSidewinderTour(target, roi, sc, inst, inittime, olapx, olapy, angle)
    Inputs:
       > target:       string name of the target body
       > roi:          matrix containing the vertices of the uncovered area
                       of the ROI polygon. The vertex points are expressed in
                       2D, in latitudinal coordinates [º]
       > sc:           string name of the spacecraft
       > inst:         string name of the instrument
       > inittime:     start time of the planning horizon, in TDB seconds past
                       J2000 epoch
       > olapx:        grid footprint overlap in the x direction (longitude),
                       in percentage (width)
       > olapy:        grid footprint overlap in the y direction (latitude),
                       in percentage (height)

    Returns:
       > topo_tour:    tour path in topographical coordinates (lat/lon on the
                       target body), in [deg]
       > inst_grid:    grid of potential observation points in instrument
                       frame coordinates
       > inst_tour:    tour path in instrument frame coordinates
       > grid_dirx:    direction of the grid along the x-axis in the
                       instrument frame
       > grid_diry:    direction of the grid along the y-axis in the
                       instrument frame
       > sweepDir1, sweepDir2: directions defining the sweep of the
                               Boustrophedon decomposition, derived according
                               to the spacecraft's position with respect to
                               the ROI

    """

    # Pre-allocate variables
    origin = np.array([0., 0.]) # initialize grid origin for grid generation
    x, y = roi[:, 0], roi[:, 1]

    if (np.isnan(x)).any():
        nanindex = np.where(np.isnan(x))[0]
        polygon_list = []
        for i in range(len(nanindex)):
            if i == 0:
                polygon_list.append(Polygon(list(zip(x[:nanindex[0]], y[:nanindex[0]]))))
            else:
                polygon_list.append(Polygon(
                    list(zip(x[nanindex[i - 1] + 1:nanindex[i]], y[nanindex[i - 1] + 1:nanindex[i]]))))
        if ~ np.isnan(x[-1]):
            polygon_list.append(Polygon(list(zip(x[nanindex[-1] + 1:], y[nanindex[-1] + 1:]))))
        polygon = MultiPolygon(polygon_list)
    else:
        polygon = Polygon((list(zip(x, y))))
    polygon = polygon.buffer(0)
    # point camera at ROI's centroid
    cx = polygon.centroid.x
    cy = polygon.centroid.y



    # Project ROI to the instrument plane
    targetArea = topo2inst(roi, cx, cy, target, sc, inst, inittime)
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
        poly_aux = MultiPolygon(polygon_list)
    else:
        poly_aux =Polygon((list(zip(targetArea[:, 0], targetArea[:, 1]))))


    #poly_aux = poly_aux.buffer(0)

    origin[0], origin[1] = poly_aux.centroid.x, poly_aux.centroid.y


    # Get minimum width direction of the footprint
    angle,_,_,_ = minimumWidthDirection(targetArea[:,0],targetArea[:,1])
    # observation angle, influencing the orientation of the observation footprints and,
    # therefore, the coverage path orientation

    # Retrieve the field of view (FOV) bounds of the instrument and calculate
    # the dimensions of a reference observation footprint
    _, _, _, bounds = mat2py_getfov(mat2py_bodn2c(inst)[0], 4) # get fovbounds in the instrument's reference frame
    maxx, minx = np.max(bounds[0, :]), np.min(bounds[0, :])
    maxy, miny = np.max(bounds[1, :]), np.min(bounds[1, :])

    # Define the footprint reference dimensions and orientation
    fpref = {
        'width': maxx - minx,
        'height': maxy - miny,
        'angle': angle
    }
    # we enforce the orientation angle of the footprint
    # to be that of its projection onto the topographical grid with respect to
    # the reference axes (east-north), given as an input. grid2D will use this
    # angle to orient the ROI according to this orientation
    # [Future work]: This angle could be given as an input to grid2D?

    # [Check]: there is no need to compute the angle because the targetArea
    # projection already accounts for that
    if fpref['width'] <= fpref['height']:
        angle = 0.
    else:
        angle = 0.
    fpref['angle'] = angle

    gt1 = np.array([0.,0.])
    gt2 = np.array([0.,0.])
    # Closest polygon side to the spacecraft's ground track position (this
    # will determine the coverage path)
    gt1[0],gt1[1] = groundtrack(sc, inittime, target) # initial ground track position
    gt2[0],gt2[1] = groundtrack(sc, inittime + 500, target) # future ground track position
    gt1 = topo2inst(gt1, cx, cy, target, sc, inst, inittime) # projected initial position
    gt2 = topo2inst(gt2, cx, cy, target, sc, inst, inittime + 500) # projected future position

    # Calculate the closest side of the target area to the spacecraft's ground track,
    # determining the observation sweep direction
    sweepDir1, sweepDir2 = closestSide(gt1, gt2, targetArea, angle)

    # Focal plane grid discretization based on the reference footprint (FOV
    # plane) and specified overlap
    inst_grid, grid_dirx, grid_diry = grid2D(fpref, olapx, olapy, origin, targetArea)

    # Boustrophedon decomposition to generate grid traversal
    inst_tour = boustrophedon(inst_grid, sweepDir1, sweepDir2)


    # Convert grid and tour from instrument frame to topographical coordinates
    topo_tour = inst2topo([inst_tour], cx, cy, target, sc, inst, inittime)[0]

    # Remove empty elements from the tour, which may result from unobservable regions
    # within the planned path 
    topo_tour = [point for point in topo_tour if point is not None]

    return topo_tour, inst_grid, inst_tour, grid_dirx, grid_diry, sweepDir1, sweepDir2


