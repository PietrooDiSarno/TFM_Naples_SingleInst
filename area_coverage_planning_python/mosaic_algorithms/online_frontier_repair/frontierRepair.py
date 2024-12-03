import numpy as np
import copy
from shapely.geometry import MultiPolygon, Polygon

from conversion_functions import mat2py_et2utc
from mosaic_algorithms.auxiliar_functions.polygon_functions.visibleroi import visibleroi
from mosaic_algorithms.auxiliar_functions.polygon_functions.interppolygon import interppolygon
from mosaic_algorithms.sidewinder.planSidewinderTour import planSidewinderTour
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.footprint import footprint
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.processObservation import processObservation
from mosaic_algorithms.auxiliar_functions.polygon_functions.sortcw import sortcw
from mosaic_algorithms.online_frontier_repair.updateGrid import updateGrid

def frontierRepair(startTime, endTime, tobs, inst, sc, target, inroi, olapx, olapy, slewRate, *args):
    """
    This function adjusts the observation grid and planning tour in response
    to new observations. It takes into account changes in observation
    geometry, incorporating new observation points, and removes points that
    are no longer necessary for covering the region of interest (ROI).
    A Boustrophedon pattern is used to ensure efficient coverage. It has been
    adapted from [1]

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022

    Usage:        [A, fpList] = frontierRepair(startTime, endTime, ...
                    tobs, inst, sc, target, inroi, olapx, olapy, slewRate, ...
                    resolution)

    Inputs:
      > startTime:    start time of the planning horizon, in TDB seconds past
                      J2000 epoch
      > endTime:      end time of the planning horizon, in TBD seconds past
                      J2000 epoch
      > tobs:         observation time, i.e. the minimum time that the
                      instrument needs to perform an observation, in seconds
      > inst:         string name of the instrument
      > sc:           string name of the spacecraft
      > target:       string name of the target body
      > inroi:        matrix containing the vertices of the ROI polygon. The
                      vertex points are expressed in 2D, in latitudinal
                      coordinates [º]
          # roi[:,0] correspond to the x values of the vertices
          # roi[:,1] correspond to the y values of the vertices
      > olapx:        grid footprint overlap in the x direction (longitude),
                      in percentage (width)
      > olapy:        grid footprint overlap in the y direction (latitude),
                      in percentage (height)
      > slewRate:     rate at which the spacecraft (or instrument platform)
                      can slew between observations, in [º/s]
      > resolution:   string 'lowres' or 'highres' that determines the
                      footprint resolution calculation. It's set to 'lowres'
                      by default

    Outputs:
      > A:            cell matrix of successive instrument observations,
                      sorted in chronological order.
                      Each observation is defined by the instrument boresight
                      projection onto the body surface, in latitudinal
                      coordinates [lon lat], in deg
      > fpList:       list of footprint structures detailing the observation
                      metadata and coverage

    [1] Shao, E., Byon, A., Davies, C., Davis, E., Knight, R., Lewellen, G.,
    Trowbridge, M. and Chien, S. (2018). Area coverage planning with 3-axis
    steerable, 2D framing sensors.
    """
    # Pre-allocate variables
    A = []  # List of observations (successive boresight ground track position)
    fpList = []
    amIntercept = False
    if len(args) == 1:
        resolution = args[0]
    else:
        resolution = 'lowres'

    # Check ROI visible area from spacecraft
    vsbroi, _, visibilityFlag = visibleroi(inroi, startTime, target, sc)  # polygon vertices of the visible area
    if visibilityFlag:
        print("ROI is not visible from the instrument")
        return A, fpList

    roi = interppolygon(vsbroi)  # interpolate polygon vertices (for improved accuracy)

    # Previous anti-meridian intersection check...
    ind = np.where(np.diff(np.sort(inroi[:, 0])) >= 180)[0]  # find the discontinuity index
    if ind.size > 0:
        amIntercept = True
        roi = copy.deepcopy(inroi)
        roi[roi[:, 0] < 0, 0] += 360  # adjust longitudes
        roi[:,0],roi[:,1] = sortcw(roi[:,0],roi[:,1])  # sort coordinates clockwise

    # [Issue]: We cannot perform visibility and anti-meridian checks
    # simultaneously. This means, either we get a sectioned ROI due to
    # visibility or anti-meridian, but both may not happen. This is because the
    # interpolation function does not work with anti-meridian intercepts.
    # [Future work]: Solve this incompatibility

    # Define target area as a polygon
    if (np.isnan(roi[:,0])).any():
        nanindex = np.where(np.isnan(roi[:,0]))[0]
        polygon_list = []
        for i in range(len(nanindex)):
            if i==0:
                polygon_list.append(Polygon(list(zip(roi[:nanindex[0],0], roi[:nanindex[0],1]))))
            else:
                polygon_list.append(Polygon(list(zip(roi[nanindex[i-1]+1:nanindex[i],0], roi[nanindex[i-1]+1:nanindex[i],1]))))
        if ~ np.isnan(roi[-1,0]):
            polygon_list.append(Polygon(list(zip(roi[nanindex[-1] + 1:, 0], roi[nanindex[-1] + 1:, 1]))))
        poly1 = MultiPolygon(polygon_list)
    else:
        poly1 = Polygon((list(zip(roi[:, 0], roi[:, 1]))))

    poly1 = poly1.buffer(0)

    cx = poly1.centroid.x
    cy = poly1.centroid.y

    ## Frontier Repair algorithm
    # The first time iteration is the starting time in the planning horizon
    t = startTime

    # Boolean that defines when to stop covering the target area
    exit = False

    while not exit:
        # Initial 2D grid layout discretization: the instrument's FOV is going
        # to be projected onto the uncovered area's centroid and the resulting
        # footprint shape is used to set the grid spatial resolution

        if (np.isnan(roi[:, 0])).any():
            nanindex = np.where(np.isnan(roi[:, 0]))[0]
            polygon_list = []
            for i in range(len(nanindex)):
                if i == 0:
                    polygon_list.append(Polygon(list(zip(roi[:nanindex[0], 0], roi[:nanindex[0], 1]))))
                else:
                    polygon_list.append(Polygon(
                        list(zip(roi[nanindex[i - 1] + 1:nanindex[i], 0], roi[nanindex[i - 1] + 1:nanindex[i], 1]))))
            if ~ np.isnan(roi[-1, 0]):
                polygon_list.append(Polygon(list(zip(roi[nanindex[-1] + 1:, 0], roi[nanindex[-1] + 1:, 1]))))
            polyroi = MultiPolygon(polygon_list)
        else:
            polyroi = Polygon((list(zip(roi[:, 0], roi[:, 1]))))

        polyroi = polyroi.buffer(0)
        gamma = [polyroi.centroid.x,polyroi.centroid.y]
        fprintc = footprint(t, inst, sc, target, resolution, gamma[0], gamma[1], 0)  # centroid footprint

        # Initialize a list of dictionaries to save footprints
        if t == startTime:
            #fpList = [{} for _ in range(len(fprintc))]
            fpList = [{}]
            for fn in fprintc.keys():
                fpList[0][fn] = []

        # Check roi visibility
        vsbroi, _, visibilityFlag = visibleroi(roi, t, target, sc)
        if visibilityFlag:
            print("ROI no longer reachable")
            break
        else:
            roi = interppolygon(vsbroi)

        # Discretize ROI area (grid) and plan Sidewinder tour based on a Boustrophedon approach
        tour, grid, itour, grid_dirx, grid_diry, dir1, dir2 = planSidewinderTour(target, roi, sc, inst, t, olapx, olapy)

        #for i in range(len(grid)):
        #    for j in range(len(grid[i])):
        #       if grid[i][j] is not None:
        #            grid[i][j] = (grid[i][j]).reshape(1,2)

        # Handle cases where the FOV projection is larger than the ROI area
        if len(tour) < 1:
            A.append(gamma)
            fpList.append(fprintc)
            exit = True
            continue

        seed = itour[0]
        while len(tour)!= 0 and t < endTime:
            # Update origin and tour
            old_seed = seed
            itour.pop(0)

            # Process each point of the tour
            A, tour, fpList, poly1, t, _  = processObservation(A, tour, fpList, poly1, t, slewRate, tobs, amIntercept, inst,
                                                           sc, target, resolution)
            if isinstance(poly1, Polygon):
                # If polygon is completely covered, break loop
                if not poly1.exterior.coords:
                    break
                # Update roi
                roi = np.array(poly1.exterior.coords)
            elif isinstance(poly1, MultiPolygon):
                for i in range(len(poly1.geoms)):
                    if i == 0:
                        roi = np.vstack((np.array(poly1.geoms[i].exterior.coords), [np.nan, np.nan]))
                    else:
                        roi = np.vstack((roi, np.array(poly1.geoms[i].exterior.coords), [np.nan, np.nan]))
                roi = roi[:-1,:]

            # Check roi visibility
            vsbroi, _, visibilityFlag = visibleroi(roi, t, target, sc)
            if visibilityFlag:
                print("ROI no longer reachable")
                break
            else:
                roi = interppolygon(vsbroi)

            if len(tour) == 0:
                break
            else:
                gamma = tour[0]  # next observation point
                seed = itour[0]  # next seed in the image plane

                # Update previous grid with the new tile reference (footprint),
                # looking for new potential tiles and/or disposable ones
                seed, grid, itour, tour = updateGrid(roi, itour, grid, grid_dirx, grid_diry, cx, cy, olapx, olapy, dir1,
                                                     dir2, seed, old_seed, gamma, t, inst, sc, target)

        # For now, the stop criteria is the end of the tour, re-starts are not
        # optimal for the purposes of the scheduling problem
        # [Future work]: automated scheduling (in-situ). Re-starts may be
        # considered, and we will need to define a criteria to prompt those.
        exit = True

    # OK message
    print('Online Frontier successfully executed')

    # Remove first element of fplist (it was just to set the struct fields)
    if len(fpList) > 0:
        fpList.pop(0)

    return A, fpList
