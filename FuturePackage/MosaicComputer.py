# Llamar al init desde el main, en checkOneROI se llamará a roi y en initObservationDataBase se llamará a precomputeOnlineFrontier, tendré que pasarle los datos que no tiene del init.
# Valorar si es mejor dejarla como una clase standalone o meterla dentro de instrument.
# AÑADIR SLEW RATE A INSTRUMENT
# AÑADIR NAME A INSTRUMENT
# mirar tobs (en runFrontierRepair) y z, no están definidas (creo)
# Mirar taboo, sweepdir, etc. En el original son variables globales. Mirar cómo afecta al código. Quizás hay que reorganizar variables.

#Cambios original vs clase:
#roi-> roi.vertices
#target -> roi.body/instrument.target  si lo meto dentro de instrument
#inst -> instrument.name
#target -> instrument.mission
#timeint -> timeinterval
import warnings
import copy
import numpy as np
import math
from geopy.distance import great_circle
from geopy.point import Point
from shapely.geometry import Polygon, MultiPolygon, LineString
from scipy.spatial import ConvexHull
import area_coverage_planning_python.conversion_functions as M2P
import PSOA as psoa

class MosaicComputer:
    def __init__(self, tcadence, olapx, olapy):
        self.tcadence = tcadence #[s] between observations
        self.olapx = olapx  # [%] of overlap in x direction
        self.olapy = olapy  # [%] of overlap in y direction
        self.fpref = None
        self.pointing0 = None
        self.sweepDir1 = None
        self.sweepDir2 = None
        self.taboo = dict()

    def precomputeOnlineFrontier(self, timeinterval, instrument, roi):
        makespan = []
        nImg = []
        resROI = []
        cov = []
        stoptime = timeinterval[-1] + 3600
        for inittime in timeinterval:
            A, fpList = self.runFrontierRepair(inittime, stoptime, instrument, roi)
            if not fpList == []:
                makespan.append(fpList[-1]['t'] + self.tcadence - inittime)
                nImg.append(len(fpList))
                resROI.append(self.computeResMosaic(fpList, instrument.ifov))
                roi_ = {'vertices': roi}
                cov_ = self.roicoverage(roi_, fpList)[0]
                if cov_ > 95:
                    cov.append(100)
                else:
                    cov.append(cov_)
            else:
                makespan.append(None)
                nImg.append(None)
                resROI.append(None)
                cov.append(None)

        return makespan, nImg, resROI, cov

    def runFrontierRepair(self, startTime, endTime, instrument, roiClass, *args):
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
          > instrument:   Instrument class containing:
            > .mission:    string name of the spacecraft
            > .name        string name of the instrument
          > roiClass:          ROI Class containing:
            > .body       string name of the target body
            > .vertices   matrix containing the vertices of the ROI polygon. The
                          vertex points are expressed in 2D, in latitudinal
                          coordinates [º]
              # roi.vertices[:,0] correspond to the x values of the vertices
              # roi.vertices[:,1] correspond to the y values of the vertices
          > self.olapx:        grid footprint overlap in the x direction (longitude),
                          in percentage (width)
          > self.olapy:        grid footprint overlap in the y direction (latitude),
                          in percentage (height)
          > self.slewRate:     rate at which the spacecraft (or instrument platform)
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
        vsbroi, _, visibilityFlag = self.computeVisibleROI(roiClass.vertices, startTime, roiClass.body, instrument.mission)  # polygon vertices of the visible area
        if visibilityFlag:
            print("ROI is not visible from the instrument")
            return A, fpList

        roiVertices = self.interppolygon(vsbroi)  # interpolate polygon vertices (for improved accuracy)

        # Previous anti-meridian intersection check...
        ind = np.where(np.diff(np.sort(roiClass.vertices[:, 0])) >= 180)[0]  # find the discontinuity index
        if ind.size > 0:
            amIntercept = True
            roiVertices = copy.deepcopy(roiClass.vertices)
            roiVertices[roiVertices[:, 0] < 0, 0] += 360  # adjust longitudes
            roiVertices[:, 0], roiVertices[:, 1] = self.sortcw(roiVertices[:, 0], roiVertices[:, 1])  # sort coordinates clockwise

        # [Issue]: We cannot perform visibility and anti-meridian checks
        # simultaneously. This means, either we get a sectioned ROI due to
        # visibility or anti-meridian, but both may not happen. This is because the
        # interpolation function does not work with anti-meridian intercepts.
        # [Future work]: Solve this incompatibility

        # Define target area as a polygon
        if (np.isnan(roiVertices[:, 0])).any():
            nanindex = np.where(np.isnan(roiVertices[:, 0]))[0]
            polygon_list = []
            for i in range(len(nanindex)):
                if i == 0:
                    polygon_list.append(Polygon(list(zip(roiVertices[:nanindex[0], 0], roiVertices[:nanindex[0], 1]))))
                else:
                    polygon_list.append(Polygon(
                        list(zip(roiVertices[nanindex[i - 1] + 1:nanindex[i], 0], roiVertices[nanindex[i - 1] + 1:nanindex[i], 1]))))
            if ~ np.isnan(roiVertices[-1, 0]):
                polygon_list.append(Polygon(list(zip(roiVertices[nanindex[-1] + 1:, 0], roiVertices[nanindex[-1] + 1:, 1]))))
            poly1 = MultiPolygon(polygon_list)
        else:
            poly1 = Polygon((list(zip(roiVertices[:, 0], roiVertices[:, 1]))))

        poly1 = poly1.buffer(0)

        cx = poly1.centroid.x
        cy = poly1.centroid.y

        ## Frontier Repair algorithm
        # The first time iteration is the starting time in the planning horizon
        t = startTime

        # Boolean that defines when to stop covering the target area
        ext = False

        while not ext:
            # Initial 2D grid layout discretization: the instrument's FOV is going
            # to be projected onto the uncovered area's centroid and the resulting
            # footprint shape is used to set the grid spatial resolution

            if (np.isnan(roiVertices[:, 0])).any():
                nanindex = np.where(np.isnan(roiVertices[:, 0]))[0]
                polygon_list = []
                for i in range(len(nanindex)):
                    if i == 0:
                        polygon_list.append(Polygon(list(zip(roiVertices[:nanindex[0], 0], roiVertices[:nanindex[0], 1]))))
                    else:
                        polygon_list.append(Polygon(
                            list(
                                zip(roiVertices[nanindex[i - 1] + 1:nanindex[i], 0], roiVertices[nanindex[i - 1] + 1:nanindex[i], 1]))))
                if ~ np.isnan(roiVertices[-1, 0]):
                    polygon_list.append(Polygon(list(zip(roiVertices[nanindex[-1] + 1:, 0], roiVertices[nanindex[-1] + 1:, 1]))))
                polyroi = MultiPolygon(polygon_list)
            else:
                polyroi = Polygon((list(zip(roiVertices[:, 0], roiVertices[:, 1]))))

            polyroi = polyroi.buffer(0)
            gamma = [polyroi.centroid.x, polyroi.centroid.y]
            fprintc = self.computeFootprint(t, instrument.name, instrument.mission, roiClass.body, resolution, gamma[0], gamma[1], 0)  # centroid footprint

            # Initialize a list of dictionaries to save footprints
            if t == startTime:
                # fpList = [{} for _ in range(len(fprintc))]
                fpList = [{}]
                for fn in fprintc.keys():
                    fpList[0][fn] = []

            # Check roi visibility
            vsbroi, _, visibilityFlag = self.computeVisibleROI(roiVertices, t, roiClass.body, instrument.mission)
            if visibilityFlag:
                print("ROI no longer reachable")
                break
            else:
                roiVertices = self.interppolygon(vsbroi)

            # Discretize ROI area (grid) and plan Sidewinder tour based on a Boustrophedon approach
            tour, grid, itour, grid_dirx, grid_diry, dir1, dir2 = self.planSidewinderTour(roiClass.body, roiVertices, instrument.mission, instrument.name, t)

            # for i in range(len(grid)):
            #    for j in range(len(grid[i])):
            #       if grid[i][j] is not None:
            #            grid[i][j] = (grid[i][j]).reshape(1,2)

            # Handle cases where the FOV projection is larger than the ROI area
            if len(tour) < 1:
                A.append(gamma)
                fpList.append(fprintc)
                ext = True
                continue

            seed = itour[0]
            while len(tour) != 0 and t < endTime:
                # Update origin and tour
                old_seed = seed
                itour.pop(0)

                # Process each point of the tour
                A, tour, fpList, poly1, t, _ = self.processObservation(A, tour, fpList, poly1, t, tobs,
                                                                  amIntercept, instrument.name,
                                                                  instrument.mission, roiClass.body, resolution)
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
                            roiVertices = np.vstack((roiVertices, np.array(poly1.geoms[i].exterior.coords), [np.nan, np.nan]))
                    roiVertices = roiVertices[:-1, :]

                # Check roi visibility
                vsbroi, _, visibilityFlag = self.computeVisibleROI(roiVertices, t, roiClass.body, instrument.mission)
                if visibilityFlag:
                    print("ROI no longer reachable")
                    break
                else:
                    roiVertices = self.interppolygon(vsbroi)

                if len(tour) == 0:
                    break
                else:
                    gamma = tour[0]  # next observation point
                    seed = itour[0]  # next seed in the image plane

                    # Update previous grid with the new tile reference (footprint),
                    # looking for new potential tiles and/or disposable ones
                    seed, grid, itour, tour = self.updateGrid(roiVertices, itour, grid, grid_dirx, grid_diry, cx, cy, dir1,
                                                         dir2, seed, old_seed, gamma, t, instrument.name, instrument.mission, roiClass.body)

            # For now, the stop criteria is the end of the tour, re-starts are not
            # optimal for the purposes of the scheduling problem
            # [Future work]: automated scheduling (in-situ). Re-starts may be
            # considered, and we will need to define a criteria to prompt those.
            ext = True

        # OK message
        print('Online Frontier successfully executed')

        # Remove first element of fplist (it was just to set the struct fields)
        if len(fpList) > 0:
            fpList.pop(0)

        return A, fpList

    def updateGrid(self,roi, inst_tour, inst_grid, grid_dirx, grid_diry, cx, cy, olapx, olapy,
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

        if self.sweepDir1 is None:
            self.sweepDir1 = insweepDir1
            self.sweepDir2 = insweepDir2

        if self.pointing0 is None:
            self.pointing0 = np.array([cx, cy])

        cin = []  # set of new observations (inside or outside from 'tour')
        cind = []  # map indices of the new observation (potential placement in the map)
        cout = []  # set of disposable observations (inside or outside from 'tour')
        N = []  # set of new tiles (outside from 'tour')
        Nind = []  # map indices of the new tiles
        X = []  # set of disposable tiles (inside of 'tour')
        epsilon = 0.02

        # Build reference tile (it's always going to be the same in subsequent calls)
        if self.fpref is None:
            # We need to craft the tile reference that we are going to use throughout the heuristic operations. To avoid
            # repetitions, we initialize fpref to None and calculate it only on the first call of the function, ensuring
            # that it is set only once.

            _, _, _, bounds = M2P.mat2py_getfov(M2P.mat2py_bodn2c(inst)[0], 4)  # get fovbounds in the instrument's reference frame
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
        targetArea = self.topo2inst(roi, cx, cy, target, sc, inst, et)
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
                polygon_list.append(
                    Polygon(list(zip(targetArea[nanindex[-1] + 1:, 0], targetArea[nanindex[-1] + 1:, 1]))))
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
        updated_seed = (self.topo2inst(np.array([np.array(gamma)]), cx, cy, target, sc, inst, et))

        if not np.isnan(updated_seed).all() and np.size(updated_seed) != 0:
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
        map = self.grid2map(inst_grid)

        # Obtain the frontier tiles in the map: points that have less than 8
        # neighbours in the grid
        frontier, indel = self.getFrontierTiles(map)

        # Update grid
        openList = copy.deepcopy(frontier)  # open list starts as frontier set F
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
            aux = np.array(o) + np.array(self.fpref['bvertices'])
            if (np.isnan(aux[:, 0])).any():
                nanindex = np.where(np.isnan(aux[:, 0]))[0]
                polygon_list = []
                for i in range(len(nanindex)):
                    if i == 0:
                        polygon_list.append(Polygon(list(zip(aux[:nanindex[0], 0], aux[:nanindex[0], 1]))))
                    else:
                        polygon_list.append(Polygon(
                            list(
                                zip(aux[nanindex[i - 1] + 1:nanindex[i], 0], aux[nanindex[i - 1] + 1:nanindex[i], 1]))))
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
                    if np.array_equal(o, inst_tour[i].reshape(1, 2)):
                        insideTour = True
                        break
                if not insideTour:
                    # If it wasn't included, then its membership changed
                    membershipChanged = True
            else:  # Otherwise (the observation's footprint falls outside the ROI)
                cout.append(o)  # Add it to the list of disposal tiles

                # Identify if the observation was already included in the planned tour
                for i in range(len(inst_tour)):
                    if np.array_equal(o, inst_tour[i].reshape(1, 2)):
                        insideTour = True
                        break
                if insideTour:
                    # If it was included, then its membership changed
                    membershipChanged = True

            if inSeed or membershipChanged:
                # Get the observation neighbor elements (diagonal and cardinal)
                n, nind = self.getNeighbours(o, currind, self.fpref['width'], self.fpref['height'], olapx, olapy, grid_dirx,
                                        grid_diry)

                # Check if the neighbors are already inside tour (in that case it is not necessary to include them in
                # the openlist for re-evaluation)
                nindel = []
                for i in range(len(n)):
                    inTour = False
                    for j in range(len(inst_tour)):
                        if np.linalg.norm(n[i] - inst_tour[j]) < 1e-5:
                            inTour = True
                            break
                    if inTour:
                        nindel.append(i)

                n = [n[i] for i in range(len(n)) if i not in nindel]
                nind = [nind[i] for i in range(len(nind)) if i not in nindel]

                for i in range(len(n)):

                    # Check if the neighbors are included in the cin or cout sets
                    in1 = False
                    in2 = False
                    in1 = any(np.linalg.norm(cin_point - n[i]) < 1e-5 for cin_point in cin)
                    in2 = any(np.linalg.norm(cout_point - n[i]) < 1e-5 for cout_point in cout)

                    # If the neighbour node is not in the cin list nor in the cout list... then add it to the openList
                    # for evaluation (if not already included)
                    if not in1 and not in2:
                        if not any(np.linalg.norm(open_point - n[i]) < 1e-5 for open_point in openList):  # if it's not
                            # in openList, then add it
                            openList.append(n[i])
                            indel.append(nind[i])

        # Identify new tiles N = Cin - Tour
        for i, c in enumerate(cin):
            if not any(np.linalg.norm(c - tour_point) < 1e-5 for tour_point in
                       inst_tour):  # if c is checked to be outside, include it in the new tiles set
                N.append(c)
                Nind.append(cind[i])

        # Check that new identified tiles are not taboo  (moving backwards in the coverage path)
        ind_row, ind_col = None, None

        for j in range(len(map[0])):
            for i in range(len(map)):
                if not ((j == 0 and i == 0) or (j == len(map[0]) - 1 and i == len(map) - 1)):
                    if not np.isnan(map[i][j]).all() and np.linalg.norm(map[i][j] - seed) < 1e-3:
                        ind_row, ind_col = i, j
                        break

        # Check that N is not coincident with old_seed...
        for i in range(len(N)):
            if np.linalg.norm(N[i] - old_seed) < 1e-5:
                N.pop(i)
                Nind.pop(i)
                break
        N, Nind = self.checkTaboo(N, Nind, map, ind_row, ind_col, self.sweepDir1, self.sweepDir2)

        # Identify tiles to remove: X = Cout - Tour
        X = []
        for c in cout:
            if any(np.linalg.norm(c - tour_point) < 1e-5 for tour_point in inst_tour):  # if c is checked to be
                # in 'tour', include it in the  disposable tiles set
                X.append(c)

        # Remove disposable tiles
        map = self.removeTiles(map, X)

        # Insert new tiles
        map = self.insertTiles(map, N, Nind)

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
        inst_grid = self.map2grid(map)
        inst_tour = self.boustrophedon(inst_grid, self.sweepDir1, self.sweepDir2)

        if inst_tour:
            topo_tour = self.inst2topo([inst_tour], cx, cy, target, sc, inst, et)[0]
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
            inst_grid = self.map2grid(map)
            inst_tour = self.boustrophedon(inst_grid, self.sweepDir1, self.sweepDir2)
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

    def checkTaboo(self,N, Nind, map, ind_row, ind_col, indir1, indir2):
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
        if 'pdir1' not in self.taboo:
            self.taboo['pdir1'] = indir1
            self.taboo['pdir2'] = indir2

        pdir1 = self.taboo['pdir1']
        pdir2 = self.taboo['pdir2']
        ## Previous checks...
        # if not N or not Nind:
        # return

        # Pre-allocate variables
        nindel = np.array([])
        grid = self.map2grid(map)
        dir1, dir2 = self.boustrophedonMod(grid, indir1, indir2)
        dir_change = False
        if pdir2 != indir2:
            dir_change = True

        # Previous checks...
        if not N or not Nind:
            return N, Nind

        # Define the grid boundaries (ending rows and columns in the map)
        for i in range(len(map) - 1, -1, -1):
            el = next((j for j, val in enumerate(map[i]) if not (np.isnan(val)).any()),
                      None)  # get non-NaN elements in the map
            if el is not None:
                Nrow = i
                break

        for i in range(len(map[0]) - 1, -1, -1):
            aux = []
            for k in range(len(map)):
                aux.append(map[k][i])
            el = next((j for j, val in enumerate(aux) if not (np.isnan(val)).any()), None)
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

        self.taboo['pdir1'] = pdir1
        self.taboo['pdir2'] = pdir2

        return N, Nind
    def getFrontierTiles(self,map):
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
                    n = self.getMapNeighbours(i, j, map)

                    # If the observation point has less than 8 planned neighbors, then
                    # it is regarded as a frontier tile
                    if len(n) < 8:
                        frontier.append(tile)
                        indel.append([i, j])

        return frontier, indel
    def processObservation(self,A, tour, fpList, poly1, t, slewRate, tobs, amIntercept, inst, sc, target, resolution):
        """
        This function handles the processing of an observation point by computing
        its footprint, updating the list of completed observations and adjusting
        the remaining area to be covered. It also accounts for anti-meridian
        interception and updates the time for the next observation based on slew
        rate and observation time.

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         09/2023

        Usage:        A, tour, fpList, poly1,t = processObservation(A, tour,
                      fpList, poly1, t, slewRate, tobs, amIntercept, inst sc, target,
                      resolution)

        Inputs:
          > A:           list of lists of successive instrument observations,
                         sorted in chronological order.
                         Each observation is defined by the instrument boresight
                         projection onto the body surface, in latitudinal
                         coordinates [lon lat], in deg
         > tour:         array of remaining observation points in the tour, in
                         latitudinal coordinates [º]
         > fpList:       list of footprint structures detailing the observation
                         metadata and coverage
         > poly1:        current polygon shape of the uncovered area on the
                         target body
         > t:            current time in ephemeris seconds past J2000 epoch
         > slewRate:     rate at which the spacecraft (or instrument platform)
                         can slew between observations, in [º/s]
         > tobs:         observation time, i.e. the minimum time that the
                         instrument needs to perform an observation, in seconds
         > amIntercept:  boolean flag indicating if the anti-meridian is
                         intercepted by the observation path
         > inst:         string name of the instrument
         > sc:           string name of the spacecraft
         > target:       string name of the target body
         > resolution:   string defining the resolution setting, affecting the
                         footprint calculation. It could be either 'lowres' or
                         'highres'. See footprint function for further
                         information

        Returns:
          > A, tour, fpList, poly1, t: updated variables

        """
        ## Previous check...
        # if len(tour) == 0:
        #    empty = True
        #   return A, tour, fpList, poly1, t, empty

        # Compute the footprint of each point in the tour successively and
        # subtract the corresponding area from the target polygon
        a = copy.deepcopy(tour[0])  # observation
        tour.pop(0)  # delete this observation from the planned tour
        empty = False

        # Check a.m. intercept...
        if a[0] > 180:
            a[0] -= 360

        # Compute the observation's footprint
        print(f"Computing {inst} FOV projection on {target} at {M2P.mat2py_et2utc(t, 'C', 0)}...")
        fprinti = self.computeFootprint(t, inst, sc, target, resolution, a[0], a[1], 0)
        # Body-fixed to inertial frame
        if np.size(fprinti['bvertices']) != 0:  # assuming 'fprinti' is a dictionary with 'bvertices' key
            print("\n")

            # Check a.m. intercept
            if amIntercept:
                aux = copy.deepcopy(fprinti)
                ind = aux['bvertices'][:, 0] < 0
                aux['bvertices'][ind, 0] += 360
                if (np.isnan(aux['bvertices'][:, 0])).any():
                    nanindex = np.where(np.isnan(aux['bvertices'][:, 0]))[0]
                    polygon_list = []
                    for i in range(len(nanindex)):
                        if i == 0:
                            polygon_list.append(Polygon(
                                list(zip(aux['bvertices'][:nanindex[0], 0], aux['bvertices'][:nanindex[0], 1]))))
                        else:
                            polygon_list.append(Polygon(list(
                                zip(aux['bvertices'][nanindex[i - 1] + 1:nanindex[i], 0],
                                    aux['bvertices'][nanindex[i - 1] + 1:nanindex[i], 1]))))
                    if ~ np.isnan(aux['bvertices'][-1, 0]):
                        polygon_list.append(Polygon(
                            list(zip(aux['bvertices'][nanindex[-1] + 1:, 0], aux['bvertices'][nanindex[-1] + 1:, 1]))))
                    poly2 = MultiPolygon(polygon_list)
                else:
                    poly2 = Polygon(aux['bvertices'])
            else:
                if (np.isnan(fprinti['bvertices'][:, 0])).any():  # create footprint polygon
                    nanindex = np.where(np.isnan(fprinti['bvertices'][:, 0]))[0]
                    polygon_list = []
                    for i in range(len(nanindex)):
                        if i == 0:
                            polygon_list.append(Polygon(list(
                                zip(fprinti['bvertices'][:nanindex[0], 0], fprinti['bvertices'][:nanindex[0], 1]))))
                        else:
                            polygon_list.append(Polygon(list(
                                zip(fprinti['bvertices'][nanindex[i - 1] + 1:nanindex[i], 0],
                                    fprinti['bvertices'][nanindex[i - 1] + 1:nanindex[i], 1]))))
                    if ~ np.isnan(fprinti['bvertices'][-1, 0]):
                        polygon_list.append(Polygon(list(zip(fprinti['bvertices'][nanindex[-1] + 1:, 0],
                                                             fprinti['bvertices'][nanindex[-1] + 1:, 1]))))
                    poly2 = MultiPolygon(polygon_list)
                else:
                    poly2 = Polygon(fprinti['bvertices'])

            poly2 = poly2.buffer(0)

            A.append(a)  # add it in the list of planned observations
            poly1 = (poly1.difference(poly2)).buffer(0)  # update uncovered area

            # Save footprint struct
            fpList.append(fprinti)

            # New time iteration
            if len(tour) != 0:
                p1 = [fprinti['olon'], fprinti['olat']]
                p2 = [tour[0][0], tour[0][1]]
                t += tobs + self.slewDur(p1, p2, t, tobs, inst, target, sc, slewRate)
        else:
            empty = True
            print(" Surface not reachable\n")

        return A, tour, fpList, poly1, t, empty

    def slewDur(self,p1, p2, t, tobs, inst, target, sc, slew_rate):
        """
        This function determines the time it takes for a spacecraft or an
        instrument mounted on a spacecraft to rotate from one pointing direction
        to another. The calculation is based on the initial and final pointing
        vectors, the slew rate of the spacecraft or instrument, and the current
        observation time.

        Assumptions:
        - We assume that the slew rate is constant
        - We assume that the spacecraft or instrument can slew at this rate in
          any direction (and simultaneously in the three directions)
        - No rotation constraints are considered

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         09/2023

        Usage:        tdur = slewDur(p1, p2, t, tobs, inst, target, sc, slew_rate)

        Inputs:
          > p1:           coordinates [longitude, latitude] of the initial
                          pointing direction, in [deg]
          > p2:           coordinates [longitude, latitude] of the final pointing
                          direction, in [deg]
          > t:            current time in ephemeris seconds past J2000 epoch
          > inst:         string name of the instrument
          > target:       string name of the target body
          > sc:           string name of the spacecraft
          > slew_rate:    slew rate of the spacecraft or instrument, in [rad/s]

        Output:
          > tdur:         duration required to complete the slew, in [sec]
        """

        # Pre-allocate variables
        maxit = 10
        epsilon = 1e-2

        # Rotation matrix corresponding to the initial pointing direction using the
        # instrument pointing information
        _, _, R1, _ = self.instpointing(inst, target, sc, t, p1[0], p1[1])

        # Initial slew duration
        init_slew = 10.
        for i in range(maxit):
            # Get final pointing matrix
            _, _, R2, _ = self.instpointing(inst, target, sc, t + tobs + init_slew, p2[0], p2[1])

            # Relative rotation matrix between the two positions
            Rdelta = np.transpose(R1) @ R2

            # Angle of rotation required
            angle = np.arccos((np.trace(Rdelta) - 1) / 2)
            if angle >= np.pi / 2:
                angle = np.pi - angle

            # Duration of the slew
            tdur = angle / slew_rate
            new_slew = tdur

            # Check if the slew duration has converged
            if abs(new_slew - init_slew) < epsilon:
                return tdur

            # Update estimate
            init_slew = new_slew

        return tdur  # Return the final computed duration
    def planSidewinderTour(self,target, roi, sc, inst, inittime):
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
        origin = np.array([0., 0.])  # initialize grid origin for grid generation
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
        targetArea = self.topo2inst(roi, cx, cy, target, sc, inst, inittime)
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
                polygon_list.append(
                    Polygon(list(zip(targetArea[nanindex[-1] + 1:, 0], targetArea[nanindex[-1] + 1:, 1]))))
            poly_aux = MultiPolygon(polygon_list)
        else:
            poly_aux = Polygon((list(zip(targetArea[:, 0], targetArea[:, 1]))))

        # poly_aux = poly_aux.buffer(0)

        origin[0], origin[1] = poly_aux.centroid.x, poly_aux.centroid.y

        # Get minimum width direction of the footprint
        angle, _, _, _ = self.minimumWidthDirection(targetArea[:, 0], targetArea[:, 1])
        # observation angle, influencing the orientation of the observation footprints and,
        # therefore, the coverage path orientation

        # Retrieve the field of view (FOV) bounds of the instrument and calculate
        # the dimensions of a reference observation footprint
        _, _, _, bounds = M2P.mat2py_getfov(M2P.mat2py_bodn2c(inst)[0], 4)  # get fovbounds in the instrument's reference frame
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

        gt1 = np.array([0., 0.])
        gt2 = np.array([0., 0.])
        # Closest polygon side to the spacecraft's ground track position (this
        # will determine the coverage path)
        gt1[0], gt1[1] = self.groundtrack(sc, inittime, target)  # initial ground track position
        gt2[0], gt2[1] = self.groundtrack(sc, inittime + 500, target)  # future ground track position
        gt1 = self.topo2inst(gt1, cx, cy, target, sc, inst, inittime)  # projected initial position
        gt2 = self.topo2inst(gt2, cx, cy, target, sc, inst, inittime + 500)  # projected future position

        # Calculate the closest side of the target area to the spacecraft's ground track,
        # determining the observation sweep direction
        sweepDir1, sweepDir2 = self.closestSide(gt1, gt2, targetArea, angle)

        # Focal plane grid discretization based on the reference footprint (FOV
        # plane) and specified overlap
        inst_grid, grid_dirx, grid_diry = self.grid2D(fpref, self.olapx, self.olapy, origin, targetArea)

        # Boustrophedon decomposition to generate grid traversal
        inst_tour = self.boustrophedon(inst_grid, sweepDir1, sweepDir2)

        # Convert grid and tour from instrument frame to topographical coordinates
        topo_tour = self.inst2topo([inst_tour], cx, cy, target, sc, inst, inittime)[0]

        # Remove empty elements from the tour, which may result from unobservable regions
        # within the planned path
        topo_tour = [point for point in topo_tour if point is not None]

        return topo_tour, inst_grid, inst_tour, grid_dirx, grid_diry, sweepDir1, sweepDir2

    def topo2inst(self,inputdata_, lon, lat, target, sc, inst, et):
        """
        This function transforms a set of points from the topographic coordinate
        system(latitude and longitude on the target body) to the instrument frame
        coordinates.

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         09/2022

        Usage:        outputData = topo2inst(inputdata, lon, lat, target, sc, inst, et)

        Inputs:
          > inputdata:     A list of lists  or ndarray of points in topographic
                           coordinates to be transformed. Each point is a row with
                           [longitude, latitude] format
          > lon:           longitude of the observation point or area center, in
                           [deg]
          > lat:           latitude of the observation point or area center, in
                           [deg]
          > target:        string name of the target body
          > sc:            string name of the spacecraft
          > inst:          string name of the instrument
          > et:            ephemeris time, TDB seconds past J2000 epoch

        Outputs:
          > outputData:    A list of lists or ndarray of the input points transformed
                           to the instrument frame coordinates. The format of the
                           output matches the input (cell array or matrix)

        """
        inputdata = copy.deepcopy(inputdata_)
        # Handle input data in list format, ensuring all empty entries are replaced
        # with [NaN, NaN]
        ii, jj, aux = [], [], []

        if isinstance(inputdata, list):
            aux = [[point if point else [np.nan, np.nan] for point in row] for row in inputdata]
            topoPoints = np.vstack([[point for point in sublist] for sublist in aux])
            for i in range(len(inputdata)):
                for j in range(len(inputdata[i])):
                    ii.append(i)
                    jj.append(j)
        else:
            if np.shape(inputdata) == (2,):
                inputdata = inputdata.reshape(1, 2)
            topoPoints = copy.deepcopy(inputdata)

        # Pre-allocate variables
        _, targetframe, _ = M2P.mat2py_cnmfrm(target)  # target frame ID in SPICE

        # Build focal plane
        fovbounds, boresight, rotmat, _ = self.instpointing(inst, target, sc, et, lon, lat)

        vertex, _ = M2P.mat2py_spkpos(sc, et, targetframe, 'NONE', target)
        point = vertex + fovbounds[:, 0]

        # Create a plane based on the boresight and a point in the focal plane
        plane = M2P.mat2py_nvp2pl(boresight, point)

        # For each topographic point, find its intersection with the focal plane
        spoint = np.zeros((topoPoints.shape[0], 3))

        for i in range(topoPoints.shape[0]):
            if not np.isnan(topoPoints[i]).any():
                dir = -(psoa.trgobsvec(topoPoints[i], et, target, sc))[0]
                found, spoint[i, :] = M2P.mat2py_inrypl(vertex, dir, plane)
                if not found:
                    print('No intersection')
            else:
                spoint[i, :] = np.full(3, np.nan)

        # Transform coordinates from body-fixed to instrument frame
        tArea = np.zeros([spoint.shape[0], 3])
        for i in range(spoint.shape[0]):
            if not np.isnan(spoint[i, :]).any():
                vpoint = -(vertex - spoint[i, :].T)  # vector from spacecraft to intersection point
                tArea[i, :] = np.linalg.inv(rotmat).dot(vpoint)  # apply inverse rotation to transform
                # to instrument frame
            else:
                tArea[i, :] = np.full(3, np.nan)

        instcoord = tArea[:, :2]  # extract 2D instrument frame coordinates
        # Prepare output data matching the format of the input,i.e., cell array or
        # matrix
        outputData = [[None for _ in row] for row in inputdata]
        if isinstance(inputdata, list):
            for k in range(len(ii)):
                if not np.isnan(instcoord[k]).any():
                    outputData[ii[k]][jj[k]] = instcoord[k, :]
        else:
            outputData = copy.deepcopy(instcoord)
        return outputData

    def inst2topo(self,grid, lon, lat, target, sc, inst, et):
        """
        This function transforms a set of points from the instrument frame to the
        topographic coordinate system (latitude and longitude on the target body)

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         09/2022

        Usage:        grid_topo = inst2topo(grid, lon, lat, target, sc, inst, et)

        Inputs:
            grid:          A list of grid points in instrument frame coordinates.
                           Each element contains a 2D point [x,y] representing a
                           location in the instrument frame, or is empty if no
                           observation point is defined
            lon:           Longitude of the observation point or area center, in [deg].
            lat:           Latitude of the observation point or area center, in [deg].
            target:        String name of the target body.
            sc:            String name of the spacecraft.
            inst:          String name of the instrument.
            et:            Ephemeris time, TDB seconds past J2000 epoch.

        Returns:
            grid_topo:     A list of the input grid points transformed to topographic
                           coordinates on the target body. Each element contains a 2D point
                           [lon,lat], in [deg]
        """

        # Pre-allocate
        _, _, rotmat, _ = self.instpointing(inst, target, sc, et, lon, lat)  # Assuming instpointing function is defined
        grid_topo = [[[] for _ in range(np.shape(grid)[1])] for _ in
                     range(np.shape(grid)[0])]  # Pre-allocate grid_topo array
        method = 'ELLIPSOID'
        _, targetframe, _ = M2P.mat2py_cnmfrm(target)  # Get target frame ID in SPICE

        # Convert grid into topographical coordinates
        for i in range(len(grid)):
            for j in range(len(grid[i])):
                instp = grid[i][j]  # retrieve current point in instrument frame
                if instp is not None and not np.all(np.isnan(instp)):  # check for empty or NaN
                    p = np.zeros(3)  # initialize 3D points for conversion
                    p[0:2] = instp  # assign [x, y] coordinates
                    p[2] = 1.  # set z to 1 for surface projection calculation
                    p_body = np.dot(rotmat, p)  # apply rotation matrix

                    # Compute surface intersection of the point on the target body
                    xpoint, _, _, found = M2P.mat2py_sincpt(method, target, et, targetframe,
                                                        'NONE', sc, targetframe, p_body)

                    if found:
                        # Convert rectangular coordinates to latitudinal
                        _, lon, lat = M2P.mat2py_reclat(xpoint)
                        grid_topo[i][j] = [lon * M2P.mat2py_dpr(), lat * M2P.mat2py_dpr()]
                    else:
                        print("Point not visible from the instrument")
                        grid_topo[i][j] = None

        return grid_topo
    def computeFootprint(self,t, inst, sc, target, res, *args):
        """
        Given the spacecraft trajectory, this function computes the FOV
        projection onto the body surface, i.e., the footprint. The spacecraft
        orientation is determined by the lon, lat, and theta angles.
        Assumption: the FOV is rectangular.
        Note: this function is contained in the framework of the automated
        scheduler that optimizes the observation plan of a mission according to
        a set of scientific objectives.
        [Warning]: limb projections on a topographical map might be very
        irregular. FUTURE WORK

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         10/2022
        Version:      2
        Last update:  12/2023

        Usage:        fp = footprint(t, inst, sc, target, res)
                      fp = footprint(t, inst, sc, target, res, lon, lat)
        Inputs:
          > t:        time epoch in TDB seconds past J2000 epoch
          > inst:     string SPICE name of the instrument
          > sc:       string SPICE name of the spacecraft
          > target:   string SPICE name of the target body
          > lon:      longitude coordinate of the target body at which the
                      instrument boresight is pointing, in [deg]
          > lat:      latitude coordinate of the target body at which the
                      instrument boresight is pointing, in [deg]

        Outputs:
          > fp:       struct containing main parameters of the footprint. In this
                      case, only the field 'vertices' is necessary.
              # inst:          string SPICE name of the instrument
              # sc:            string SPICE name of the spacecraft
              # target:        string SPICE name of the target body
              # t:             time epoch in TDB seconds past J2000 epoch
              # bvertices:     matrix (N, 2) containing the N boundary vertices
                               of the footprint polygon, in latitudinal
                               coordinates, in [deg].
                      + fp may contain more than one polygon, split because the
                        actual footprint intersects the anti-meridian of the body
                        surface. In this case, the polygons are separated in
                        'vertices' by [NaN, NaN] Example:
                                area1 = [350 30; 360 30; 360  0; 350 0];
                                area2 = [ 0 30; 10 30; 10  0; 0  0];
                                area = [area1; [NaN NaN]; area2];
                                figure plot(polyshape(area(:,1), area(:,2));
              # olon:          longitude coordinate of the target body at which
                               the instrument boresight is pointing, in [deg]
              # olat:          latitude coordinate of the target body at which
                               the instrument boresight is pointing, in [deg]
               # fovbsight:    FOV boresight in the body-fixed reference frame
                               centered at the spacecraft position at time t
              # fovbounds:     FOV bounds in the body-fixed reference frame
                               centered at the spacecraft position at time t

       """

        # Pre-allocate variables
        _, targetframe, _ = M2P.mat2py_cnmfrm(target)  # target frame ID in SPICE
        abcorr = 'LT+S'  # one-way light time aberration correction parameter
        method = 'ELLIPSOID'  # assumption: ray intercept function is going to model the target
        # body as a tri-axial ellipsoid
        surfPoints = np.array([])  # matrix that saves the rectangular coordinates of the intercept points between
        # the FOV perimeter and the body surface
        count = 0  # counter of found intercept points
        geom = False
        ckpointing, lon, lat = [], [], []

        # Define instrument pointing (3-axis steerable or constrained)
        if len(args) == 0:  # instrument pointing is provided by (and retrieve from) a ck
            ckpointing = True
        elif len(args) >= 1:  # instrument is considered 3-axis steerable
            lon, lat = args[0], args[1]
            ckpointing = False
            if len(args) == 3:
                geom = args[2]

        # Definition of footprint resolution
        if res == 'lowres':  # footprint vertices resolution
            N = 10  # number of intercept searches per side
        elif res == 'highres':
            N = 500
        else:
            raise ValueError("Invalid resolution method")

        # Initialize footprint dictionary keys
        fp = {
            'inst': inst,
            'sc': sc,
            'target': target,
            't': t,
            'bvertices': np.array([]),  # footprint boundary vertices, in latitudinal coordinates, in [deg]
            'olon': np.nan,  # longitude value of FOV boresight projection onto the body surface, in [deg]
            'olat': np.nan,  # latitude value of FOV boresight projection onto the body surface, in [deg]
            'fovbsight': np.array([]),  # FOV boresight in target frame centered at the spacecraft position
            'fovbounds': np.array([]),  # FOV bounds in target frame centered at the spacecraft position
            'limb': 'none',  # boolean that defines if the FOV projects onto the planetary body's limb
            'recVertices': np.array([]),
            # matrix that contains the rectangular coordinates of the footprint boundary vertices
            'angle': np.nan,
            'width': np.nan,
            'height': np.nan
        }

        # Calculate instrument orientation
        if ckpointing:  # constrained pointing
            bounds, boresight, pointingRotation, found, lon, lat = self.instpointing(inst, target, sc, t)
        else:
            bounds, boresight, pointingRotation, found = self.instpointing(inst, target, sc, t, lon, lat)

        if not found:
            return fp  # the point is not visible from the instrument's FOV, therefore
            # the function is exited and the footprint is returned empty

        fp['fovbounds'] = bounds  # save fov bounds in the body-fixed reference frame
        fp['fovbsight'] = boresight  # save fov boresight vector in the body-fixed reference frame

        # Update boresight pointing (in latitudinal coordinates)
        fp['olon'] = lon  # longitude value of FOV boresight projection onto the body surface, in [deg]
        fp['olat'] = lat  # latitude value of FOV boresight projection onto the body surface, in [deg]

        # Project instrument FOV onto body surface

        # Retrieve FOV parameters
        _, _, _, bounds = M2P.mat2py_getfov((M2P.mat2py_bodn2c(inst))[0], 4)  # instrument FOV's boundary
        # vectors in the instrument frame
        minx = np.min(bounds[0, :])  # minimum x focal plane
        maxx = np.max(bounds[0, :])  # maximum x focal plane
        miny = np.min(bounds[1, :])  # minimum y focal plane
        maxy = np.max(bounds[1, :])  # maximum y focal plane
        z = bounds[2, 0]  # z-coordinate of the boundary vectors

        boundPoints = np.zeros((3, max(bounds.shape)))  # intercept points of the FOV
        # boundary, in Cartesian coordinates, and in the body-fixed reference frame
        # In its simplest form, the footprint should have the same number of
        # vertices as boundaries has the instrument's FOV

        fp['limb'] = 'none'  # boolean to define if the instrument FOV projection is not
        # enclosed in the body surface, i.e. at least one of the boundary vectors
        # do not intercept the body surface

        intsec = False  # boolean that indicates if at least one of the boundary
        # vectors intercepts the body surface

        for i in range(max(fp['fovbounds'].shape)):
            boundPoints[:, i], _, _, found = M2P.mat2py_sincpt(method, target, t, targetframe, abcorr, sc, targetframe,
                                                           fp['fovbounds'][:, i])
            # If the FOV boundary does not intercept the
            # object's surface, then we're seeing (at least, partially)
            # the limb
            if not found:
                fp['limb'] = 'partial'
            else:
                intsec = True

        # Assume that, if there are no intercepts of the FOV boundaries, the FOV
        # projection is likely to contain the total limb of the body
        # In the next step, we will find out if this assumption is correct
        if not intsec:
            fp['limb'] = 'total'

        # When the footprint is likely to contain the limb, perform a more
        # accurate search in order to conclude if the FOV intercept the
        # body at some point
        if fp['limb'] == 'total':
            fp, maxx, minx, maxy, miny, pointingRotation, method, target, t, targetframe, abcorr, sc = self.refineFOVsearch(
                fp, maxx, minx, maxy, miny, pointingRotation, method, target, t, targetframe, abcorr, sc)
            # Compute footprint
        elif fp['limb'] == 'none':
            # FOV projects entirely on the body surface
            boundPoints, N, surfPoints, count = self.inFOVprojection(boundPoints, N, surfPoints, count)
        elif fp['limb'] == 'partial':
            # FOV contains partially the limb
            surfPoints, pointingRotation, minx, maxx, miny, maxy, z, N, t, method, target, targetframe, abcorr, sc, res, count = self.plimbFOVprojection(
                surfPoints, pointingRotation, minx, maxx, miny, maxy, z, N, t, method, target, targetframe, abcorr, sc,
                res, count)
        else:
            # FOV contains the whole body (total limb)
            surfPoints, target, t, targetframe, sc = self.tlimbFOVprojection(surfPoints, target, t, targetframe, sc)

        if surfPoints.size == 0:  # the FOV does not intercept with the object at any point
            # of its focal plane
            return fp

            # Save values
        fp['recVertices'] = surfPoints
        # Conversion from rectangular to latitudinal coordinates of the polygon vertices
        # and geometry computation
        surfPoints, t, target, sc, fp, inst = self.footprint2map(surfPoints, t, target, sc, fp, inst)
        return fp
    @staticmethod
    def instpointing(inst, target, sc, t, *args):
        """
        This function sets the instrument's orientation, provided a target and
        the latitudinal coordinates of the point the instrument should be aiming
        at. It also checks if this point is actually visible from the FOV (could
        be on the dark side of the object as seen from the instrument).

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         01/2023

        Usage:        [fovbounds, boresight, rotmat, visible] = instorient(inst,
                          target, sc, t, lon, lat)
                      [fovbounds, boresight, rotmat, visible, lon, lat] = instorient(inst,
                          target, sc, t)

        Inputs:
          > inst:       string name of the instrument
          > target:     string name of the target body
          > sc:         string name of the spacecraft
          > t:          observation time, i.e. the minimum time that the
                        instrument needs to perform an observation, in seconds
          > lon:        longitude coordinate of the target body at which the
                        instrument boresight is pointing, in [deg]
          > lat:        latitude coordinate of the target body at which the
                        instrument boresight is pointing, in [deg]

        Returns:
          > fovbounds:  FOV bounds in the body-fixed reference frame centered at
                        the spacecraft position at time t
          > boresight:  FOV boresight in the body-fixed reference frame
                        centered at the spacecraft position at time t
          > rotmat:     rotation matrix (numpy.ndarray) from instrument frame to target frame
          > visible:    boolean that determines if the point is visible from the
                        instrument's FOV
          > lon:        longitude coordinate of the target body at which the
                        instrument boresight is pointing, in [deg]
          > lat:        latitude coordinate of the target body at which the
                        instrument boresight is pointing, in [deg]
        """
        # Pre-allocate variables
        axis3 = False  # boolean variable that indicates if the spacecraft is 3-axis steerable
        lon = 0
        lat = 0
        if len(args) > 0:
            lon = args[0]
            lat = args[1]
            axis3 = True

        method = 'ELLIPSOID'  # assumption: ray intercept function is going to
        # model the target body as a tri-axial ellipsoid
        _, targetframe, _ = M2P.mat2py_cnmfrm(target)  # target frame ID in SPICE
        abcorr = 'LT'  # one-way light time aberration correction parameter.

        # Retrieve FOV parameters
        shape, instframe, boresight, bounds = M2P.mat2py_getfov((M2P.mat2py_bodn2c(inst))[0], 4)  # instrument FOV's boundary
        # vectors in the instrument frame
        if shape in ["CIRCLE", "ELLIPSE"]:
            raise ValueError("Circular and ellipsoidal FOV shapes have not been implemented yet")

        fovbounds = np.zeros((3, max(np.shape(bounds))))
        rotmat = np.zeros((3, 3))
        visible = False

        # Pointing matrix
        if axis3:  # 3-axis steerable
            lon = lon * M2P.mat2py_rpd()
            lat = lat * M2P.mat2py_rpd()

            # Boresight of the instrument must point at the target point
            recpoint = M2P.mat2py_srfrec(M2P.mat2py_bodn2c(target)[0], lon, lat)  # rectangular
            # coordinates of the target point in the body-fixed reference frame

            instpos, _ = M2P.mat2py_spkpos(sc, t, targetframe, abcorr, target)  # rectangular coordinates
            # of the instrument in the body-fixed reference frame
            v1 = recpoint - instpos  # distance vector to the target point from the instrument in the body-fixed reference frame
            boresight = v1 / np.linalg.norm(v1)  # boresight of the instrument in the
            # body-fixed reference frame

            # The z-axis of the instrument is the boresight
            rotmat[:, 2] = boresight

            # Define a consistent reference vector, e.g., [0, 0, 1] (celestial north)
            reference_vector = np.array([0, 0, 1])

            # Check if boresight is aligned or anti-aligned with reference vector
            if abs(np.dot(boresight, reference_vector)) > 0.999:
                # Adjust reference vector if aligned/anti-aligned to avoid singularity
                reference_vector = np.array([0, 1, 0])

            # Define y-axis using cross product to ensure perpendicularity
            yinst = np.cross(boresight, reference_vector)
            yinst = yinst / np.linalg.norm(yinst)  # normalize yinst

            # Define x-axis using cross product between boresight and yinst
            xinst = np.cross(yinst, boresight)
            xinst = xinst / np.linalg.norm(xinst)  # normalize xinst

            # Assign to rotation matrix
            rotmat[:, 0] = xinst
            rotmat[:, 1] = yinst
            rotmat[:, 2] = boresight

        else:  # Pointing constrained by ckernel
            rotmat = M2P.mat2py_pxform(instframe, targetframe, t)
            xpoint, _, _, found = M2P.mat2py_sincpt(method, target, t, targetframe, abcorr, sc, instframe, boresight)

            if found:
                _, lon, lat = M2P.mat2py_reclat(xpoint)
                lon = lon * M2P.mat2py_dpr()
                lat = lat * M2P.mat2py_dpr()
            else:
                # not visible
                print(f"On {M2P.mat2py_et2utc(t, 'C', 0)}, {inst} is not pointing at {target}")
                if len(args) > 0:
                    return fovbounds, boresight, rotmat, visible
                else:
                    return fovbounds, boresight, rotmat, visible, lon, lat

            # Boresight of the instrument must point at the target point
            recpoint = xpoint  # rectangular coordinates of the target point in the body-fixed reference frame
            instpos, _ = M2P.mat2py_spkpos(sc, t, targetframe, abcorr,
                                       target)  # rectangular coordinates of the instrument in the body-fixed reference frame
            v1 = recpoint - instpos  # distance vector to the target point from the
            # instrument in the body-fixed reference frame

        # Transform coordinates
        fovbounds = np.zeros([3, max(np.shape(bounds))])
        for i in range(max(np.shape(bounds))):
            fovbounds[:, i] = np.dot(rotmat, bounds[:, i])  # instrument FOV's boundary vectors in the target frame

        # Check if the point is visible as seen from the instrument
        if np.dot(v1, recpoint) > 0:  # check if the point is visible as seen from the instrument
            if len(args) > 0:
                return fovbounds, boresight, rotmat, visible
            else:
                return fovbounds, boresight, rotmat, visible, lon, lat
        else:
            visible = True

        # Output values
        if len(args) > 0:
            return fovbounds, boresight, rotmat, visible
        else:
            return fovbounds, boresight, rotmat, visible, lon, lat

    def minimumWidthDirection(self,x_, y_):
        """
        This function computes the orientation (angle) at which the width of the
        polygon is minimized. It also calculates the minimum width, the height
        (assuming the minimum width as the base), and the direction vector (axes)
        of the minimum width

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         09/2022

        Usage:        [thetamin, minwidth, height, axes] = minimumWidthDirection(x, y)

        Inputs:
          > x:            array of x-coordinates of the polygon vertices
          > y:            array of y-coordinates of the polygon vertices

        Outputs:
          > thetamin:     angle at which polygon's width is minimized, in [deg]
          > minwidth:     minimum width of the polygon
          > height:       height of the polygon, at the orientation specified by
                          thetamin and assuming minwidth as the base
          > axes:         unit vector representing the direction of the polygon
                          axes (width, height)
       """
        x = copy.deepcopy(x_)
        y = copy.deepcopy(y_)

        # Check if the polygon is divided in two (a.m. intersection)...
        ind = np.where(np.isnan(x))[0]
        if ind.size > 0:
            x[:ind[0]] += 360
            x = np.delete(x, ind)
            y = np.delete(y, ind)

        # Find centroid
        poly_aux = Polygon(list(zip(x, y)))
        cx = poly_aux.centroid.x
        cy = poly_aux.centroid.y

        vertices = np.zeros([np.size(x), 2])
        # Sort the vertices in clockwise direction
        vertices[:, 0], vertices[:, 1] = self.sortcw(x, y)

        # Minimum width direction
        npoints = 361
        angle = np.linspace(0, 360, npoints)
        cosa = np.cos(np.deg2rad(angle))
        sina = np.sin(np.deg2rad(angle))
        rotMats = np.zeros((2, 2, npoints))
        for i in range(npoints):
            rotMats[:, :, i] = np.array([[cosa[i], sina[i]], [-sina[i], cosa[i]]])

        minwidth = np.inf
        mini = np.nan
        for i in range(npoints):
            # Compute rotation matrix to orient the convex hull with the main x-y
            # axis (this way it is easier to find the bounding box with the min and
            # max values of x/y)
            rotVertices = (rotMats[:, :, i] @ (vertices - [cx, cy]).T).T

            # Obtain width length of the rotated polygon
            maxx = np.max(rotVertices[:, 0])
            minx = np.min(rotVertices[:, 0])
            l = maxx - minx
            if l < minwidth:
                minwidth = l
                mini = i

        # Return minimum width direction
        thetamin = angle[mini]
        polygon = Polygon(np.column_stack((x, y))).buffer(0)
        area = polygon.area
        height = area / minwidth

        # Constrain angle between 0º and 180º
        if thetamin >= 180:
            thetamin -= 180

        axes = np.array([np.cos(np.deg2rad(thetamin)), np.sin(np.deg2rad(thetamin))])

        return thetamin, minwidth, height, axes
    @staticmethod
    def refineFOVsearch(fp, maxx, minx, maxy, miny, pointingRotation, method, target, t, targetframe, abcorr,
                        sc):
        # Perform a perimetral search of the FOV to find out if the FOV
        # contains totally or partially the body limb
        Nl = 20
        found = False
        for ii in range(2):
            # Vertical sweep of the focal perimeter
            x = (maxx - minx) * ii + minx
            for jj in range(Nl + 1):
                y = (maxy - miny) * jj / Nl + miny
                vec = np.array([x, y, z]).reshape([3, 1])
                vec = np.dot(pointingRotation, vec)  # transform vector coordinates to target frame
                _, _, _, found = M2P.mat2py_sincpt(method, target, t, targetframe, abcorr, sc, targetframe,
                                               vec)
                if found:
                    fp['limb'] = 'partial'
                    break
            if found:
                break

        # If intercept still has not been found...
        if not found:
            for ii in range(Nl + 1):
                # Horizontal sweep of the focal perimeter
                x = (maxx - minx) * ii / Nl + minx
                for jj in range(2):
                    y = (maxy - miny) * jj / Nl + miny
                    vec = np.array([x, y, z]).reshape([3, 1])
                    vec = np.dot(pointingRotation, vec)  # transform vector coordinates to target frame
                    _, _, _, found = M2P.mat2py_sincpt(method, target, t, targetframe, abcorr, sc,
                                                   targetframe, vec)
                    if found:
                        fp['limb'] = 'partial'
                        break
                if found:
                    break
        return fp, maxx, minx, maxy, miny, pointingRotation, method, target, t, targetframe, abcorr, sc

    @staticmethod
    def plimbFOVprojection(surfPoints,pointingRotation, minx, maxx, miny, maxy, z, N, t, method, target, targetframe, abcorr, sc,
                           res,count):
        # Warning message
        if res == 'lowres':
            warnings.warn(
                "Warning: It is likely that the footprint contains the limb, low resolution method may lead to significant inaccuracies")

        # Initialize variables
        maxfx, maxfy = minx, miny
        minfx, minfy = maxx, maxy
        #count = 0

        old_found = False
        old_surfPoint = None

        # For those cases where the FOV does not completely contain the target,
        # a more refined search is going to be performed in order to define the
        # limits of the footprint

        for i in range(N + 1):
            # Vertical sweep of the focal plane
            x = (maxx - minx) * i / N + minx
            for j in range(N + 1):
                y = (maxy - miny) * j / N + miny
                vec = np.array([x, y, z]).reshape([3,1])
                vec = np.dot(pointingRotation,vec) # transform vector coordinates to
                # target frame

                corloc = 'SURFACE POINT' # since alt is close to 0, there
                                         # shoul not be a significant difference between the target and
                                         # surface point correction locus (see spice.tangpt)
                found = False # found intercept

                _,alt,_, aux,_,_ = M2P.mat2py_tangpt(method, target, t, targetframe, abcorr, corloc, sc, targetframe,
                                          vec)
                if alt < 15:
                    # When the footprint contains the limb, its intercept is
                    # irregular, meaning that the boundary is not a smooth
                    # curve (the limb) but a set of scattered points with
                    # certain deviation around the limb. Besides this, since
                    # the research consists of a set of discretized points, we
                    # may not always find the limb intercept point. This
                    # depends on the resolution of the discretized refined
                    # mesh. To avoid incurring in excessive computational
                    # demands, instead of calculating the intercept point by
                    # refining the mesh, we calculate the tangent point. This
                    # is the closest surface point of the surface to the
                    # "intercepting" ray. When the ray actually intercepts the
                    # surface, the parameter 'alt', which is the distance
                    # between the tangent points and the surface, is equal to
                    # 0. We may find the limb "intercept" by finding those
                    # points where 'alt' is close or equal to 0.
                    # Future work: to be more precise... try to minimize alt
                    # along the vertical sweep line
                    found = True

                if found and (y == miny or y == maxy or x == minx or x == maxx):
                    # if the vector intercepts the surface and is at the focal
                    # plane boundary...
                    count += 1
                    if np.size(surfPoints) == 0:
                        surfPoints = aux
                    else:
                     surfPoints = np.vstack((surfPoints,aux))

                    if y == miny:
                        minfy = miny
                    elif y == maxy:
                        maxfy = maxy

                    if x == minx:
                        minfx = minx
                    elif x == maxx:
                        maxfx = maxx

                elif j > 0 and found != old_found:
                    # if the vector intercept status changes from the previous
                    # one, we're sweeping across the object's limb
                    count += 1
                    if old_found:
                        if np.size(surfPoints) == 0:
                            surfPoints = old_surfPoint # save the previous intercept
                        else:
                            surfPoints = np.vstack((surfPoints, old_surfPoint)) # save the previous intercept
                    else:
                        if np.size(surfPoints) == 0:
                            surfPoints = aux # save the current intercept
                        else:
                            surfPoints = np.vstack((surfPoints, aux)) # save the current intercept

                    if x < minfx:
                        minfx = x
                    if y < minfy:
                        minfy = y
                    if x > maxfx:
                        maxfx = x
                    if y > maxfy:
                        maxfy = y

                old_found = found # save element intercept status
                old_surfPoint = aux # save element intercept point
        return surfPoints, pointingRotation, minx, maxx, miny, maxy, z, N, t, method, target, targetframe, abcorr, sc, res, count

    @staticmethod
    def inFOVprojection(self,boundPoints,N,surfPoints,count):
        """
        The FOV projection is enclosed in the target surface.
        """
        # Close polygon
        boundPoints = np.hstack((boundPoints, boundPoints[:,[0]]))
        #count = 0
        surfPoints = np.zeros([N*(max(boundPoints.shape) - 1),3])
        # high resolution
        for i in range(max(boundPoints.shape) - 1):
            # linear (approximation) interpolation between vertices to define
            # the boundary of the footprint
            v = boundPoints[:, i + 1] - boundPoints[:, i]
            lambda_vals = np.linspace(0, 1, N)  # line parametrization
            for l in range(N):
                surfPoints[count,0] = boundPoints[0, i] + v[0] * lambda_vals[l]
                surfPoints[count,1] = boundPoints[1, i] + v[1] * lambda_vals[l]
                surfPoints[count,2] = boundPoints[2, i] + v[2] * lambda_vals[l]
                count += 1
        return boundPoints,N,surfPoints,count
    @staticmethod
    def tlimbFOVprojection(surfPoints, target, t, targetframe, sc):
        # Compute limb with SPICE function (easier)
        # Parameters for spice.limbpt function
        lbmethod = 'TANGENT/ELLIPSOID'
        abcorr = 'XLT+S'
        corloc = 'CENTER'
        refvec = np.array([0, 0, 1]).reshape(3, 1)  # first of the sequence of cutting half-planes
        ncuts = int(2e3)  # number of cutting half-planes
        delrol = M2P.mat2py_twopi() / ncuts  # angular step by which to roll the
        # cutting half-planes about the observer-target vector
        schstp = 1.0e-6  # search angular step size
        soltol = 1.0e-10  # solution convergence tolerance

        # Limb calculation with spice.limbpt function
        _, limb, _, _ = M2P.mat2py_limbpt(lbmethod, target, t, targetframe, abcorr, corloc, sc, refvec, delrol, ncuts,
                                      schstp, soltol, ncuts)
        # limb points expressed in targetframe ref frame
        surfPoints = limb.T
        return surfPoints, target, t, targetframe, sc


    def footprint2map(self,surfPoints,t,target,sc,fp,inst):

        # Pre-allocate variables
        vertices = np.zeros([max(np.shape(surfPoints)),2]) # matrix that saves the
        # latitudinal coordinates of the intercept points between the FOV
        # perimeter and the body surface

        # Sort points
        surfPoints[:,0],surfPoints[:,1],surfPoints[:,2] = self.sortcw(surfPoints[:,0],surfPoints[:,1],surfPoints[:,2])
        # sort polygon boundary vertices in clockwise order (for representation)

        for i in range(max(np.shape(surfPoints))):
            _, auxlon, auxlat = M2P.mat2py_reclat(surfPoints[i, :].T)  # rectangular to
            # latitudinal coordinates
            vertices[i, 0] = auxlon * M2P.mat2py_dpr()  # longitude in [deg]
            vertices[i, 1] = auxlat * M2P.mat2py_dpr()  # latitude in [deg]
        # Future work: surfPoints does not need to be saved, we could convert
        # from rectangular to latitudinal inside the first loop, instead of
        # doing separately. The reason why it is not is because we need to sort
        # the vertices in clockwise order, and the sortcw algorithm for 2D does
        # not work with non-convex polygons...

        if fp['limb'] == 'total':
            lblon = copy.deepcopy(vertices[:,0])
            lblat = copy.deepcopy(vertices[:,1])

            # We need to discern between two different limbs:
            # 1.- Sub-spacecraft point is located at the equator (observer-to-pole line is
            # perpendicular to the normal vector at the poles). In this case, the limb's
            # longitude cannot be > 180º.
            # 2.- Sub-spacecraft point is not located at the equator. In this case, the limb's
            # longitude may be > 180º (and includes the north/south poles).

            northpole = False
            southpole = False
            # Check north-pole
            srfpoint = np.array([0, 90])
            angle = psoa.emissionang(srfpoint, t, target, sc)
            if angle < 90:
                northpole = True
            # Check south-pole
            srfpoint = np.array([0, -90])
            angle = psoa.emissionang(srfpoint, t, target, sc)
            if angle < 90:
                southpole = True

            # Case 1.
            if not northpole and not southpole:
                # Check a.m. split
                ind2 = np.where(np.diff(np.sort(lblon)) >= 180)[0]  # find the discontinuity
                if ind2.size > 0:
                    lblon, lblat = self.amsplit(lblon, lblat)
                # Check if we are keeping the correct polygon (full disk polygons may be
                # misleading, we can only guarantee through emission angle check)
                ext = False

                while not ext:
                    randPoint = np.array([np.random.randint(-180, 181), np.random.randint(-90, 91)])
                    point = Point(randPoint)
                    if (np.isnan(lblon)).any():
                        nanindex = np.where(np.isnan(lblon))[0]
                        polygon_list = []
                        for i in range(len(nanindex)):
                            if i == 0:
                                polygon_list.append(Polygon(list(zip(lblon[:nanindex[0]], lblat[:nanindex[0]]))))
                            else:
                                polygon_list.append(Polygon(
                                    list(zip(lblon[nanindex[i - 1] + 1:nanindex[i]],
                                             lblat[nanindex[i - 1] + 1:nanindex[i]]))))
                        if ~ np.isnan(lblon[-1]):
                            polygon_list.append(Polygon(list(zip(lblon[nanindex[-1] + 1:], lblat[nanindex[-1] + 1:]))))
                        polyaux = MultiPolygon(polygon_list)
                    else:
                        polyaux = Polygon((list(zip(lblon, lblat))))
                    polyaux = polyaux.buffer(0)

                    if polyaux.intersects(point):
                        angle = psoa.emissionang(randPoint, t, target, sc)
                        if angle < 85:
                            ext = True
                    else:
                        angle =  psoa.emissionang(randPoint, t, target, sc)
                        if angle < 85:
                            ext = True
                            # This calculation is approximated, we should find a better way
                            # to find the complementary
                            # [Future work]
                            lonmap = [-180, -180, 180, 180]
                            latmap = [-90, 90, 90, -90]
                            polymap = Polygon(list(zip(lonmap, latmap)))
                            if (np.isnan(lblon)).any():
                                nanindex = np.where(np.isnan(lblon))[0]
                                polygon_list = []
                                for i in range(len(nanindex)):
                                    if i == 0:
                                        polygon_list.append(Polygon(list(zip(lblon[:nanindex[0]], lblat[:nanindex[0]]))))
                                    else:
                                        polygon_list.append(Polygon(
                                            list(zip(lblon[nanindex[i - 1] + 1:nanindex[i]],
                                                     lblat[nanindex[i - 1] + 1:nanindex[i]]))))
                                if ~ np.isnan(lblon[-1]):
                                    polygon_list.append(
                                        Polygon(list(zip(lblon[nanindex[-1] + 1:], lblat[nanindex[-1] + 1:]))))
                                poly1 = MultiPolygon(polygon_list)
                            else:
                                poly1 = Polygon((list(zip(lblon, lblat))))
                            poly1 = poly1.buffer(0)
                            poly1 = polymap.difference(poly1)
                            poly1 = poly1.buffer(0)

                            if isinstance(poly1, Polygon):
                                lblon, lblat = np.array(poly1.exterior.coords.xy)
                            elif isinstance(poly1, MultiPolygon):
                                for i in range(len(poly1.geoms)):
                                    lblonaux, lblataux = np.array(poly1.geoms[i].exterior.coords.xy)
                                    if i == 0:
                                        lblon = np.append(lblonaux, np.nan)
                                        lblat = np.append(lblataux, np.nan)
                                    else:
                                        lblon = np.append(lblon, np.append(lblonaux, np.nan))
                                        lblat = np.append(lblat, np.append(lblataux, np.nan))
                                lblon = lblon[:-1]
                                lblat = lblat[:-1]
            else:
                # Case 2.
                lblon, indsort = np.sort(lblon), np.argsort(lblon)
                lblat = lblat[indsort]

                if northpole or southpole:
                    # Include northpole to close polygon
                    auxlon, auxlat = copy.deepcopy(lblon), copy.deepcopy(lblat)
                    lblon = np.zeros(len(auxlon) + 2)
                    lblat = np.zeros(len(auxlat) + 2)
                    if northpole:
                        lblon[0], lblat[0] = -180, 90
                        lblon[-1], lblat[-1] = 180, 90
                    else:
                        lblon[0], lblat[0] = -180, -90
                        lblon[-1], lblat[-1] = 180, -90
                    lblon[1:-1] = auxlon
                    lblat[1:-1] = auxlat
            fp['bvertices'] = np.hstack((lblon.reshape(len(lblon),1),lblat.reshape(len(lblat),1)))
        elif fp['limb'] == 'partial':
            lblon = copy.deepcopy(vertices[:,0])
            lblat = copy.deepcopy(vertices[:,1])
            # Check north-pole visibility:
            northpole = M2P.mat2py_fovray(inst, target, sc, t, 0, 90, fp['olon'], fp['olat'])
            # Check south-pole visibility:
            southpole = M2P.mat2py_fovray(inst, target, sc, t, 0, -90, fp['olon'], fp['olat'])
            # Case 1.
            if not northpole and not southpole:
                lblon, lblat = self.amsplit(lblon, lblat)
            else:
                # Case 2.
                lblon, indsort = np.sort(lblon), np.argsort(lblon)
                lblat = lblat[indsort]

                if northpole or southpole:
                    # Include northpole to close polygon
                    auxlon, auxlat = copy.deepcopy(lblon), copy.deepcopy(lblat)
                    lblon = np.zeros(len(auxlon) + 2)
                    lblat = np.zeros(len(auxlat) + 2)
                    if northpole:
                        lblon[0], lblat[0] = -180, 90
                        lblon[-1], lblat[-1] = 180, 90
                    else:
                        lblon[0], lblat[0] = -180, -90
                        lblon[-1], lblat[-1] = 180, -90
                    lblon[1:-1] = auxlon
                    lblat[1:-1] = auxlat
            fp['bvertices'] = np.hstack((lblon.reshape(len(lblon), 1), lblat.reshape(len(lblat), 1)))
        else:
            # Check if the footprint intersects the anti-meridian
            # To ease the footprint representation on the topography map, we must
            # consider the case where the footprint intercepts with the anti-meridian.
            # If it does, we split the footprint in two polygons, cleaved by the line
            # that the original footprint is crossing (a.m.)

            col1, col2 = self.amsplit(vertices[:, 0], vertices[:, 1])  # save footprint vertices
            fp['bvertices'] = np.hstack((col1.reshape(len(col1), 1), col2.reshape(len(col2), 1)))

            if geom:
                # Get minimum width direction and size
                angle, width, height, _ = self.minimumWidthDirection(fp['bvertices'][:, 0], fp['bvertices'][:, 1])
                fp['angle'] = angle
                fp['width'] = width
                fp['height'] = height

        return surfPoints,t,target,sc,fp,inst

    def computeResMosaic(self, fpList, ifov):
        r = []
        for fp in fpList:
            srfpoint = [fp['olon'], fp['olat']]
            t = fp['t']
            target = fp['target']
            obs = fp['sc']
            r.append(psoa.pointres(ifov, srfpoint, t, target, obs))
        res = np.mean(r)
        return res

    def roicoverage(self, roi, fplist):

        """
         Provided a list of footprints on a body surface and ROIs, this function
         computes the cumulative coverage and overlap on each ROI, in [%]
         Overlap accounts for the % of surface covered at least more than once
         [Note:] This needs a revision, only works for a single roi!!

         Programmers:  Paula Betriu (UPC/ESEIAAT)
         Date:         12/2023
         Revision:     1

         Usage:        [coverage, overlap] = roicoverage(target, roi, fplist)

         Inputs:
           > target:   SPICE ID (int) or name (string) of the body
           > roi:      roi struct that contains, at least, its boundary vertices,
                       in latitudinal coordinates (in [deg])
           > fplist:   list of footprints (struct). See footprint for further
                       information

         Output:
           > coverage: Percentage of surface that the list of provided footprints
                       collectively cover with respect to the ROI surface,
                       in [%]
           > overlap:  Percentage of overlap of the list of provided footprints
                       with respect to the ROI surface, in [%]
        """

        # Pre-allocate variables

        coverage = 0  # initialize coverage
        overlap = 0  # initialize overlap

        if not fplist:
            return coverage, overlap

        # In case roi is not input as a struct but as a matrix (vertices)...
        if not isinstance(roi, dict):
            col1, col2 = self.amsplit(roi.vertices[:, 0], roi.vertices[:, 1])
            vertices = np.hstack((col1.reshape(len(col1), 1), col2.reshape(len(col2), 1)))
            del roi
            roi = {'vertices': self.interppolygon(vertices)}

        # Cumulative coverage
        for footprint in fplist:
            if 'polyJ' not in roi:
                roi['polyJ'] = Polygon()
                roi['polyO'] = Polygon()

            # Footprint polyshape
            if len(footprint['bvertices']) == 0:
                continue

            x, y = footprint['bvertices'][:, 0], footprint['bvertices'][:, 1]

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
                polyFP = MultiPolygon(polygon_list)
            else:
                polyFP = Polygon(list(zip(x,y)))
            polyFP = polyFP.buffer(0)

            # Intersect footprint to cumulative coverage
            polyI = roi['polyJ'].intersection(polyFP)
            roi['polyO'] = roi['polyO'].union(polyI)

            # Join footprint
            roi['polyJ'] = roi['polyJ'].union(polyFP)

        # Calculate ROI specific coverage
        # Get total surface of the ROI
        lon, lat = roi['vertices'][:, 0], roi['vertices'][:, 1]
        # RA, _ = regionarea(target, lon, lat)

        if (np.isnan(lon)).any():
            nanindex = np.where(np.isnan(lon))[0]
            polygon_list = []
            for i in range(len(nanindex)):
                if i == 0:
                    polygon_list.append(Polygon(list(zip(lon[:nanindex[0]], lat[:nanindex[0]]))))
                else:
                    polygon_list.append(Polygon(
                        list(zip(lon[nanindex[i - 1] + 1:nanindex[i]], lat[nanindex[i - 1] + 1:nanindex[i]]))))
            if ~ np.isnan(lon[-1]):
                polygon_list.append(Polygon(list(zip(lon[nanindex[-1] + 1:], lat[nanindex[-1] + 1:]))))
            polyROI = MultiPolygon(polygon_list)
        else:
            polyROI = Polygon(list(zip(lon, lat)))
        polyROI = polyROI.buffer(0)

        # Intersect ROI with cumulative footprints
        if (isinstance(roi['polyJ'], Polygon) and roi['polyJ'].exterior.coords) or isinstance(roi['polyJ'], MultiPolygon):
            polyI = polyROI.intersection(roi['polyJ'])
            if (isinstance(polyI, Polygon) and polyI.exterior.coords) or isinstance(polyI, MultiPolygon):
                inter = polyROI.difference(polyI)
                coverage = ((polyROI.area - inter.area) * 100) / polyROI.area

        # Intersect ROI with overlap
        if (isinstance(roi['polyO'], Polygon) and roi['polyO'].exterior.coords) or isinstance(roi['polyO'], MultiPolygon):
            polyI = polyROI.intersection(roi['polyO'])
            if (isinstance(polyI, Polygon) and polyI.exterior.coords) or isinstance(polyI, MultiPolygon):
                inter = polyROI.difference(polyI)
                overlap = ((polyROI.area - inter.area) * 100) / polyROI.area

        return coverage, overlap

    def amsplit(self,x_, y_):
        """
        Provided the vertices of a polygon in latitudinal coordinates, this
        function analyzes if the polygon intercepts the anti-meridian line and,
        in that case, divides the polygon by this line.

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         10/2022

        Usage:        [xf, yf] = amsplit(x, y)

        Inputs:
          > x:        array of longitude values in [deg]. x ∈ [-180, 180]
          > y:        array of latitude values in [deg]. y ∈ [-90, 90]

        Outputs:
          > x:        array of longitude values in [deg]. In case the polygon
                      intercepts the a.m. line, then the longitude values of the
                      polygon are separated by a NaN value. x ∈ [-180, 180]
          > y:        array of latitude values in [deg]. In case the polygon
                      intercepts the a.m. line, then the latitude values of the
                      polygon are separated by a NaN value. y ∈ [-90, 90]
        """
        # [To be resolved]: If polygon longitude size is 180º, the function does
        # nothing.
        x = copy.deepcopy(x_)
        y = copy.deepcopy(y_)
        if (np.max(x) - np.min(x)) == 180:
            print("Full longitude. This function does nothing, in this case." +
                  "The user must check if the polygon is well defined!")

        # Previous check...
        if (np.max(x) - np.min(x)) <= 180:
            xf = x
            yf = y
            return xf, yf

        # In case the polygon indeed intercepts this line, then we're going to
        # calculate the intercept points by computing the intersection between two
        # polygons: the input polygon and another one which corresponds to the
        # anti-meridian line (actually, it's a small polygon because the 'intersect'
        # function does not operate well with lines)
        x[x < 0] += 360
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
            poly1 = MultiPolygon(polygon_list)
        else:
            poly1 = Polygon(zip(x, y))

        poly1 = poly1.buffer(0)

        vpoly2 = np.vstack((np.column_stack((180. * np.ones(20), np.linspace(-90., 90., 20))),
                            np.column_stack((181. * np.ones(20), np.linspace(90., -90., 20)))))
        poly2 = Polygon(vpoly2)

        # Compute the intersection points
        polyinter = (poly1.intersection(poly2)).buffer(0)

        if isinstance(polyinter, Polygon):
            xinter, yinter = np.array(polyinter.exterior.coords.xy)
        elif isinstance(polyinter, MultiPolygon):
            for i in range(len(polyinter.geoms)):
                xinteraux, yinteraux = np.array(polyinter.geoms[i].exterior.coords.xy)
                if i == 0:
                    xinter = np.append(xinteraux, np.nan)
                    yinter = np.append(yinteraux, np.nan)
                else:
                    xinter = np.append(xinter, np.append(xinteraux, np.nan))
                    yinter = np.append(yinter, np.append(yinteraux, np.nan))
            xinter = xinter[:-1]
            yinter = yinter[:-1]

        # Only keep the anti-meridian intercepts
        yi = yinter[np.abs(xinter - 180) < 1e-2]
        xi = 180 * np.ones(len(yi))

        # Define the new polygon
        P = np.column_stack((x, y))

        # Add the intersection points to the polygon vertices
        P = np.vstack((P, np.column_stack((xi, yi))))

        # Split the polygon in two, cleaved by the anti-meridian line
        x = np.sort(P[:, 0])
        P = P[np.argsort(P[:, 0])]
        poly1 = P[x >= 180, :]  # this is the polygon that falls in negative lon.
        poly2 = P[x <= 180, :]  # this is the polygon that falls in positive lon.

        # Sort the polygon vertices in clockwise order...
        poly1[:, 0], poly1[:, 1] = self.sortcw(poly1[:, 0], poly1[:, 1])
        poly2[:, 0], poly2[:, 1] = self.sortcw(poly2[:, 0], poly2[:, 1])

        # Retrieve the original values (longitude values cannot be >180º in our system) and output the final vertices
        xf = poly1[:, 0] - 360
        xf = np.append(xf, [np.nan])
        xf = np.append(xf, poly2[:, 0])

        yf = poly1[:, 1]
        yf = np.append(yf, [np.nan])
        yf = np.append(yf, poly2[:, 1])

        # [Future work]: Consider polygons that are split by the a.m. but also are
        # discontinuous in latitude (i.e., more than two polygons)
        # yi = sorted(yi, reverse=True)
        # if len(xi) / 2 > 1:
        #     for j in range(len(yi)):
        #         # Get inner intersections
        #         if yi[j] > min(yi) and yi[j] < max(yi):
        #             # To determine whether the intersection should be included in the left or
        #             # right sides of the a.m. we need to check the immediately upper and
        #             # lower points
        #             poly3 = Polygon()
        #             if Polygon(poly2).contains(Point(180, yi[j] + 1)):
        #                 # Split poly1
        #                 ind = (x >= 180) & (y >= yi[j])  # poly1[:, 1] >= yi[j]
        #                 spoly1 = sortcw(poly1[ind, 0], poly1[ind, 1])
        #                 poly3 = addboundary(poly3, spoly2)
        #             else:
        #                 # Split poly2
        #                 pass

        return xf, yf

    def interppolygon(self,roi0):
        """
        This function interpolates a polygon defined by longitude and latitude
        points. The interpolation distance is defined according to the minimum
        Euclidean distance between the points that enclose the polygon

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         1/2024

        Usage:        roi = interppolygon(roi0)

        Input:
          > roi0:     A Nx2 matrix where each row represents a point in 2D space,
                      typically [longitude, latitude]

        Output:
          > roi:      Updated roi with interpolated coordinates
        """

        # Definition of maximum allowable distance
        # Initialize variables to store longitude and latitude from the polygon.
        lon = roi0[:, 0]
        lat = roi0[:, 1]
        epsilon = np.inf  # initialize epsilon to infinity. This will be used to find the minimum distance between points.

        # Loop through each pair of points to find the minimum non-zero distance
        for i in range(len(lon) - 1):
            p1 = np.array([lon[i], lat[i]])  # current point
            p2 = np.array([lon[i + 1], lat[i + 1]])  # next point
            dist = np.linalg.norm(p2 - p1)  # calculate Euclidean distance between points
            if dist == 0:
                continue  # skip if points are identical
            if dist < epsilon:
                epsilon = dist  # update minimum distance

        # Find number of regions of polygon
        indnan = np.where(np.isnan(roi0[:, 0]))[0]
        newlat = []
        newlon = []
        if len(indnan) > 0:
            # indnan[-1] =  roi0.shape[0]
            indnan = np.append(indnan, roi0.shape[0])
            from_idx = 0
            for i in range(len(indnan)):
                to = indnan[i]
                # Perform great circles interpolation of the latitude and longitude based
                # on the minimum distance found
                latd = lat[from_idx:to]
                lond = lon[from_idx:to]
                auxlat, auxlon = self.interpm(latd, lond, math.ceil(epsilon / 2), 'gc')
                newlat.extend(auxlat)
                newlat.append(np.nan)
                newlon.extend(auxlon)
                newlon.append(np.nan)
                from_idx = to + 1
        else:
            # Perform great circles interpolation of the latitude and longitude based
            # on the minimum distance found
            newlat, newlon = self.interpm(lat, lon, math.ceil(epsilon / 2), 'gc')

        # Update the roi array with the interpolated latitude and longitude values
        roi = np.column_stack((newlon, newlat))

        # Sort coordinates in clockwise order
        # [roi[:, 0], roi[:, 1]] = sortcw(roi[:, 0], roi[:, 1])
        return roi

    def interpm(self,lat, lon, maxdiff, method='gc'):

        maxdiff = maxdiff * (math.pi / 180) * 6371.009

        if method != 'gc':
            raise ValueError("This example can only handle 'gc' interpolation (Great Circle).")

        latout = [lat[0]]
        lonout = [lon[0]]

        for i in range(1, len(lat)):
            start = Point(lat[i - 1], lon[i - 1])
            end = Point(lat[i], lon[i])

            dist = great_circle(start, end).kilometers

            if dist <= maxdiff:
                latout.append(lat[i])
                lonout.append(lon[i])
                continue

            num_points = int(dist // maxdiff)
            if round(dist % maxdiff) == 0 and dist % maxdiff <= 0.05:
                num_points = num_points - 1

            for j in range(1, num_points + 1):
                fraction = j / (num_points + 1)
                azimuth = self.calculate_azimuth(start, end)
                intermediate_point = great_circle(kilometers=fraction * dist).destination(
                    (start.latitude, start.longitude), azimuth)
                latout.append(intermediate_point.latitude)
                lonout.append(intermediate_point.longitude)

            latout.append(lat[i])
            lonout.append(lon[i])

        return latout, lonout

    def computeVisibleROI(self,roi_, et, target, obs):
        """
            Given a target area (region-of-interest) on a planetary body surface and
            an observer, this function calculates the portion of the former that is
            visible.
            Note: if roi intercepts the anti-meridian line, this function returns
            the visible roi polygon divided by this line (see amsplit function).

            Programmers:  Paula Betriu (UPC/ESEIAAT)
                          Diego Andía  (UPC/ESEIAAT)
            Date:         06/2023

            Usage:        vroi = visibleroi(roi, et, target, obs)

            Inputs:
              > roi:     matrix containing the vertices of the ROI polygon. The
                         vertex points are expressed in 2D, in latitudinal
                         coordinates [º]
                  # roi[:,1] correspond to the x values of the vertices
                  # roi[:,2] correspond to the y values of the vertices
              > et:      limb projection time, in seconds past J2000 epoch
              > target:  string name of the target body (SPICE ID)
              > obs:     string name of the observer (SPICE ID)

            Output:
              > vroi:    matrix containing the vertices of the intersection between
                         the input ROI polygon and the limb projection, i.e., the
                         visible portion of the ROI on the body surface as seen from
                         the observer
            """
        roi = copy.deepcopy(roi_)
        # Previous anti-meridian intersection check...
        ind1 = np.where(np.diff(np.sort(roi[:, 0])) >= 180)[0]  # find the discontinuity index
        if ind1.size > 0:
            col1, col2 = self.amsplit(roi[:, 0], roi[:, 1])
            roi = np.hstack((col1.reshape(len(col1), 1), col2.reshape(len(col2), 1)))

        # Parameters for mat2py_limbpt function
        flag = False  # assume the target area is visible from the instrument
        method = 'TANGENT/ELLIPSOID'
        _, targetframe, _ = M2P.mat2py_cnmfrm(target)  # body-fixed frame
        abcorr = 'XLT+S'
        corloc = 'CENTER'
        refvec = np.array([0, 0, 1])  # first of the sequence of cutting half-planes
        ncuts = int(1e3)  # number of cutting half-planes
        delrol = M2P.mat2py_twopi() / ncuts  # angular step by which to roll the
        # cutting half-planes about the observer-target vector
        schstp = 1.0e-4  # search angular step size
        soltol = 1.0e-7  # solution convergence tolerance

        # Limb calculation with mat2py_limbpt function
        _, limb, _, _ = M2P.mat2py_limbpt(method, target, et, targetframe, abcorr,
                                      corloc, obs, refvec, delrol, ncuts, schstp, soltol,
                                      ncuts)  # limb points expressed in targetframe ref frame
        _, lblon, lblat = M2P.mat2py_reclat(limb)  # conversion from rectangular to latitudinal coordinates
        lblon = lblon * M2P.mat2py_dpr()
        lblat = lblat * M2P.mat2py_dpr()

        ind2 = np.where(np.diff(np.sort(lblon)) >= 180)[0]
        if ind2.size > 0:
            lblon, lblat = self.amsplit(lblon, lblat)

        # We need to discern between two different limb:
        # 1.- Sub-spacecraft point is located at equator (observer-to-pole line is
        # perpendicular to normal vector at the poles). In this case, limb's
        # longitude cannot be > 180º
        # 2.- Sub-spacecraft point is not located at equator. In this case, limb's
        # longitude may be > 180º (and includes the north/south poles).

        northpole = False
        southpole = False
        # Check north-pole
        srfpoint = np.array([0, 90])
        angle = psoa.emissionang(srfpoint, et, target, obs)
        if angle < 90:
            northpole = True
        # Check south-pole
        srfpoint = np.array([0, -90])
        angle = psoa.emissionang(srfpoint, et, target, obs)
        if angle < 90:
            southpole = True

        # Case 1
        if not northpole and not southpole:
            # Check a.m. split
            ind2 = np.where(np.diff(np.sort(lblon)) >= 180)[0]  # find the discontinuity
            if ind2.size > 0:
                lblon, lblat = self.amsplit(lblon, lblat)
            # Check if we are keeping the correct polygon (full disk polygons may be
            # misleading, we can only guarantee through emission angle check)
            ext = False

            while not ext:
                randPoint = np.array([np.random.randint(-180, 181), np.random.randint(-90, 91)])
                point = Point(randPoint)
                if (np.isnan(lblon)).any():
                    nanindex = np.where(np.isnan(lblon))[0]
                    polygon_list = []
                    for i in range(len(nanindex)):
                        if i == 0:
                            polygon_list.append(Polygon(list(zip(lblon[:nanindex[0]], lblat[:nanindex[0]]))))
                        else:
                            polygon_list.append(Polygon(
                                list(zip(lblon[nanindex[i - 1] + 1:nanindex[i]],
                                         lblat[nanindex[i - 1] + 1:nanindex[i]]))))
                    if ~ np.isnan(lblon[-1]):
                        polygon_list.append(Polygon(list(zip(lblon[nanindex[-1] + 1:], lblat[nanindex[-1] + 1:]))))
                    polyaux = MultiPolygon(polygon_list)
                else:
                    polyaux = Polygon((list(zip(lblon, lblat))))
                polyaux = polyaux.buffer(0)

                if polyaux.intersects(point):
                    angle = psoa.emissionang(randPoint, et, target, obs)
                    if angle < 85:
                        ext = True
                else:
                    angle = psoa.emissionang(randPoint, et, target, obs)
                    if angle < 85:
                        ext = True
                        # This calculation is approximated, we should find a better way
                        # to find the complementary
                        # [Future work]
                        lonmap = [-180, -180, 180, 180]
                        latmap = [-90, 90, 90, -90]
                        polymap = Polygon(list(zip(lonmap, latmap)))
                        if (np.isnan(lblon)).any():
                            nanindex = np.where(np.isnan(lblon))[0]
                            polygon_list = []
                            for i in range(len(nanindex)):
                                if i == 0:
                                    polygon_list.append(Polygon(list(zip(lblon[:nanindex[0]], lblat[:nanindex[0]]))))
                                else:
                                    polygon_list.append(Polygon(
                                        list(zip(lblon[nanindex[i - 1] + 1:nanindex[i]],
                                                 lblat[nanindex[i - 1] + 1:nanindex[i]]))))
                            if ~ np.isnan(lblon[-1]):
                                polygon_list.append(
                                    Polygon(list(zip(lblon[nanindex[-1] + 1:], lblat[nanindex[-1] + 1:]))))
                            poly1 = MultiPolygon(polygon_list)
                        else:
                            poly1 = Polygon((list(zip(lblon, lblat))))
                        poly1 = poly1.buffer(0)
                        poly1 = polymap.difference(poly1)
                        poly1 = poly1.buffer(0)

                        if isinstance(poly1, Polygon):
                            lblon, lblat = np.array(poly1.exterior.coords.xy)
                        elif isinstance(poly1, MultiPolygon):
                            for i in range(len(poly1.geoms)):
                                lblonaux, lblataux = np.array(poly1.geoms[i].exterior.coords.xy)
                                if i == 0:
                                    lblon = np.append(lblonaux, np.nan)
                                    lblat = np.append(lblataux, np.nan)
                                else:
                                    lblon = np.append(lblon, np.append(lblonaux, np.nan))
                                    lblat = np.append(lblat, np.append(lblataux, np.nan))
                            lblon = lblon[:-1]
                            lblat = lblat[:-1]
        else:
            # Case 2.
            lblon, indsort = np.sort(lblon), np.argsort(lblon)
            lblat = lblat[indsort]

            if northpole or southpole:
                # Include northpole to close polygon
                auxlon, auxlat = copy.deepcopy(lblon), copy.deepcopy(lblat)
                lblon = np.zeros(len(auxlon) + 2)
                lblat = np.zeros(len(auxlat) + 2)
                if northpole:
                    lblon[0], lblat[0] = -180, 90
                    lblon[-1], lblat[-1] = 180, 90
                else:
                    lblon[0], lblat[0] = -180, -90
                    lblon[-1], lblat[-1] = 180, -90
                lblon[1:-1] = auxlon
                lblat[1:-1] = auxlat

        # roi and limb intersection
        if (np.isnan(lblon)).any():
            nanindex = np.where(np.isnan(lblon))[0]
            polygon_list = []
            for i in range(len(nanindex)):
                if i == 0:
                    polygon_list.append(Polygon(list(zip(lblon[:nanindex[0]], lblat[:nanindex[0]]))))
                else:
                    polygon_list.append(Polygon(
                        list(zip(lblon[nanindex[i - 1] + 1:nanindex[i]], lblat[nanindex[i - 1] + 1:nanindex[i]]))))
            if ~ np.isnan(lblon[-1]):
                polygon_list.append(Polygon(list(zip(lblon[nanindex[-1] + 1:], lblat[nanindex[-1] + 1:]))))
            poly1 = MultiPolygon(polygon_list)
        else:
            poly1 = Polygon((list(zip(lblon, lblat))))
        poly1 = poly1.buffer(0)

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
            poly2 = MultiPolygon(polygon_list)
        else:
            poly2 = Polygon((list(zip(roi[:, 0], roi[:, 1]))))
        poly2 = poly2.buffer(0)

        inter = poly1.intersection(poly2)
        inter = inter.buffer(0)

        # output visible roi
        if isinstance(inter, Polygon):
            vroi = np.array(inter.exterior.coords)
        elif isinstance(inter, MultiPolygon):
            for i in range(len(inter.geoms)):
                if i == 0:
                    vroi = np.vstack((np.array(inter.geoms[i].exterior.coords), [np.nan, np.nan]))
                else:
                    vroi = np.vstack((vroi, np.array(inter.geoms[i].exterior.coords), [np.nan, np.nan]))
            vroi = vroi[:-1, :]

        # visibility flag
        if vroi.size == 0:
            flag = True

        return vroi, inter, flag

    @staticmethod
    def sortcw(*args):
        """
        Given a set of polygon vertices, this function sorts them in clockwise order.

        Programmers: Paula Betriu (UPC/ESEIAAT)
        Date: 10/2022

        Usage:
            x, y = sortcw(x, y)
            x, y, z = sortcw(x, y, z)

        Inputs:
        - x:    x rectangular coordinate (for 3 inputs) or longitude values (for 2 inputs).
                Units are irrelevant as long as the 2 or 3 arrays are consistent
        - y:    y rectangular coordinate (for 3 inputs) or latitude values (for 2 inputs).
                Units are irrelevant as long as the 2 or 3 arrays are consistent
        - z:    z rectangular coordinate (for 3 inputs)

        Returns:
        - Sorted inputs in clockwise order

        Note: 2D sorting algorithm does not work with concave algorithms. In such
        case, 3D rectangular coordinates are recommended.
        """
        if len(args) < 3:
            #  This algorithm consists of dividing the space in 4 quadrants,
            #  centered at the polygon centroid (any other inner point would be also
            #  valid). We calculate iteratively the angle with respect to that point
            #  and sort the vertices according to their respective angle value.
            #  This algorithm does not work with non-convex polygons
            x = copy.deepcopy(args[0])
            y = copy.deepcopy(args[1])

            if np.size(x) < 3:
                return x, y

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

            if polygon.is_empty:
                unique_points = list(set(list(zip(x, y))))
                xuni, yuni = zip(*unique_points)
                if (np.isnan(xuni)).any():
                    nanindex = np.where(np.isnan(xuni))[0]
                    polygon_list = []
                    for i in range(len(nanindex)):
                        if i == 0:
                            polygon_list.append(Polygon(list(zip(xuni[:nanindex[0]], yuni[:nanindex[0]]))))
                        else:
                            polygon_list.append(Polygon(
                                list(
                                    zip(xuni[nanindex[i - 1] + 1:nanindex[i]], yuni[nanindex[i - 1] + 1:nanindex[i]]))))
                    if ~ np.isnan(xuni[-1]):
                        polygon_list.append(Polygon(list(zip(xuni[nanindex[-1] + 1:], yuni[nanindex[-1] + 1:]))))
                    polygon = MultiPolygon(polygon_list)
                else:
                    polygon = Polygon((list(zip(xuni, yuni))))

            cx = polygon.centroid.x
            cy = polygon.centroid.y  # polygon centroid

            angle = np.arctan2(np.array(y) - cy, np.array(x) - cx)  # obtain angle
            ind = np.argsort(angle)[::-1]  # sort angles and get the indices
            x_sorted = np.array(x)[ind]  # sort the longitude values according to the angle order
            y_sorted = np.array(y)[ind]  # sort the latitude values according to the angle order
            return x_sorted, y_sorted  # save output

        elif len(args) == 3:
            # Algorithm extracted from
            # https://stackoverflow.com/questions/47949485/
            # sorting-a-list-of-3d-points-in-clockwise-order
            x = copy.deepcopy(args[0])
            y = copy.deepcopy(args[1])
            z = copy.deepcopy(args[2])

            innerpoint = np.array(
                [np.mean(np.array(x)), np.mean(np.array(y)), np.mean(np.array(z))])  # find an inner point (this
            # is not the centroid but still works)

            i, j, k = np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1])

            # Compute the cross products (perpendicular to the surface normal)
            pn = [np.cross(i, innerpoint), np.cross(j, innerpoint), np.cross(k, innerpoint)]

            # Choose the largest cross product to get vector p
            p = max(pn, key=np.linalg.norm)

            # Compute a second perpendicular vector q in the plane
            q = np.cross(innerpoint, p)

            # Now we have two perpendicular reference vectors in the plane given by
            # the normal. Take triple products of those, and these will be the sine
            # and the cosine of an angle that can be used for sorting
            angles = []
            for xi, yi, zi in zip(x, y, z):
                rmc = np.array([xi, yi, zi]) - innerpoint
                t = np.dot(innerpoint, np.cross(rmc, p))
                u = np.dot(innerpoint, np.cross(rmc, q))
                angle = np.arctan2(u, t)
                angles.append(angle)

            angles = np.array(angles)
            # Sort vertices by angles
            ind = np.argsort(angles)
            x_sorted = np.array(x)[ind]
            y_sorted = np.array(y)[ind]
            z_sorted = np.array(z)[ind]
            return x_sorted, y_sorted, z_sorted
        else:
            raise ValueError("Too many input arguments")

    @staticmethod
    def calculate_azimuth(start, end):
        lat1 = math.radians(start.latitude)
        lon1 = math.radians(start.longitude)
        lat2 = math.radians(end.latitude)
        lon2 = math.radians(end.longitude)

        dlon = lon2 - lon1
        x = math.sin(dlon) * math.cos(lat2)
        y = math.cos(lat1) * math.sin(lat2) - (math.sin(lat1) * math.cos(lat2) * math.cos(dlon))
        initial_bearing = math.atan2(x, y)
        initial_bearing = math.degrees(initial_bearing)
        compass_bearing = (initial_bearing + 360) % 360

        return compass_bearing

    @staticmethod
    def inFOVprojection(boundPoints, N, surfPoints, count):
        """
        The FOV projection is enclosed in the target surface.
        """
        # Close polygon
        boundPoints = np.hstack((boundPoints, boundPoints[:, [0]]))
        # count = 0
        surfPoints = np.zeros([N * (max(boundPoints.shape) - 1), 3])
        # high resolution
        for i in range(max(boundPoints.shape) - 1):
            # linear (approximation) interpolation between vertices to define
            # the boundary of the footprint
            v = boundPoints[:, i + 1] - boundPoints[:, i]
            lambda_vals = np.linspace(0, 1, N)  # line parametrization
            for l in range(N):
                surfPoints[count, 0] = boundPoints[0, i] + v[0] * lambda_vals[l]
                surfPoints[count, 1] = boundPoints[1, i] + v[1] * lambda_vals[l]
                surfPoints[count, 2] = boundPoints[2, i] + v[2] * lambda_vals[l]
                count += 1
        return boundPoints, N, surfPoints, count

    @staticmethod
    def groundtrack(obs, t, target):
        """
        This function returns the spacecraft ground track across the target
        surface, at time t.

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         09/2023

        Usage:        gtlon, gtlat = groundtrack(obs, t, target)

        Inputs:
          > obs:      string SPICE name of the observer body
          > t:        time epoch in TDB seconds past J2000 epoch. It can be
                      either a single point in time or a discretized vector of
                      different time values
          > target:   string SPICE name of the target body

        Outputs:
          > gtlon:    longitude coordinate of the observer ground track, in
                      [deg]
          > gtlat:    latitude coordinate of the observer ground track, in
                      [deg]
        """
        # Pre-allocate variables
        method = 'INTERCEPT/ELLIPSOID'  # target modeling
        abcorr = 'NONE'  # aberration correction
        _, tframe, _ = M2P.mat2py_cnmfrm(target)  # body-fixed frame

        # Get ground track
        sctrack, _, _ = M2P.mat2py_subpnt(method, target, t, tframe, abcorr, obs)  # sub-spacecraft point
        _, gtlon, gtlat = M2P.mat2py_reclat(sctrack)  # convert to latitudinal coordinates

        gtlon = gtlon * M2P.mat2py_dpr()  # [rad] to [deg]
        gtlat = gtlat * M2P.mat2py_dpr()  # [rad] to [deg]

        return gtlon, gtlat

    def closestSide(gt1, gt2, targetArea, angle):
        """
        Given a region-of-interest, this function calculates the spacecraft
        ground track position with respect to the edges of the target area. It
        also determines the direction in which the spacecraft is moving relative
        to the target area. The target area is rotated based on the specified
        angle to align with the instrument's observation footprint.

        Inputs:
          > gt1:          initial ground track position ([lon, lat]) of the
                          spacecraft, in [deg]
          > gt2:          subsequent ground track position ([lon, lat]) of the
                          spacecraft, in [deg]
          > targetArea:   matrix containing the vertices of the target area. The
                          vertex points are expressed in 2D latitudinal coord.
          > angle:        rotation angle of the observation footprint axes
                          relative to the target area, in [deg]

        Outputs:
          > dir1:         closest side of the target area to the spacecraft's
                          initial position ('north', 'south', 'east', 'west')
          > dir2:         direction of the spacecraft's movement relative to the
                          target area ('north', 'south', 'east', 'west')
        """

        # Rotate target area according to the footprint's angle
        angle = -angle * M2P.mat2py_rpd()
        rotmat = np.array([[math.cos(angle), -math.sin(angle)],
                           [math.sin(angle), math.cos(angle)]])
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
            poly_aux = MultiPolygon(polygon_list)
        else:
            poly_aux = Polygon((list(zip(targetArea[:, 0], targetArea[:, 1]))))

        cx, cy = poly_aux.centroid.x, poly_aux.centroid.y

        roi = np.zeros((max(np.shape(targetArea)), 2))
        for j in range(max(np.shape(targetArea))):
            roi[j, :] = np.dot(rotmat, (targetArea[j] - np.array([cx, cy]))) + np.array([cx, cy])

        # Adjust ground track position for spacecraft movement analysis with respect to the oriented area
        # Assumption: Tri-axial ellipsoid to model the target surface
        # Initial position
        sclon, sclat = np.dot(rotmat, ((np.array(gt1)).reshape(2, ) - np.array([cx, cy]))) + np.array([cx, cy])
        # Subsequent position
        sclon_, sclat_ = np.dot(rotmat, ((np.array(gt2)).reshape(2, ) - np.array([cx, cy]))) + np.array([cx, cy])

        roiaux = roi[~np.isnan(roi).all(axis=1)]

        # Find the 4 boundary vertices of the rotated target area
        maxlon, minlon = np.max(roiaux[:, 0]), np.min(roiaux[:, 0])
        maxlat, minlat = np.max(roiaux[:, 1]), np.min(roiaux[:, 1])

        # ROI's boundary box
        xlimit = [minlon, maxlon]
        ylimit = [minlat, maxlat]
        xbox = np.array([xlimit[0], xlimit[0], xlimit[1], xlimit[1], xlimit[0]])
        ybox = np.array([ylimit[0], ylimit[1], ylimit[1], ylimit[0], ylimit[0]])
        boundary_box = Polygon(zip(xbox, ybox)).buffer(0)

        # Define line between the centroid and the ground track
        line = LineString([(sclon, sclat), (cx, cy)])
        intersection = line.intersection(boundary_box)  # check if the line intersects with the boundary box to determine
        # closest side

        if not intersection.is_empty:
            if intersection.geom_type == 'Point':
                xi, yi = intersection.x, intersection.y
            elif intersection.geom_type == 'MultiPoint':
                xi, yi = zip(*[(pt.x, pt.y) for pt in intersection])
                # if math.sqrt((xi[-1] - xi[0]) ** 2 + (yi[-1] - yi[0]) ** 2) < 0.026:
                xi = xi[0]
                yi = yi[0]

            elif intersection.geom_type == 'LineString':
                xi, yi = zip(*intersection.coords)
                # if math.sqrt((xi[1] - xi[0]) ** 2 + (yi[1] - yi[0]) ** 2) < 0.026:
                xi = xi[0]
                yi = yi[0]

        # Determine closest side based on the intersection points
        if intersection.is_empty:  # The ground track is inside the ROI's boundary box (no intersection)
            # Define corner points of the boundary box
            p1 = np.array([maxlon, maxlat])
            p2 = np.array([maxlon, minlat])
            p3 = np.array([minlon, minlat])
            p4 = np.array([minlon, maxlat])

            # Calculate distance to the mid-points of the 4 edges of the boundary box
            midp = [
                0.5 * (p1 + p4),  # midup
                0.5 * (p2 + p3),  # middown
                0.5 * (p3 + p4),  # midleft
                0.5 * (p1 + p2)  # midright
            ]

            # Find minimum distance to midpoints
            mindist = float('inf')
            for i, midpoint in enumerate(midp):
                dist = np.linalg.norm(np.array([sclon, sclat]) - midpoint)
                if dist < mindist:
                    mindist = dist
                    minidx = i

            # Determine closest side
            if minidx == 1:
                dir1 = 'north'
            elif minidx == 2:
                dir1 = 'south'
            elif minidx == 3:
                dir1 = 'west'
            elif minidx == 4:
                dir1 = 'east'
        else:
            # Determine closest side based on the intersection point
            if abs(xi - maxlon) < 1e-2:
                dir1 = 'east'
            elif abs(xi - minlon) < 1e-2:
                dir1 = 'west'
            elif abs(yi - minlat) < 1e-2:
                dir1 = 'south'
            elif abs(yi - maxlat) < 1e-2:
                dir1 = 'north'

            # Determine if the spacecraft is moving towards or away from the cside
            dist = np.linalg.norm([sclon - cx, sclat - cy])
            dist_ = np.linalg.norm([sclon_ - cx, sclat_ - cy])
            if dist_ < dist:
                if dir1 == 'north':
                    dir1 = 'south'
                elif dir1 == 'south':
                    dir1 = 'north'
                elif dir1 == 'east':
                    dir1 = 'west'
                elif dir1 == 'west':
                    dir1 = 'east'

        # Calculate how the spacecraft ground track position is moving along the map.
        # The coverage path does not only depend on the spacecraft position itself but also its velocity direction
        # Determine direction of coverage path based on the spacecraft's movement

        # Determine direction of coverage path based on the spacecraft's movement
        if dir1 in {'north', 'south'}:  # Horizontal sweep
            if (sclon - sclon_) >= 0:  # sc is moving leftwards
                dir2 = 'west'  # if the spacecraft is moving left (in
                # the topography map) then the coverage path should start
                # at the position furthest to the right (right -> left
                # direction)
            else:  # sc is moving to the right
                dir2 = 'east'  # if the spacecraft is moving right (in
                # the topography map) then the coverage path should start
                # at the position furthest to the left (left -> right
                # direction)
        elif dir1 in {'east', 'west'}:  # Vertical sweep
            if (sclat - sclat_) < 0:  # sc is moving upwards
                dir2 = 'north'  # if the spacecraft is moving up (in the
                # topography map) then the coverage path should start at
                # the position furthest to the bottom (down -> top
                # direction)
            else:  # sc is moving down
                dir2 = 'south'  # if the spacecraft is moving down (in the
                # topography map) then the coverage path should start at
                # the position furthest to the top (top -> down direction)

        return dir1, dir2


    def grid2D(self,fpref, olapx, olapy, gamma_, targetArea):
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
                        list(zip(targetArea[nanindex[i - 1] + 1:nanindex[i], 0],
                                 targetArea[nanindex[i - 1] + 1:nanindex[i], 1]))))
            if ~ np.isnan(targetArea[-1, 0]):
                polygon_list.append(
                    Polygon(list(zip(targetArea[nanindex[-1] + 1:, 0], targetArea[nanindex[-1] + 1:, 1]))))
            polygon = MultiPolygon(polygon_list)
        else:
            polygon = Polygon(targetArea)

        cx, cy = polygon.centroid.x, polygon.centroid.y

        orientedArea = np.zeros([max(np.shape(targetArea)), 2])
        for j in range(max(np.shape(targetArea))):
            orientedArea[j, :] = np.array([cx, cy]) + rotmat @ (targetArea[j, :] - np.array([cx, cy]))

        gamma = np.array([cx, cy]) + rotmat @ (np.array(gamma) - np.array([cx, cy]))

        # If the area is divided in smaller regions, then we get the convex polygon
        # that encloses all of them (flood-fill)
        aux = np.column_stack(
            (orientedArea[~np.isnan(orientedArea[:, 0]), 0], orientedArea[~np.isnan(orientedArea[:, 1]), 1]))
        # convhull does not accept NaN nor Inf

        k = (ConvexHull(aux)).vertices  # boundary vertices that constitute the convex polygon
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
        gridPoints, _ = self.floodFillAlgorithm(fpref['width'], fpref['height'], olapx, olapy, gamma, orientedArea, periArea,
                                           np.array([]), np.array([]), '4fill')

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
            sortedGrid = np.array(
                sorted(gridPoints, key=lambda x: -x[1]))  # the elements of gridPoints are sorted by latitude (+ to -)
            uniqueLat = np.unique([pt[1] for pt in sortedGrid])  # get the different latitude values
            ind = np.abs(
                np.diff(uniqueLat)) < 1e-5  # double check that there are no "similar" latitude values (it may happen)
            ind = np.append(ind, False)
            uniqueLat = uniqueLat[~ind]
            uniqueLon = np.unique([pt[0] for pt in sortedGrid])  # get the different longitude values unique check
            ind = np.abs(
                np.diff(uniqueLon)) < 1e-5  # double check that there are no "similar" longitude values (it may happen)
            ind = np.append(ind, False)
            uniqueLon = uniqueLon[~ind]

            # Sort and rotate the grid points and insert them in the grid matrix
            matrixGrid = [[None for _ in range(max(np.shape(uniqueLon)))] for _ in range(max(np.shape(uniqueLat)))]
            for i in range(max(np.shape(uniqueLat))):
                # We will sweep across the grid by, first, latitude and, second, longitude
                lat = uniqueLat[max(np.shape(uniqueLat)) - 1 - i]
                indlat = np.abs(sortedGrid[:, 1] - lat) < 1e-5
                mrow = sortedGrid[indlat]
                mrow = np.sort(mrow[:, 0], axis=0)
                for j in range(max(np.shape(mrow))):
                    indlon = np.abs(uniqueLon - mrow[j]) < 1e-5
                    lon = mrow[j]
                    for k in range(max(np.shape(uniqueLon))):
                        if indlon[k]:
                            matrixGrid[i][k] = np.array([cx, cy]) + rotmat.T @ (
                                        np.array([lon, lat]) - np.array([cx, cy]))

        return matrixGrid, dirx, diry

    def floodFillAlgorithm(self,w, h, olapx, olapy, gamma, targetArea, perimeterArea, gridPoints_, vPoints_, method):
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

        if isinstance(gridPoints, np.ndarray):
            gridPoints = list(gridPoints)
        if isinstance(vPoints, np.ndarray):
            vPoints = list(vPoints)

        # Variables pre-allocation
        inside = False
        ovlapx = olapx * w / 100;
        ovlapy = olapy * h / 100  # convert overlaps from
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
                    polygon_list.append(
                        Polygon(list(zip(perimeterArea[:nanindex[0], 0], perimeterArea[:nanindex[0], 1]))))
                else:
                    polygon_list.append(Polygon(
                        list(zip(perimeterArea[nanindex[i - 1] + 1:nanindex[i], 0],
                                 perimeterArea[nanindex[i - 1] + 1:nanindex[i], 1]))))
            if ~ np.isnan(perimeterArea[-1, 0]):
                polygon_list.append(
                    Polygon(list(zip(perimeterArea[nanindex[-1] + 1:, 0], perimeterArea[nanindex[-1] + 1:, 1]))))
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
                polygon_list.append(
                    Polygon(list(zip(targetArea[nanindex[-1] + 1:, 0], targetArea[nanindex[-1] + 1:, 1]))))
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
            gridPoints, vPoints = self.floodFillAlgorithm(w, h, olapx, olapy,
                                                     np.array([gamma[0] - w + ovlapx, gamma[1]]),
                                                     targetArea, perimeterArea,
                                                     gridPoints, vPoints, method)
            # South
            gridPoints, vPoints = self.floodFillAlgorithm(w, h, olapx, olapy,
                                                     np.array([gamma[0], gamma[1] - h + ovlapy]),
                                                     targetArea, perimeterArea,
                                                     gridPoints, vPoints, method)
            # North
            gridPoints, vPoints = self.floodFillAlgorithm(w, h, olapx, olapy,
                                                     np.array([gamma[0], gamma[1] + h - ovlapy]),
                                                     targetArea, perimeterArea,
                                                     gridPoints, vPoints, method)
            # East
            gridPoints, vPoints = self.floodFillAlgorithm(w, h, olapx, olapy,
                                                     np.array([gamma[0] + w - ovlapx, gamma[1]]),
                                                     targetArea, perimeterArea,
                                                     gridPoints, vPoints, method)

            # If method is '8fill', check diagonal neighbors
            if method == '8fill':
                # Northwest
                gridPoints, vPoints = self.floodFillAlgorithm(w, h, olapx, olapy,
                                                         np.array([gamma[0] - w + ovlapx, gamma[1] + h - ovlapy]),
                                                         targetArea, perimeterArea,
                                                         gridPoints, vPoints, method)
                # Southwest
                gridPoints, vPoints = self.floodFillAlgorithm(w, h, olapx, olapy,
                                                         np.array([gamma[0] - w + ovlapx, gamma[1] - h + ovlapy]),
                                                         targetArea, perimeterArea,
                                                         gridPoints, vPoints, method)
                # Northeast
                gridPoints, vPoints = self.floodFillAlgorithm(w, h, olapx, olapy,
                                                         np.array([gamma[0] + w - ovlapx, gamma[1] + h - ovlapy]),
                                                         targetArea, perimeterArea,
                                                         gridPoints, vPoints, method)
                # Southeast
                gridPoints, vPoints = self.floodFillAlgorithm(w, h, olapx, olapy,
                                                              np.array([gamma[0] + w - ovlapx, gamma[1] - h + ovlapy]),
                                                              targetArea, perimeterArea,
                                                              gridPoints, vPoints, method)

        return np.array(gridPoints), np.array(vPoints)

    @staticmethod
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
            tour = grid[0]

        sweep = dir2 in ['east', 'south']

        # Plan tour over the grid discretization
        # The origin of the coverage path depends on the spacecraft ground track position
        if dir1 in ['north', 'south']:  # Horizontal sweep
            if sweep:
                bearing = True  # left -> right
            else:
                bearing = False  # right -> left

            tour = [[] for _ in range(
                np.count_nonzero([item is not None for row in grid for item in row]))]  # list of planned observations
            ii = 0
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
                        x, y = grid[irow][icol][0], grid[irow][icol][1]
                        tour[ii] = np.array([x, y])  # Save it in the coverage tour
                        ii += 1

                bearing = not bearing  # Switch coverage direction after each row sweeping, i.e. left (highest lon) to right
                # (lowest lon) or vice versa

        elif dir1 in ['east', 'west']:  # Vertical sweep
            if sweep:
                bearing = True  # top -> down
            else:
                bearing = False  # down -> top

            tour = [[] for _ in range(
                np.count_nonzero([item is not None for row in grid for item in row]))]  # list of planned observations
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
    @staticmethod
    def grid2map(grid):
        """
        This function creates a map from a given grid. It adds a border of NaN
        values around the entire grid. Within the grid, any lists that are empty
        or have been excluded from observation (because their coverage is deemed
        too small) are also filled with NaN values.

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         09/2022

        Usage:        map = grid2map(grid)

        Input:
          > grid:          list of lists that contains a 2-element vector for grid
                           points within the ROI, or None for points outside
                           the ROI or excluded from coverage.

        Returns:
          > map:           list of lists representing the map. It is the input grid
                           with an added border of NaN values. Inside the grid,
                           lists with None are replaced with [NaN NaN]

        """

        # Initialize map with additional rows and columns to place the NaN boundaries
        map_height = len(grid) + 2
        map_width = len(grid[0]) + 2
        map = [[[] for _ in range(map_width)] for _ in range(map_height)]

        # Populate the map with data from the grid and NaN borders
        for i in range(map_height):
            for j in range(map_width):
                if i == 0 or j == 0 or i == map_height - 1 or j == map_width - 1 or grid[i - 1][j - 1] is None:
                    map[i][j] = np.array([np.nan, np.nan])
                else:
                    map[i][j] = grid[i - 1][j - 1]

        return map

    @staticmethod
    def getMapNeighbours(indrow, indcol, map, search='all'):
        """
        Given a grid of points (2D list) and an element, this function outputs the neighboring points
        of the current point in the matrix.

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         09/2022
        Last Rev.:    06/2023

        Usage:        n = getMapNeighbours(indrow, indcol, map)
                      n = getMapNeighbours(indrow, indcol, map, search)

        Inputs:
          > indrow:       int row index of the matrix element (grid point)
          > indcol:       int column index of the matrix element (grid point)
          > map:          2D list of grid points. In order to avoid
                          mapping boundaries, map is bounded by NaN rows and
                          columns (first and last)
          > search:       string that defines if the function shall differentiate
                          between 'cardinal' and 'diagonal' searches. Otherwise,
                          the 8 adjacent points are visited

        Returns:
          > n:            2D list with the non-NaN neighbouring points in
                          the map
        """

        # Previous checks...
        # Searching element is not in the boundaries
        if indrow == 0 or indrow == len(map) - 1 or indcol == 0 or indcol == len(map[0]) - 1:
            raise ValueError("Searching element cannot be in the map boundaries")
        # Future work: check if first and last rows and columns are NaN

        # Search neighbors of the given element in the map
        if search == 'all':
            aux_n = [None] * 8
            if map[indrow][indcol] is not None and not np.isnan(map[indrow][indcol]).any():
                aux_n[0] = map[indrow - 1][indcol + 1]  # northeast
                aux_n[1] = map[indrow][indcol + 1]  # east
                aux_n[2] = map[indrow + 1][indcol + 1]  # southeast
                aux_n[3] = map[indrow - 1][indcol]  # north
                aux_n[4] = map[indrow + 1][indcol]  # south
                aux_n[5] = map[indrow - 1][indcol - 1]  # northwest
                aux_n[6] = map[indrow][indcol - 1]  # west
                aux_n[7] = map[indrow + 1][indcol - 1]  # southwest

        elif search == 'cardinal':
            aux_n = [None] * 4
            aux_n[0] = map[indrow - 1][indcol]  # north
            aux_n[1] = map[indrow][indcol + 1]  # east
            aux_n[2] = map[indrow + 1][indcol]  # south
            aux_n[3] = map[indrow][indcol - 1]  # west

        elif search == 'diagonal':
            aux_n = [None] * 4
            aux_n[0] = map[indrow - 1][indcol - 1]  # northwest
            aux_n[1] = map[indrow - 1][indcol + 1]  # northeast
            aux_n[2] = map[indrow + 1][indcol + 1]  # southeast
            aux_n[3] = map[indrow + 1][indcol - 1]  # southwest

        # Output neighbours (not empy nor NaN)
        n = [neighbor for neighbor in aux_n if neighbor is not None and not np.isnan(neighbor).any()]

        return n

    @staticmethod
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
        n[1] = gamma + (w - ovlapx) * dx  # east
        n[2] = gamma + (w - ovlapx) * dx + (-h + ovlapy) * dy  # southeast
        n[3] = gamma + (h - ovlapy) * dy  # north
        n[4] = gamma + (-h + ovlapy) * dy  # south
        n[5] = gamma + (-w + ovlapx) * dx + (h - ovlapy) * dy  # northwest
        n[6] = gamma + (-w + ovlapx) * dx  # west
        n[7] = gamma + (-w + ovlapx) * dx + (-h + ovlapy) * dy  # southwest

        # Cardinal and diagonal neighbouring points grid indices
        nind = [None] * 8
        ind = np.array(ind)
        nind[0] = ind + np.array([-1, 1])
        nind[1] = ind + np.array([0, 1])
        nind[2] = ind + np.array([1, 1])
        nind[3] = ind + np.array([-1, 0])
        nind[4] = ind + np.array([1, 0])
        nind[5] = ind + np.array([-1, -1])
        nind[6] = ind + np.array([0, -1])
        nind[7] = ind + np.array([1, -1])

        return n, nind

    @staticmethod
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
    @staticmethod
    def insertTiles(*args):
        """
        This function includes new observation points in a planned tour
        observations

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         06/2023

        Usage:        [tour, map] = insertTiles(map, newp, indp)

        Inputs:
          > map:       list of lists of grid points. In order to avoid
                       mapping boundaries, map is bounded by NaN rows and
                       columns (first and last)
          > newp:      list of lists of the new observation points to be
                       included in 'tour'
          > indp:      list of lists of the new observation points locations to
                       be included in 'map'

        Returns:
          > map:       updated list of lists of grid points

        """

        map = args[0]
        newp = args[1]
        indp = args[2]

        # Insert elements in map
        offcol = 0
        offrow = 0

        for i in range(len(newp)):
            aux_map = copy.deepcopy(map)  # save the map at its current state (for potential relocation purposes)
            indel = indp[i]

            # Update index position
            indel[0] += offrow
            indel[1] += offcol

            # If the index is in the boundaries (or even further) of the map, then
            # we will have to relocate the elements in map and create a bigger grid

            offrow0 = copy.deepcopy(offrow)
            offcol0 = copy.deepcopy(offcol)

            if indel[0] >= (len(map) - 1):  # last row or more
                nrows = 2 + indel[0] - len(map)  # number of additional rows in the grid (last rows)

                # Map relocation
                map = aux_map + [[np.array([np.nan, np.nan]) for _ in range(len(map[0]))] for _ in range(nrows)]

            elif indel[0] <= 0:  # First row or less
                nrows = 1 - indel[0]  # Number of additional rows in the grid (first rows)
                offrow += nrows  # rows offset

                # Map relocation
                map = [[np.array([np.nan, np.nan]) for _ in range(len(map[0]))] for _ in range(nrows)] + aux_map

            aux_map = copy.deepcopy(map)

            if indel[1] >= (len(map[0]) - 1):  # last column or more
                ncols = 2 + indel[1] - len(map[0])  # number of additional columns
                # in the grid (last columns)

                # Map relocation
                map = [row + [np.array([np.nan, np.nan])] * ncols for row in aux_map]

            elif indel[1] <= 0:  # First column or less
                ncols = 1 - indel[1]  # number of additional columns in the grid (first columns)
                offcol += ncols  # columns offset

                # Map relocation
                map = [[np.array([np.nan, np.nan])] * ncols + row for row in aux_map]

            # Update the current element index
            indel[0] += (offrow - offrow0)
            indel[1] += (offcol - offcol0)

            # Include elements in map
            map[indel[0]][indel[1]] = newp[i]

        return map

    # The following is a MATLAB code that must be translated in Python. Since it is a comment, it is not used so it is not
    # translated for now.

    # [Future work]: # # Insert elements in tour?
    # switch heuristics
    #     case 'nearest neighbour'
    #
    #         # Find the nearest neighbour in the tour and insert the new
    #         # observation before or after this element
    #         for i=1:length(newp)
    #
    #             # Find the nearest element in the tour
    #             mindist = inf;
    #             for j=1:length(tour)
    #                 dist = norm(tour{j} - newp{i});
    #                 if dist < mindist
    #                     mindist = dist;
    #                     idx     = j;
    #                 end
    #             end
    #
    #             # Evaluate before and after elements distance to decide
    #             # whether the new point should be put before or after its
    #             # nearest neighbour
    #             aux_tour = tour;
    #             if idx ~= 1 && idx < length(tour)
    #                 dist1 = norm(tour{idx - 1} - newp{i});
    #                 dist2 = norm(tour{idx + 1} - newp{i});
    #                 if dist1 < dist2
    #                     tour(idx)     = newp(i);
    #                     tour(idx+1:end+1) = aux_tour(idx:length(tour));
    #                 else
    #                     tour(idx+1)     = newp(i);
    #                     tour(idx+2:length(tour)+1) = aux_tour(idx+1: ...
    #                         length(tour));
    #                 end
    #
    #             # Particular cases: when the nearest neighbour is the first
    #             # or last element in the tour
    #             elseif idx == 1
    #                 dist = norm(tour{2} - newp{i});
    #                 if dist < mindist
    #                     tour(2) = newp(i);
    #                     tour(3:length(tour)+1) = aux_tour(2:length(tour));
    #                 else
    #                     tour(1) = newp(i);
    #                     tour(2:length(tour)+1) = aux_tour(:);
    #                 end
    #             else
    #                 dist = norm(tour{idx - 1} - newp{i});
    #                 if dist < mindist
    #                     tour(idx)     = newp(i);
    #                     tour(idx + 1) = aux_tour(idx);
    #                 else
    #                     tour(idx + 1) = newp(i);
    #                 end
    #             end
    #         end
    #     case 'manhattan'
    #
    #         # Analyze the tour length when inserting the point at each
    #         # position in the tour
    #         for i=1:length(newp)
    #             mindelta = inf;
    #             for j=1:length(tour)
    #                 if j ~= length(tour)
    #                     %oldd = norm(tour{j+1} - tour{j});
    #                     oldd = manhattan(tour{j+1}, tour{j});
    #                     %newd = norm(tour{j+1} - newp{i}) + ...
    #                     %    norm(tour{j} - newp{i});
    #                     newd = manhattan(tour{j+1}, newp{i}) + ...
    #                         manhattan(tour{j}, newp{i});
    #                     delta = newd - oldd; % increase of tour length when
    #                     % adding the new insertion point
    #                 else
    #                     %delta = norm(tour{j} - newp{i});
    #                     delta = manhattan(tour{j}, newp{i});
    #                 end
    #                 if delta < mindelta
    #                     mindelta = delta;
    #                     idx = j;
    #                 end
    #             end
    #
    #             # Insert element in the position that minimizes the tour
    #             # deviation length
    #             if idx ~= length(tour)
    #                 aux_tour = tour;
    #                 tour(idx+1) = newp(i);
    #                 tour(idx+2:end+1) = aux_tour(idx+1:end);
    #             else
    #                 tour(end+1) = newp(i);
    #             end
    #
    #         end
    #     case 'simulated annealing'
    #
    # end
    @staticmethod
    def map2grid(map):
        """
        This function takes a map, defined as a list of lists where the outermost
        rows and columns are filled with NaNs to denote boundaries or areas
        outside of interest, and converts it into a grid by removing these NaN
        borders. Within the resulting grid, any list containing NaN values is
        emptied, signifying that it does not contain useful data or represents an
        area outside the region of interest.

        Programmers:  Paula Betriu (UPC/ESEIAAT)
        Date:         09/2022

        Usage:        grid = map2grid(map)

        Inputs:
          > map:          list of lists representing the map. It is the input grid
                          with an added border of NaN values. Inside the grid,
                          empty cells are replaced with [NaN NaN]

        Returns:
          > grid:         list of lists that contains a 2-element vector for grid
                          points within the ROI, or is empty for points outside
                          the ROI or excluded from coverage.
        """
        grid = [row[1:-1] for row in map[1:-1]]
        for i in range(len(grid)):
            for j in range(len(grid[i])):
                if np.size(grid[i][j]) == 1 and grid[i][j] == None:
                    grid[i][j] = grid[i][j]
                elif (np.isnan(grid[i][j])).all():
                    grid[i][j] = None

        return grid

    @staticmethod
    def boustrophedonMod(grid, dir1, dir2):

        # Previous check...
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
                bearing = True  # left -> right
            else:
                bearing = False  # right -> left

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
