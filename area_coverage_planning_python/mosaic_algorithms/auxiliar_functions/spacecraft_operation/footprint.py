import copy
import numpy as np
from conversion_functions import *
from mosaic_algorithms.auxiliar_functions.observation_geometry.emissionang import emissionang
from mosaic_algorithms.auxiliar_functions.polygon_functions.amsplit import amsplit
from mosaic_algorithms.auxiliar_functions.polygon_functions.minimumWidthDirection import minimumWidthDirection
from mosaic_algorithms.auxiliar_functions.polygon_functions.sortcw import sortcw
from mosaic_algorithms.auxiliar_functions.polygon_functions.amsplit import amsplit
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.instpointing import instpointing
from science_opportunity_main.queries.geometric.fovray import fovray
import warnings
from shapely.geometry import MultiPolygon, Polygon, Point

def footprint(t, inst, sc, target, res, *args):
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
    _,targetframe,_ = mat2py_cnmfrm(target)  # target frame ID in SPICE
    abcorr = 'LT+S'  # one-way light time aberration correction parameter
    method = 'ELLIPSOID'  # assumption: ray intercept function is going to model the target
    # body as a tri-axial ellipsoid
    surfPoints = np.array([])  # matrix that saves the rectangular coordinates of the intercept points between
    # the FOV perimeter and the body surface
    count = 0  # counter of found intercept points
    geom = False
    ckpointing,lon,lat = [],[],[]

    # Define instrument pointing (3-axis steerable or constrained)
    if len(args) == 0:  # instrument pointing is provided by (and retrieve from) a ck
        ckpointing = True
    elif len(args) >= 1:  # instrument is considered 3-axis steerable
        lon, lat = args[0], args[1]
        ckpointing = False
        if len(args) == 3:
            geom = args[2]

    # Definition of footprint resolution
    if res == 'lowres': # footprint vertices resolution
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
        'bvertices':np.array([]), # footprint boundary vertices, in latitudinal coordinates, in [deg]
        'olon': np.nan, # longitude value of FOV boresight projection onto the body surface, in [deg]
        'olat': np.nan, # latitude value of FOV boresight projection onto the body surface, in [deg]
        'fovbsight': np.array([]), # FOV boresight in target frame centered at the spacecraft position
        'fovbounds': np.array([]), # FOV bounds in target frame centered at the spacecraft position
        'limb': 'none', # boolean that defines if the FOV projects onto the planetary body's limb
        'recVertices': np.array([]), # matrix that contains the rectangular coordinates of the footprint boundary vertices
        'angle': np.nan,
        'width': np.nan,
        'height': np.nan
    }

    # Calculate instrument orientation
    if ckpointing: # constrained pointing
        bounds, boresight, pointingRotation, found, lon, lat = instpointing(inst, target, sc, t)
    else:
        bounds, boresight, pointingRotation, found = instpointing(inst, target, sc, t, lon, lat)

    if not found:
        return fp  # the point is not visible from the instrument's FOV, therefore
                   # the function is exited and the footprint is returned empty

    fp['fovbounds'] = bounds # save fov bounds in the body-fixed reference frame
    fp['fovbsight'] = boresight # save fov boresight vector in the body-fixed reference frame

    # Update boresight pointing (in latitudinal coordinates)
    fp['olon'] = lon # longitude value of FOV boresight projection onto the body surface, in [deg]
    fp['olat'] = lat # latitude value of FOV boresight projection onto the body surface, in [deg]

    # Project instrument FOV onto body surface

    # Retrieve FOV parameters
    _, _, _, bounds = mat2py_getfov((mat2py_bodn2c(inst))[0], 4) # instrument FOV's boundary
    # vectors in the instrument frame
    minx = np.min(bounds[0, :]) # minimum x focal plane
    maxx = np.max(bounds[0, :]) # maximum x focal plane
    miny = np.min(bounds[1, :]) # minimum y focal plane
    maxy = np.max(bounds[1, :]) # maximum y focal plane
    z = bounds[2, 0] # z-coordinate of the boundary vectors

    boundPoints = np.zeros((3, max(bounds.shape))) # intercept points of the FOV
    # boundary, in Cartesian coordinates, and in the body-fixed reference frame
    # In its simplest form, the footprint should have the same number of
    # vertices as boundaries has the instrument's FOV

    fp['limb'] = 'none' # boolean to define if the instrument FOV projection is not
    # enclosed in the body surface, i.e. at least one of the boundary vectors
    # do not intercept the body surface

    intsec = False # boolean that indicates if at least one of the boundary
    # vectors intercepts the body surface

    for i in range(max(fp['fovbounds'].shape)):
        boundPoints[:, i], _, _, found = mat2py_sincpt(method, target, t, targetframe, abcorr, sc, targetframe, fp['fovbounds'][:, i])
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

    def refineFOVsearch(fp,maxx,minx,maxy,miny,pointingRotation, method, target, t, targetframe, abcorr,sc):
        # Perform a perimetral search of the FOV to find out if the FOV
        # contains totally or partially the body limb
        Nl = 20
        found = False
        for ii in range(2):
            # Vertical sweep of the focal perimeter
            x = (maxx - minx) * ii + minx
            for jj in range(Nl + 1):
                y = (maxy - miny) * jj/Nl + miny
                vec = np.array([x, y, z]).reshape([3,1])
                vec = np.dot(pointingRotation,vec)  # transform vector coordinates to target frame
                _, _, _, found = mat2py_sincpt(method, target, t, targetframe, abcorr, sc, targetframe,
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
                    vec = np.array([x, y, z]).reshape([3,1])
                    vec = np.dot(pointingRotation,vec)  # transform vector coordinates to target frame
                    _, _, _, found =mat2py_sincpt(method, target, t, targetframe, abcorr, sc,
                                                   targetframe, vec)
                    if found:
                        fp['limb'] = 'partial'
                        break
                if found:
                    break
        return  fp,maxx,minx,maxy,miny,pointingRotation, method, target, t, targetframe, abcorr,sc

    # When the footprint is likely to contain the limb, perform a more
    # accurate search in order to conclude if the FOV intercept the
    # body at some point
    if fp['limb'] == 'total':
        fp,maxx,minx,maxy,miny,pointingRotation, method, target, t, targetframe, abcorr,sc = refineFOVsearch(fp,maxx,
                                            minx,maxy,miny,pointingRotation, method, target, t, targetframe, abcorr,sc)

    def inFOVprojection(boundPoints,N,surfPoints,count):
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

                _,alt,_, aux,_,_ = mat2py_tangpt(method, target, t, targetframe, abcorr, corloc, sc, targetframe,
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

    def tlimbFOVprojection(surfPoints,target,t,targetframe,sc):
        # Compute limb with SPICE function (easier)
        # Parameters for spice.limbpt function
        lbmethod = 'TANGENT/ELLIPSOID'
        abcorr = 'XLT+S'
        corloc = 'CENTER'
        refvec = np.array([0,0,1]).reshape(3,1) # first of the sequence of cutting half-planes
        ncuts=int(2e3) # number of cutting half-planes
        delrol = mat2py_twopi() / ncuts # angular step by which to roll the
        # cutting half-planes about the observer-target vector
        schstp = 1.0e-6 # search angular step size
        soltol = 1.0e-10 # solution convergence tolerance

        # Limb calculation with spice.limbpt function
        _,limb,_,_ = mat2py_limbpt(lbmethod,target,t,targetframe,abcorr,corloc,sc,refvec,delrol,ncuts,schstp,soltol,ncuts)
        # limb points expressed in targetframe ref frame
        surfPoints = limb.T
        return surfPoints,target,t,targetframe,sc

    # Compute footprint
    if fp['limb'] == 'none':
        # FOV projects entirely on the body surface
       boundPoints, N, surfPoints, count = inFOVprojection(boundPoints, N, surfPoints, count)
    elif fp['limb'] == 'partial':
        #FOV contains partially the limb
        surfPoints,pointingRotation, minx, maxx, miny, maxy, z, N, t, method, target, targetframe, abcorr, sc,res,count = plimbFOVprojection(surfPoints,pointingRotation, minx, maxx, miny, maxy, z, N, t, method, target, targetframe, abcorr, sc,res,count)
    else:
        #FOV contains the whole body (total limb)
        surfPoints,target,t,targetframe,sc = tlimbFOVprojection(surfPoints,target,t,targetframe,sc)

    if surfPoints.size == 0: # the FOV does not intercept with the object at any point
        # of its focal plane
        return fp

    # Save values
    fp['recVertices'] = surfPoints

    def footprint2map(surfPoints,t,target,sc,fp,inst):

        # Pre-allocate variables
        vertices = np.zeros([max(np.shape(surfPoints)),2]) # matrix that saves the
        # latitudinal coordinates of the intercept points between the FOV
        # perimeter and the body surface

        # Sort points
        surfPoints[:,0],surfPoints[:,1],surfPoints[:,2] = sortcw(surfPoints[:,0],surfPoints[:,1],surfPoints[:,2])
        # sort polygon boundary vertices in clockwise order (for representation)

        for i in range(max(np.shape(surfPoints))):
            _, auxlon, auxlat = mat2py_reclat(surfPoints[i, :].T)  # rectangular to
            # latitudinal coordinates
            vertices[i, 0] = auxlon * mat2py_dpr()  # longitude in [deg]
            vertices[i, 1] = auxlat * mat2py_dpr()  # latitude in [deg]
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
            angle = emissionang(srfpoint, t, target, sc)
            if angle < 90:
                northpole = True
            # Check south-pole
            srfpoint = np.array([0, -90])
            angle = emissionang(srfpoint, t, target, sc)
            if angle < 90:
                southpole = True

            # Case 1.
            if not northpole and not southpole:
                # Check a.m. split
                ind2 = np.where(np.diff(np.sort(lblon)) >= 180)[0]  # find the discontinuity
                if ind2.size > 0:
                    lblon, lblat = amsplit(lblon, lblat)
                # Check if we are keeping the correct polygon (full disk polygons may be
                # misleading, we can only guarantee through emission angle check)
                exit = False

                while not exit:
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
                        angle = emissionang(randPoint, t, target, sc)
                        if angle < 85:
                            exit = True
                    else:
                        angle =  emissionang(randPoint, t, target, sc)
                        if angle < 85:
                            exit = True
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
            northpole = fovray(inst, target, sc, t, 0, 90, fp['olon'], fp['olat'])
            # Check south-pole visibility:
            southpole = fovray(inst, target, sc, t, 0, -90, fp['olon'], fp['olat'])
            # Case 1.
            if not northpole and not southpole:
                lblon, lblat = amsplit(lblon, lblat)
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

            col1, col2 = amsplit(vertices[:, 0], vertices[:, 1])  # save footprint vertices
            fp['bvertices'] = np.hstack((col1.reshape(len(col1), 1), col2.reshape(len(col2), 1)))

            if geom:
                # Get minimum width direction and size
                angle, width, height = minimumWidthDirection(fp['bvertices'][:, 0], fp['bvertices'][:, 1])
                fp['angle'] = angle
                fp['width'] = width
                fp['height'] = height

        return surfPoints,t,target,sc,fp,inst

    # Conversion from rectangular to latitudinal coordinates of the polygon vertices
    # and geometry computation
    surfPoints,t,target,sc,fp,inst = footprint2map(surfPoints,t,target,sc,fp,inst)

    return fp
