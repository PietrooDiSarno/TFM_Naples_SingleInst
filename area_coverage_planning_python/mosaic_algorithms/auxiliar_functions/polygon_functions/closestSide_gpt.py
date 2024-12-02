import numpy as np
from shapely.geometry import Polygon, LineString
# Import mat2py wrapper functions corresponding to MATLAB's SPICE functions
from conversion_functions import *


def closestSide(target, sc, t, targetArea, angle):
    """
    Given a region-of-interest, this function defines what is the spacecraft
    ground track position with respect to the edges of the target area.

    Usage: dir1, dir2 = closestSide(target, sc, t, targetArea, angle)

    Inputs:
        target:     SPICE string name of the target body.
        sc:         SPICE string name of the spacecraft.
        t:          Time in TDB seconds past J2000 epoch.
        targetArea: Numpy array containing the vertices of the ROI polygon.
                    The vertex points are expressed in 2D latitudinal coordinates.
                    targetArea[:,0] corresponds to the x values of the vertices.
                    targetArea[:,1] corresponds to the y values of the vertices.
        angle:      Footprint's angle in degrees.

    Outputs:
        dir1:       String defining the spacecraft ground track position
                    with respect to the ROI ('north', 'south', 'east', 'west').
        dir2:       String defining the spacecraft movement direction.
    """

    # Parameters
    # Equivalent to MATLAB: [~, targetFrame, ~] = cspice_cnmfrm(target);
    # This obtains the target body-fixed reference frame
    _, targetFrame, _ = mat2py_cnmfrm(target)

    # Rotate ROI according to the footprint's angle
    # Equivalent to MATLAB: angle = -angle*cspice_rpd;
    angle = -angle * mat2py_rpd()  # Convert angle to radians and invert sign

    # Create rotation matrix
    # Equivalent to MATLAB rotmat calculation
    rotmat = np.array([[np.cos(angle), -np.sin(angle)],
                       [np.sin(angle),  np.cos(angle)]])

    # Compute the centroid of the ROI polygon
    # Equivalent to MATLAB: [cx, cy] = centroid(polyshape(targetArea(:,1), targetArea(:,2)));
    polygon = Polygon(targetArea)
    centroid = polygon.centroid
    cx, cy = centroid.x, centroid.y

    # Rotate each point in the ROI about the centroid
    # Initialize rotated ROI array
    roi = np.zeros_like(targetArea)
    for j in range(len(targetArea)):
        # Shift point to origin (relative to centroid)
        point = targetArea[j, :] - np.array([cx, cy])
        # Rotate point
        rotated_point = rotmat @ point
        # Shift back to original position
        roi[j, :] = np.array([cx, cy]) + rotated_point

    # Compute spacecraft sub-observer point to see what side of the target area is closer to it
    # Assumption: Tri-axial ellipsoid to model the target surface
    # Equivalent to MATLAB: subobs = cspice_subpnt('NEAR POINT/ELLIPSOID', target, t, targetFrame, 'NONE', sc);
    subobs = mat2py_subpnt('NEAR POINT/ELLIPSOID', target, t, targetFrame, 'NONE', sc)
    # Convert rectangular coordinates to latitudinal coordinates
    # Equivalent to MATLAB: [~, sclon, sclat] = cspice_reclat(subobs);
    _, sclon, sclat = mat2py_reclat(subobs)
    # Convert radians to degrees
    sclon *= mat2py_dpr()
    sclat *= mat2py_dpr()

    # Rotate spacecraft sub-observer point coordinates
    point = np.array([sclon, sclat]) - np.array([cx, cy])
    aux = rotmat @ point + np.array([cx, cy])
    sclon, sclat = aux[0], aux[1]

    # Compute spacecraft sub-observer point 5 minutes later to determine movement direction
    subobs_ = mat2py_subpnt('NEAR POINT/ELLIPSOID', target, t + 5*60, targetFrame, 'NONE', sc)
    _, sclon_, sclat_ = mat2py_reclat(subobs_)
    sclon_ *= mat2py_dpr()
    sclat_ *= mat2py_dpr()

    # Rotate future spacecraft sub-observer point coordinates
    point_ = np.array([sclon_, sclat_]) - np.array([cx, cy])
    aux = rotmat @ point_ + np.array([cx, cy])
    sclon_, sclat_ = aux[0], aux[1]

    # Find the 4 boundary vertices of the rotated ROI
    maxlon = np.max(roi[:, 0])
    minlon = np.min(roi[:, 0])
    maxlat = np.max(roi[:, 1])
    minlat = np.min(roi[:, 1])

    # ROI's boundary box
    xlimit = [minlon, maxlon]
    ylimit = [minlat, maxlat]
    xbox = [xlimit[0], xlimit[0], xlimit[1], xlimit[1], xlimit[0]]
    ybox = [ylimit[0], ylimit[1], ylimit[1], ylimit[0], ylimit[0]]

    # Define line between the centroid and the ground track position
    x = [sclon, cx]
    y = [sclat, cy]
    line = LineString([(sclon, sclat), (cx, cy)])
    bbox = Polygon(zip(xbox, ybox))

    # Find intersection between the line and the bounding box
    intersection = line.intersection(bbox)

    if intersection.is_empty:
        # The ground track is inside the ROI's boundary box (no intersection)

        # Define the four corners of the bounding box
        p1 = np.array([maxlon, maxlat])
        p2 = np.array([maxlon, minlat])
        p3 = np.array([minlon, minlat])
        p4 = np.array([minlon, maxlat])

        # Calculate mid-points of the 4 edges
        midp = np.zeros((4, 2))
        midp[0, :] = 0.5 * (p1 + p4)  # midup (north)
        midp[1, :] = 0.5 * (p2 + p3)  # middown (south)
        midp[2, :] = 0.5 * (p3 + p4)  # midleft (west)
        midp[3, :] = 0.5 * (p1 + p2)  # midright (east)

        # Find the minimum distance between spacecraft position and mid-points
        mindist = np.inf
        minidx = -1
        for i in range(4):
            dist = np.linalg.norm(np.array([sclon, sclat]) - midp[i, :])
            if dist < mindist:
                mindist = dist
                minidx = i

        # Determine the closest side based on minidx
        if minidx == 0:
            dir1 = 'north'
        elif minidx == 1:
            dir1 = 'south'
        elif minidx == 2:
            dir1 = 'west'
        elif minidx == 3:
            dir1 = 'east'
    else:
        # There is an intersection between the line and the bounding box
        # Extract intersection point(s)
        if intersection.geom_type == 'MultiPoint':
            # Multiple intersection points, take the first one
            xi, yi = list(intersection.geoms)[0].x, list(intersection.geoms)[0].y
        elif intersection.geom_type == 'Point':
            xi, yi = intersection.x, intersection.y
        else:
            # Intersection is not a point, take centroid
            xi, yi = intersection.centroid.x, intersection.centroid.y

        # Determine dir1 based on which side the intersection occurs
        if abs(xi - maxlon) < 1e-2:
            dir1 = 'east'
        elif abs(xi - minlon) < 1e-2:
            dir1 = 'west'
        elif abs(yi - minlat) < 1e-2:
            dir1 = 'south'
        elif abs(yi - maxlat) < 1e-2:
            dir1 = 'north'

        # Determine if the spacecraft is moving towards or away from the closest side
        dist = np.linalg.norm(np.array([sclon, sclat]) - np.array([cx, cy]))
        dist_ = np.linalg.norm(np.array([sclon_, sclat_]) - np.array([cx, cy]))
        if dist_ < dist:
            # Spacecraft is moving towards the centroid, flip dir1
            if dir1 == 'north':
                dir1 = 'south'
            elif dir1 == 'south':
                dir1 = 'north'
            elif dir1 == 'east':
                dir1 = 'west'
            elif dir1 == 'west':
                dir1 = 'east'

    # Calculate how the spacecraft ground track position is moving along the map
    # The coverage path depends on the spacecraft's velocity direction
    if dir1 in ['north', 'south']:
        # Horizontal sweep
        if (sclon - sclon_) >= 0:
            # Spacecraft is moving leftwards
            dir2 = 'west'
        else:
            # Spacecraft is moving rightwards
            dir2 = 'east'
    elif dir1 in ['east', 'west']:
        if (sclat - sclat_) < 0:
            # Spacecraft is moving upwards
            dir2 = 'north'
        else:
            # Spacecraft is moving downwards
            dir2 = 'south'

    return dir1, dir2
