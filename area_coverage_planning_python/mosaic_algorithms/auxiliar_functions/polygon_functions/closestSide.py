import numpy as np
from shapely.geometry import MultiPolygon, Polygon, LineString
from math import cos, sin, radians
import math
from conversion_functions import *

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
    angle = -angle*mat2py_rpd()
    rotmat = np.array([[cos(angle), -sin(angle)],
                       [sin(angle), cos(angle)]])
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
        poly_aux = Polygon((list(zip(targetArea[:, 0], targetArea[:, 1]))))

    cx,cy = poly_aux.centroid.x, poly_aux.centroid.y

    roi = np.zeros((max(np.shape(targetArea)), 2))
    for j in range(max(np.shape(targetArea))):
        roi[j, :] = np.dot(rotmat, (targetArea[j] - np.array([cx, cy]))) + np.array([cx, cy])

    # Adjust ground track position for spacecraft movement analysis with respect to the oriented area
    # Assumption: Tri-axial ellipsoid to model the target surface
    # Initial position
    sclon, sclat = np.dot(rotmat, ((np.array(gt1)).reshape(2,) - np.array([cx, cy]))) + np.array([cx, cy])
    # Subsequent position
    sclon_, sclat_ = np.dot(rotmat, ((np.array(gt2)).reshape(2,) - np.array([cx, cy]))) + np.array([cx, cy])

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
    intersection = line.intersection(boundary_box) #check if the line intersects with the boundary box to determine
    # closest side

    if not intersection.is_empty:
        if intersection.geom_type == 'Point':
            xi, yi = intersection.x, intersection.y
        elif intersection.geom_type == 'MultiPoint':
            xi, yi = zip(*[(pt.x, pt.y) for pt in intersection])
            #if math.sqrt((xi[-1] - xi[0]) ** 2 + (yi[-1] - yi[0]) ** 2) < 0.026:
            xi = xi[0]
            yi = yi[0]

        elif intersection.geom_type == 'LineString':
            xi, yi = zip(*intersection.coords)
            #if math.sqrt((xi[1] - xi[0]) ** 2 + (yi[1] - yi[0]) ** 2) < 0.026:
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
