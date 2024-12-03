# Import necessary modules
import numpy as np
from shapely.geometry import Polygon, LineString, Point
from shapely.geometry.collection import GeometryCollection
from shapely.geometry.multilinestring import MultiLineString
from shapely.geometry.multipoint import MultiPoint

from conversion_functions.mat2py_rpd import mat2py_rpd  # Wrapper for radians per degree conversion


def closestSide2(gt1, gt2, targetArea, angle):
    """
    Given a region-of-interest (ROI), this function defines what is the spacecraft
    ground track position with respect to the edges of the target area.

    Usage:
        dir1, dir2 = closestSide2(gt1, gt2, targetArea, angle)

    Inputs:
        gt1 (array-like): Spacecraft ground track position at initial time [longitude, latitude].
        gt2 (array-like): Spacecraft ground track position at later time [longitude, latitude].
        targetArea (array-like): Nx2 array containing the vertices of the ROI polygon.
                                 Each row corresponds to a vertex [x, y].
        angle (float): Footprint's angle in degrees.

    Outputs:
        dir1 (str): Direction ('north', 'south', 'east', or 'west') defining the spacecraft
                    ground track position with respect to the ROI.
        dir2 (str): Direction ('north', 'south', 'east', or 'west') indicating the movement
                    direction of the spacecraft, which determines the starting point of the coverage path.
    """

    # Conversion notes:
    # - MATLAB indices start at 1; Python indices start at 0.
    # - MATLAB uses column vectors; NumPy uses row vectors by default.
    # - The function 'polyxpoly' in MATLAB computes intersection points between polylines.
    #   In Python, we use Shapely's 'intersection' method.

    # Rotate ROI according to the footprint's angle
    # MATLAB equivalent: angle = -angle * cspice_rpd;
    angle_rad = -angle * mat2py_rpd()  # Convert angle to radians and negate

    # Build rotation matrix
    # MATLAB equivalent:
    # rotmat = [cos(angle)   -sin(angle);
    #           sin(angle)    cos(angle)];
    rotmat = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],
                       [np.sin(angle_rad),  np.cos(angle_rad)]])

    # Compute centroid of the targetArea polygon
    # MATLAB equivalent: [cx, cy] = centroid(polyshape(targetArea(:,1), targetArea(:,2)));
    target_polygon = Polygon(targetArea)
    centroid = target_polygon.centroid
    cx, cy = centroid.x, centroid.y  # Coordinates of the centroid

    # Rotate the ROI around its centroid
    # MATLAB equivalent:
    # roi  = zeros(length(targetArea), 2);
    # for j=1:length(targetArea)
    #     roi(j, :) = [cx, cy]' + rotmat*(targetArea(j, :)' - [cx, cy]');
    # end
    roi = np.zeros_like(targetArea)
    for j in range(len(targetArea)):
        roi[j, :] = np.array([cx, cy]) + rotmat @ (targetArea[j, :] - np.array([cx, cy]))

    # Compute spacecraft sub-observer point to see which side of the target area is closer
    # Assumption: Tri-axial ellipsoid to model the target surface
    # Rotate gt1 (spacecraft ground track position at initial time)
    sclon, sclat = gt1[0], gt1[1]
    aux = np.array([cx, cy]) + rotmat @ (np.array([sclon, sclat]) - np.array([cx, cy]))
    sclon, sclat = aux[0], aux[1]

    # Rotate gt2 (spacecraft ground track position at later time)
    sclon_, sclat_ = gt2[0], gt2[1]
    aux = np.array([cx, cy]) + rotmat @ (np.array([sclon_, sclat_]) - np.array([cx, cy]))
    sclon_, sclat_ = aux[0], aux[1]

    # Find the boundary vertices
    # MATLAB equivalent:
    # maxlon = max(roi(:, 1)); minlon = min(roi(:, 1));
    # maxlat = max(roi(:, 2)); minlat = min(roi(:, 2));
    maxlon = np.max(roi[:, 0])
    minlon = np.min(roi[:, 0])
    maxlat = np.max(roi[:, 1])
    minlat = np.min(roi[:, 1])

    # Define ROI's bounding box
    # MATLAB equivalent:
    # xlimit = [minlon maxlon];
    # ylimit = [minlat  maxlat];
    # xbox = xlimit([1 1 2 2 1]);
    # ybox = ylimit([1 2 2 1 1]);
    xlimit = [minlon, maxlon]
    ylimit = [minlat, maxlat]
    xbox = [xlimit[0], xlimit[0], xlimit[1], xlimit[1], xlimit[0]]
    ybox = [ylimit[0], ylimit[1], ylimit[1], ylimit[0], ylimit[0]]

    # Define line between the spacecraft position and the centroid
    # MATLAB equivalent:
    # x = [sclon cx];
    # y = [sclat cy];
    # [xi, yi] = polyxpoly(x, y, xbox, ybox);
    line = LineString([(sclon, sclat), (cx, cy)])
    bbox_polygon = Polygon(zip(xbox, ybox))
    intersection = line.intersection(bbox_polygon)

    # Determine if the ground track is inside the ROI's bounding box
    if intersection.is_empty:
        # The ground track is inside the ROI's boundary box (no intersection)
        # Calculate distance to the mid-points of the 4 edges of the boundary box

        # Define the four corners of the bounding box
        p1 = np.array([maxlon, maxlat])
        p2 = np.array([maxlon, minlat])
        p3 = np.array([minlon, minlat])
        p4 = np.array([minlon, maxlat])

        # Calculate mid-points of the edges
        midp = np.zeros((4, 2))
        midp[0, :] = 0.5 * (p1 + p4)  # midup
        midp[1, :] = 0.5 * (p2 + p3)  # middown
        midp[2, :] = 0.5 * (p3 + p4)  # midleft
        midp[3, :] = 0.5 * (p1 + p2)  # midright

        # Find the minimum distance to the mid-points
        mindist = np.inf
        minidx = -1
        for i in range(4):
            dist = np.linalg.norm(np.array([sclon, sclat]) - midp[i, :])
            if dist < mindist:
                mindist = dist
                minidx = i

        # Determine the closest side based on minimum index
        if minidx == 0:
            dir1 = 'north'
        elif minidx == 1:
            dir1 = 'south'
        elif minidx == 2:
            dir1 = 'west'
        elif minidx == 3:
            dir1 = 'east'
        else:
            dir1 = 'unknown'  # Default case
    else:
        # The ground track intersects with the ROI's boundary box
        # Extract intersection point(s)
        # Handle different intersection types
        if isinstance(intersection, Point):
            xi, yi = intersection.x, intersection.y
        elif isinstance(intersection, LineString):
            # Take the midpoint of the line as the intersection point
            xi, yi = intersection.interpolate(0.5, normalized=True).x, intersection.interpolate(0.5, normalized=True).y
        elif isinstance(intersection, (MultiPoint, MultiLineString, GeometryCollection)):
            # Take the first point
            xi, yi = list(intersection.geoms)[0].x, list(intersection.geoms)[0].y
        else:
            xi, yi = centroid.x, centroid.y  # Default to centroid if type is unexpected

        # Determine which side the intersection point is on
        if abs(xi - maxlon) < 1e-2:
            dir1 = 'east'
        elif abs(xi - minlon) < 1e-2:
            dir1 = 'west'
        elif abs(yi - minlat) < 1e-2:
            dir1 = 'south'
        elif abs(yi - maxlat) < 1e-2:
            dir1 = 'north'
        else:
            dir1 = 'unknown'  # Default case

        # Determine if the spacecraft is moving towards or away from the closest side
        dist = np.linalg.norm(np.array([sclon, sclat]) - np.array([cx, cy]))
        dist_ = np.linalg.norm(np.array([sclon_, sclat_]) - np.array([cx, cy]))
        if dist_ < dist:
            # Spacecraft is moving towards the centroid; reverse direction
            if dir1 == 'north':
                dir1 = 'south'
            elif dir1 == 'south':
                dir1 = 'north'
            elif dir1 == 'east':
                dir1 = 'west'
            elif dir1 == 'west':
                dir1 = 'east'
            else:
                dir1 = 'unknown'  # Default case

    # Calculate the spacecraft ground track movement direction
    # This determines the starting point of the coverage path
    if dir1 in ['north', 'south']:  # Horizontal sweep
        if (sclon - sclon_) >= 0:  # Spacecraft is moving leftwards
            dir2 = 'west'  # Start from the rightmost position
        else:  # Spacecraft is moving rightwards
            dir2 = 'east'  # Start from the leftmost position
    elif dir1 in ['east', 'west']:  # Vertical sweep
        if (sclat - sclat_) < 0:  # Spacecraft is moving upwards
            dir2 = 'north'  # Start from the bottom position
        else:  # Spacecraft is moving downwards
            dir2 = 'south'  # Start from the top position
    else:
        dir2 = 'unknown'  # Default case

    return dir1, dir2
