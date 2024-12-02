import numpy as np
from shapely.geometry import Polygon


def minimumWidthDirection(x, y):
    """
    This function computes the orientation (angle) at which the width of the
    polygon is minimized. It also calculates the minimum width, the height
    (assuming the minimum width as the base), and the direction vector (axes)
    of the minimum width.

    Usage:        thetamin, minwidth, height, axes = minimumWidthDirection(x, y)

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

    # Ensure x and y are numpy arrays for efficient numerical operations
    x = np.array(x)
    y = np.array(y)

    # Check if the polygon is divided in two (a.m. intersection)
    # MATLAB code:
    # ind = find(isnan(x));
    # if ~isempty(ind)
    #     x(1:ind) = x(1:ind) + 360;
    #     x(ind) = [];
    #     y(ind) = [];
    # end
    #
    # Conversion notes:
    # - Use np.isnan and np.where to find NaN indices.
    # - Assumed at most one NaN in x, as per original code logic.

    ind = np.where(np.isnan(x))[0]

    if ind.size > 0:
        ind = ind[0]  # Take the first occurrence
        # Adjust x values before the NaN index by adding 360
        x[:ind] = x[:ind] + 360
        # Remove the NaN entries from x and y
        x = np.delete(x, ind)
        y = np.delete(y, ind)

    # Find centroid of the polygon
    # MATLAB code:
    # [cx, cy] = centroid(polyshape(x, y));
    #
    # Conversion notes:
    # - Used shapely.geometry.Polygon for polygon creation and centroid calculation.

    polygon = Polygon(zip(x, y))
    centroid = polygon.centroid
    cx, cy = centroid.x, centroid.y

    # Sort the vertices in clockwise direction
    # MATLAB code:
    # [vertices(:,1), vertices(:,2)] = sortcw(x, y);
    #
    # Conversion notes:
    # - Implemented sorting by calculating angles from centroid to each point.

    # Stack x and y into an array of points
    vertices = np.column_stack((x, y))

    # Compute angles from centroid to each point
    angles = np.arctan2(y - cy, x - cx)
    # Sort indices in descending order (clockwise direction)
    idx = np.argsort(-angles)
    vertices = vertices[idx]

    # Minimum width direction
    npoints = 361
    angle = np.linspace(0, 360, npoints)  # Angles from 0 to 360 degrees

    # Compute cosine and sine of angles in degrees
    cosa = np.cos(np.radians(angle))
    sina = np.sin(np.radians(angle))

    # Initialize rotation matrices for each angle
    rotMats = np.zeros((npoints, 2, 2))
    for i in range(npoints):
        rotMats[i] = np.array([[cosa[i], sina[i]],
                               [-sina[i], cosa[i]]])

    # Initialize variables to track minimum width and corresponding angle index
    minwidth = np.inf
    mini = None

    # Loop through each rotation angle to find the minimum width
    for i in range(npoints):
        rotMat = rotMats[i]
        # Center vertices at centroid
        centered_vertices = vertices - [cx, cy]
        # Rotate the vertices
        rotVertices = (rotMat @ centered_vertices.T).T

        # Obtain width length of the rotated polygon (max x - min x)
        maxx = np.max(rotVertices[:, 0])
        minx = np.min(rotVertices[:, 0])
        l = maxx - minx  # Current width

        # Update minimum width and index if a smaller width is found
        if l < minwidth:
            minwidth = l
            mini = i

    # Retrieve the angle corresponding to the minimum width
    thetamin = angle[mini]

    # Compute height as area divided by minimum width
    height = polygon.area / minwidth

    # Constrain angle between 0º and 180º
    if thetamin >= 180:
        thetamin -= 180

    # Compute axes as unit vector in the direction of thetamin
    axes = [np.cos(np.radians(thetamin)), np.sin(np.radians(thetamin))]

    return thetamin, minwidth, height, axes
