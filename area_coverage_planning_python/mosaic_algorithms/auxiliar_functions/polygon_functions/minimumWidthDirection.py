import copy

import numpy as np
from scipy.spatial import ConvexHull
from shapely.geometry import Polygon

from mosaic_algorithms.auxiliar_functions.polygon_functions.sortcw import sortcw


def minimumWidthDirection(x_, y_):
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
    poly_aux = Polygon(list(zip(x,y)))
    cx = poly_aux.centroid.x
    cy = poly_aux.centroid.y

    vertices=np.zeros([np.size(x),2])
    # Sort the vertices in clockwise direction
    vertices[:,0],vertices[:,1] = sortcw(x,y)

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