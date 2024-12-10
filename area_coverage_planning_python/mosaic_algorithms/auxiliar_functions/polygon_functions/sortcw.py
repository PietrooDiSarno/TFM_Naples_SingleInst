import copy

import numpy as np
from shapely.geometry import MultiPolygon, Polygon


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
            return x,y

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
            polygon = Polygon((list(zip(x,y))))

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
                            list(zip(xuni[nanindex[i - 1] + 1:nanindex[i]], yuni[nanindex[i - 1] + 1:nanindex[i]]))))
                if ~ np.isnan(xuni[-1]):
                    polygon_list.append(Polygon(list(zip(xuni[nanindex[-1] + 1:], yuni[nanindex[-1] + 1:]))))
                polygon = MultiPolygon(polygon_list)
            else:
                polygon = Polygon((list(zip(xuni, yuni))))

        cx = polygon.centroid.x
        cy = polygon.centroid.y # polygon centroid


        angle = np.arctan2(np.array(y) - cy, np.array(x) - cx)  # obtain angle
        ind = np.argsort(angle)[::-1]  # sort angles and get the indices
        x_sorted = np.array(x)[ind] # sort the longitude values according to the angle order
        y_sorted = np.array(y)[ind] # sort the latitude values according to the angle order
        return x_sorted, y_sorted # save output

    elif len(args) == 3:
        # Algorithm extracted from
        # https://stackoverflow.com/questions/47949485/
        # sorting-a-list-of-3d-points-in-clockwise-order
        x = copy.deepcopy(args[0])
        y = copy.deepcopy(args[1])
        z = copy.deepcopy(args[2])

        innerpoint = np.array([np.mean(np.array(x)), np.mean(np.array(y)), np.mean(np.array(z))])  # find an inner point (this
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
            
        angles=np.array(angles)
        # Sort vertices by angles
        ind = np.argsort(angles)
        x_sorted = np.array(x)[ind]
        y_sorted = np.array(y)[ind]
        z_sorted = np.array(z)[ind]
        return x_sorted, y_sorted, z_sorted
    else:
        raise ValueError("Too many input arguments")