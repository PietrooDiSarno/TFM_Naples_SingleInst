import numpy as np
import math
from geopy.distance import geodesic
from mosaic_algorithms.auxiliar_functions.polygon_functions.sortcw import  sortcw

def interppolygon(roi0):
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
        dist = np.linalg.norm(p2 - p1) # calculate Euclidean distance between points
        if dist == 0:
            continue  # skip if points are identical
        if dist < epsilon:
            epsilon = dist  # update minimum distance

    # Find number of regions of polygon
    indnan = np.where(np.isnan(roi0[:, 0]))[0]
    newlat = []
    newlon = []
    if len(indnan) > 0:
        #indnan[-1] =  roi0.shape[0]
        indnan = np.append(indnan,roi0.shape[0])
        from_idx = 0
        for i in range(len(indnan)):
            to = indnan[i]
            # Perform great circles interpolation of the latitude and longitude based
            # on the minimum distance found
            latd = lat[from_idx:to]
            lond = lon[from_idx:to]
            auxlat, auxlon = interpm(latd, lond, math.ceil(epsilon / 2), 'gc')
            newlat.extend(auxlat)
            newlat.append(np.nan)
            newlon.extend(auxlon)
            newlon.append(np.nan)
            from_idx = to + 1
    else:
        # Perform great circles interpolation of the latitude and longitude based
        # on the minimum distance found
        newlat, newlon = interpm(lat, lon, math.ceil(epsilon / 2), 'gc')

    # Update the roi array with the interpolated latitude and longitude values
    roi = np.column_stack((newlon, newlat))

    # Sort coordinates in clockwise order
    # [roi[:, 0], roi[:, 1]] = sortcw(roi[:, 0], roi[:, 1])
    return roi

from geopy.distance import great_circle
from geopy.point import Point

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


def interpm(lat, lon, maxdiff, method='gc'):

    maxdiff = maxdiff*(math.pi/180)*6371.009

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
        if round(dist % maxdiff) == 0 and dist % maxdiff <=0.05:
            num_points = num_points - 1

        for j in range(1, num_points + 1):
            fraction = j / (num_points + 1)
            azimuth = calculate_azimuth(start,end)
            intermediate_point = great_circle(kilometers=fraction * dist).destination((start.latitude,start.longitude), azimuth)
            latout.append(intermediate_point.latitude)
            lonout.append(intermediate_point.longitude)

        latout.append(lat[i])
        lonout.append(lon[i])

    return latout, lonout


