import numpy as np


def amsplit(lon, lat):
    """
    Splits a polygon crossing the anti-meridian (±180 degrees longitude) into
    a format suitable for polygon operations.

    Inputs:
        lon: array of longitudes (degrees)
        lat: array of latitudes (degrees)

    Outputs:
        lon_out: array of adjusted longitudes
        lat_out: array of latitudes (same as input)

    This function is equivalent to the MATLAB 'amsplit' function used in the original code.
    It adjusts the longitude values to ensure that the polygon does not have discontinuities
    when crossing the anti-meridian.
    """
    # Copy the input arrays
    lon = np.array(lon, copy=True)
    lat = np.array(lat, copy=True)

    # Wrap longitudes into the range [-180, 180]
    lon = ((lon + 180) % 360) - 180

    # Find where the difference between consecutive longitudes is greater than 180 degrees
    lon_diff = np.abs(np.diff(lon))
    idx = np.where(lon_diff > 180)[0]

    if idx.size > 0:
        # Adjust the longitudes to avoid crossing the anti-meridian
        for i in idx:
            if lon[i + 1] > lon[i]:
                lon[i + 1:] -= 360
            else:
                lon[i + 1:] += 360

    return lon, lat
