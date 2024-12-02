import numpy as np
from conversion_functions import *
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.instpointing_gpt import instpointing


def inst2topo(grid, lon, lat, target, sc, inst, et):
    """
    Transforms a set of points from the instrument frame to the
    topographic coordinate system (latitude and longitude on the target body).

    Parameters:
        grid (list of lists): Grid points in instrument frame coordinates.
            Each element is a 2D point [x, y], or None if no observation point is defined.
        lon (float): Longitude of the observation point or area center, in degrees.
        lat (float): Latitude of the observation point or area center, in degrees.
        target (str): Name of the target body.
        sc (str): Name of the spacecraft.
        inst (str): Name of the instrument.
        et (float): Ephemeris time, TDB seconds past J2000 epoch.

    Returns:
        grid_topo (list of lists): The input grid points transformed to
            topographical coordinates on the target body. Each element is a 2D point [lon, lat], in degrees.
    """

    # Pre-allocate
    # Call instpointing to get rotation matrix rotmat
    # Note: instpointing is assumed to be implemented in Python with the same interface
    fovbounds, boresight, rotmat, visible, lon_out, lat_out = instpointing(inst, target, sc, et, lon, lat)
    # Initialize grid_topo as a list of lists, same size as grid
    grid_topo = [[None for _ in row] for row in grid]
    method = 'ELLIPSOID'
    # Get target frame name in SPICE
    _, targetframe, _ = mat2py_cnmfrm(target)
    targetframe = targetframe[0][0]

    # Convert grid into topographical coordinates
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            instp = grid[i][j]  # Retrieve current point in instrument frame
            if instp is not None and not np.any(np.isnan(instp)):
                # Initialize 3D point for conversion
                p = np.zeros(3)
                p[0:2] = instp  # Assign [x, y] coordinates
                p[2] = 1  # Set z to 1 for surface projection calculation
                # Transform point to body-fixed frame
                p_body = np.dot(rotmat, p)
                # Compute surface intersection of the point on the target body
                xpoint, _, _, found = mat2py_sincpt(
                    method, target, et, targetframe, 'NONE', sc, targetframe, p_body)
                if found:
                    # Convert rectangular coordinates to latitudinal
                    _, lon_rad, lat_rad = mat2py_reclat(xpoint)
                    # Convert radians to degrees
                    lon_deg = lon_rad * mat2py_dpr()
                    lat_deg = lat_rad * mat2py_dpr()
                    grid_topo[i][j] = [lon_deg, lat_deg]
                else:
                    print(f"Point not visible from the instrument at grid position ({i}, {j})")
                    grid_topo[i][j] = None
            else:
                grid_topo[i][j] = None  # Maintain the same structure if point is None or NaN

    return grid_topo
