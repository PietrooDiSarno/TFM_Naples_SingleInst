import numpy as np

from conversion_functions.mat2py_dpr import mat2py_dpr
from conversion_functions.mat2py_reclat import mat2py_reclat
from conversion_functions.mat2py_sincpt import mat2py_sincpt
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.instpointing import instpointing
from conversion_functions.mat2py_cnmfrm import mat2py_cnmfrm

def inst2topo(grid, lon, lat, target, sc, inst, et):
    """
    This function transforms a set of points from the instrument frame to the
    topographic coordinate system (latitude and longitude on the target body)

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022

    Usage:        grid_topo = inst2topo(grid, lon, lat, target, sc, inst, et)

    Inputs:
        grid:          A list of grid points in instrument frame coordinates.
                       Each element contains a 2D point [x,y] representing a
                       location in the instrument frame, or is empty if no
                       observation point is defined
        lon:           Longitude of the observation point or area center, in [deg].
        lat:           Latitude of the observation point or area center, in [deg].
        target:        String name of the target body.
        sc:            String name of the spacecraft.
        inst:          String name of the instrument.
        et:            Ephemeris time, TDB seconds past J2000 epoch.

    Returns:
        grid_topo:     A list of the input grid points transformed to topographic
                       coordinates on the target body. Each element contains a 2D point
                       [lon,lat], in [deg]
    """

    # Pre-allocate
    _, _, rotmat, _ = instpointing(inst, target, sc, et, lon, lat)  # Assuming instpointing function is defined
    grid_topo = [[[] for _ in range(np.shape(grid)[1])] for _ in range(np.shape(grid)[0])]  # Pre-allocate grid_topo array
    method = 'ELLIPSOID'
    _, targetframe, _ = mat2py_cnmfrm(target)  # Get target frame ID in SPICE

    # Convert grid into topographical coordinates
    for i in range(len(grid)):
        for j in range(len(grid[i])):
            instp = grid[i][j]  # retrieve current point in instrument frame
            if instp is not None and not np.all(np.isnan(instp)):  # check for empty or NaN
                p = np.zeros(3)  # initialize 3D points for conversion
                p[0:2] = instp  # assign [x, y] coordinates
                p[2] = 1.  # set z to 1 for surface projection calculation
                p_body = np.dot(rotmat, p)  # apply rotation matrix

                # Compute surface intersection of the point on the target body
                xpoint, _, _, found = mat2py_sincpt(method,target,et,targetframe,
                                            'NONE',sc,targetframe,p_body)

                if found:
                    # Convert rectangular coordinates to latitudinal
                    _, lon, lat = mat2py_reclat(xpoint)
                    grid_topo[i][j] = [lon*mat2py_dpr(), lat*mat2py_dpr()]
                else:
                    print("Point not visible from the instrument")
                    grid_topo[i][j] = None

    return grid_topo