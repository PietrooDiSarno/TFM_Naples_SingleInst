import numpy as np
from conversion_functions import *
from mosaic_algorithms.auxiliar_functions.grid_functions.topo2inst import topo2inst
from mosaic_algorithms.auxiliar_functions.plot.trgobsvec import trgobsvec

def fovray(inst, target, obs, et, lon, lat, *varargin):
    # Pre-allocate variables
    _, rframe, _ = mat2py_cnmfrm(target)  # target frame ID in SPICE

    if len(varargin) > 0:
        # Determine spacecraft attitude according to lon, lat boresight pointing
        blon = varargin[0]
        blat = varargin[1]
        if np.isscalar(lon) and np.isscalar(lat):
            coords = np.array([[lon, lat]])
        else:
            coords = np.column_stack((lon, lat))

        # Project latitudinal coordinates on body surface to focal plane
        instpoint = topo2inst(coords, blon, blat, target, obs, inst, et)
        if np.isnan(instpoint).any() or instpoint.size == 0:
            # Point is not visible
            visible = False
            return visible

        # If point is visible, let's check if it is inside the FOV limits
        # Retrieve FOV parameters
        _, _, _, bounds = mat2py_getfov(mat2py_bodn2c(inst)[0],4)  # instrument FOV's boundary vectors in the
        # instrument frame
        # Get min-max FOV boundaries in the focal plane
        maxx = np.max(bounds[0, :])
        minx = np.min(bounds[0, :])
        maxy = np.max(bounds[1, :])
        miny = np.min(bounds[1, :])

        if minx <= instpoint[0][0] <= maxx and miny <= instpoint[0][1] <= maxy:
            # Point is visible and in FOV
            visible = True
        else:
            # Point is visible but not in FOV
            visible = False
    else:
        # Check point in FOV by retrieving spacecraft attitude from CK
        abcorr = 'NONE'
        raydir = -trgobsvec([lon, lat], et, target, obs)
        visible = mat2py_fovray(inst, raydir, rframe, abcorr, obs, et)

    return visible
