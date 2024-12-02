# Import external modules
import numpy as np

# Import the mat2py wrapper functions for SPICE routines
from conversion_functions import *

# Import local modules
from mosaic_algorithms.auxiliar_functions.plot.trgobsvec_gpt import trgobsvec


def emissionang(srfpoint, t, target, obs):
    """
    Computes the phase angle between the target's surface normal and the
    distance vector to the observer from the target surface point at given time(s).

    Usage:
        angle = emissionang(srfpoint, t, target, obs)

    Inputs:
        srfpoint: Target surface point.
                  - If input as latitudinal coordinates: list or array with 2 elements [lon, lat] in degrees.
                  - If input as Cartesian coordinates: list or array with 3 elements [x, y, z] in kilometers.
        t:        Time epoch in TDB seconds past J2000 epoch.
                  - Can be a single float value or a NumPy array of time values.
        target:   String SPICE name of the target body.
        obs:      String SPICE name of the observer body.

    Output:
        angle:    Angle between the normal at the surface point and the
                  vector to the observer, in degrees.
                  - If `t` is a single value, `angle` is a float.
                  - If `t` is an array, `angle` is a NumPy array of angles.
    """

    # Parameters
    method = 'ELLIPSOID'  # Assumption: target body modeled as a tri-axial ellipsoid
    _, targetframe, _ = mat2py_cnmfrm(target)  # Retrieve the target frame ID in SPICE
    targetframe = targetframe[0][0]

    # Convert srfpoint to rectangular coordinates if provided in latitudinal coordinates
    if len(srfpoint) == 2:
        # Convert degrees to radians
        rpd = mat2py_rpd()  # Radians per degree
        lon_rad = srfpoint[0] * rpd
        lat_rad = srfpoint[1] * rpd
        # Convert to rectangular coordinates
        target_code = mat2py_bodn2c(target)
        srfpoint_rect = mat2py_srfrec(target_code[0], lon_rad, lat_rad)
    elif len(srfpoint) == 3:
        # Ensure srfpoint is a column vector
        srfpoint_rect = np.array(srfpoint).reshape(3, 1) if len(srfpoint.shape) == 1 else srfpoint.T
    else:
        raise ValueError("srfpoint must be either 2 elements (lon, lat) or 3 elements (x, y, z).")

    # Compute the observer position as seen from the surface point
    # Note: trgobsvec is a user-defined function and needs to be translated separately.
    # Here, we assume that trgobsvec is already implemented and imported.
    # If not, you need to provide the MATLAB code for trgobsvec for translation.
    obsvec, _ = trgobsvec(srfpoint_rect, t, target, obs)  # [3, N] array

    # Obtain the outward surface normal vector at the surface point
    # Initialize normal vectors array
    if np.isscalar(t):
        nrmvec = mat2py_srfnrm(method, target, t, targetframe, srfpoint_rect).reshape(3, 1)
    else:
        nrmvec = np.zeros((3, len(t)))
        for i in range(len(t)):
            nrmvec[:, i] = mat2py_srfnrm(method, target, t[i], targetframe, srfpoint_rect).flatten()

    # Compute the angle between the observer vector and the normal vector
    angle_rad = mat2py_vsep(obsvec, nrmvec)  # Angular separation in radians
    angle_deg = angle_rad * mat2py_dpr()  # Convert radians to degrees

    return angle_deg
