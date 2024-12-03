# Import external modules
import numpy as np

# Import the mat2py wrapper functions for SPICE routines
from conversion_functions import *


def trgobsvec(srfpoint, t, target, obs, frame=None):
    """
    Distance vector between the target point P and the observer at time t.

    Usage:
        obsvec, dist = trgobsvec(srfpoint, t, target, obs)
        obsvec, dist = trgobsvec(srfpoint, t, target, obs, frame)

    Inputs:
        srfpoint: Target surface point. It can be input either in latitudinal
                  coordinates (list or array of [lon, lat] in degrees) or
                  Cartesian coordinates (list or array of [x, y, z] in km)
                  with respect to the body-fixed reference frame.
        t:        Time epoch in TDB seconds past J2000 epoch. It can be
                  either a single scalar value or a numpy array of different
                  time values.
        target:   String SPICE name of the target body.
        obs:      String SPICE name of the observer body.
        frame:    (Optional) String SPICE name of the reference frame with respect
                  to which the vector is going to be expressed. If not provided,
                  the body-fixed reference frame is used by default.

    Outputs:
        obsvec:   Observer position vector as seen from the target surface
                  point in the target body-fixed reference frame, in [km].
                  If multiple time points are provided, obsvec will be a
                  2D array with shape (3, N), where N is the number of time points.
        dist:     Distance between the observer and the target surface point, in [km].
                  If multiple time points are provided, dist will be a 1D array of length N.
    """

    # Convert inputs to numpy arrays for consistency
    srfpoint = np.array(srfpoint, dtype=float)
    t = np.array(t, dtype=float)

    # Retrieve the target frame ID in SPICE
    _, frame_id, _ = mat2py_cnmfrm(target)
    frame_id = frame_id[0][0]
    abcorr = 'NONE'  # Geometric position, no light aberration

    # Convert srfpoint from latitudinal to rectangular coordinates if necessary
    if srfpoint.size == 2:
        # srfpoint provided as [lon, lat] in degrees
        rpd = mat2py_rpd()  # Radians per degree
        lon_rad = srfpoint[0] * rpd
        lat_rad = srfpoint[1] * rpd
        target_code = mat2py_bodn2c(target)
        srfpoint_rect = mat2py_srfrec(target_code[0], lon_rad, lat_rad)
    elif srfpoint.size == 3:
        # srfpoint provided as [x, y, z] in km
        if srfpoint.shape[0] != 3:
            srfpoint = srfpoint.flatten()
            if srfpoint.size != 3:
                raise ValueError("srfpoint must be a 3-element vector if provided in Cartesian coordinates.")
        srfpoint_rect = srfpoint
    else:
        raise ValueError("srfpoint must be either a 2-element [lon, lat] or a 3-element [x, y, z] vector.")

    # Compute the observer position as seen from the target body in the body-fixed frame
    # mat2py_spkpos returns the position vector and the light time; we ignore light time as abcorr='NONE'
    obspos, _ = mat2py_spkpos(obs, t, frame_id, abcorr, target)
    # obspos shape: (3, N) if t has multiple elements, else (3,)

    # Compute the observer vector relative to the surface point
    # Ensure srfpoint_rect has shape (3, 1) if t has multiple elements
    if t.size > 1:
        srfpoint_matrix = np.tile(srfpoint_rect.reshape(3, 1), (1, t.size))
        obsvec = obspos - srfpoint_matrix
    else:
        obsvec = obspos - srfpoint_rect.reshape(3, 1)

    # If a different reference frame is requested, perform the transformation
    if frame is not None:
        from_frame = frame_id
        to_frame = mat2py_bodn2c(target) if isinstance(frame, str) else frame  # Assuming frame is a string SPICE name
        to_frame = str(to_frame[0])
        rotmat = mat2py_pxform(from_frame, to_frame, t)

        if t.size > 1:
            # rotmat is expected to have shape (3, 3, N)
            transformed_obsvec = np.zeros_like(obsvec)
            for i in range(t.size):
                transformed_obsvec[:, i] = rotmat[:, :, i] @ obsvec[:, i]
            obsvec = transformed_obsvec
        else:
            # rotmat has shape (3, 3)
            obsvec = rotmat @ obsvec

    # Compute the distance as the Euclidean norm of the observer vector
    dist = np.linalg.norm(obsvec, axis=0)  # Shape: (N,) if t has multiple elements, else scalar

    return obsvec, dist
