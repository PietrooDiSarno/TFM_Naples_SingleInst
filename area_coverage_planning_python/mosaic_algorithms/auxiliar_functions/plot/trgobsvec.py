from conversion_functions import *
import numpy as np


def trgobsvec(srfpoint, t, target, obs, inputframe=None):
    """
    Distance vector between the target point P and the observer at time t

    :param srfpoint: target surface point. It can be input either in latitudinal coordinates (in [deg])
                     or Cartesian coordinates (in [km]) with respect to the body-fixed reference frame
    :param t: time epoch in TDB seconds past J2000 epoch. It can be either a single point in time or
              a discretized vector of different time values
    :param target: string SPICE name of the target body
    :param obs: string SPICE name of the observer body
    :param frame: string SPICE name of the reference frame with respect to which the vector is going to be
                  expressed. If this variable is not input, the body-fixed reference frame is used by default
    :return: obsvec: observer position vector as seen from the target surface point in the target body-fixed
                     reference frame, in [km]
             dist: distance between the observer and the surface point
    """

    # Target frame
    _, frame, _ = mat2py_cnmfrm(target)

    abcorr = 'NONE'  # Assumption: geometric positions, no light aberrations

    # Convert latitudinal coordinates to rectangular if needed
    if len(srfpoint) == 2:
        srfpoint = np.radians(srfpoint)  # [deg] to [rad]
        srfpoint = mat2py_srfrec(mat2py_bodn2c(target)[0], srfpoint[0], srfpoint[1])
    else:
            srfpoint = np.array(srfpoint).reshape(3,)

    # Compute the observer position as seen from the srfpoint
    obspos, _ = mat2py_spkpos(obs, t, frame, abcorr, target)
    if np.size(t) == 1:
        obsvec = obspos - srfpoint
    else:
        obsvec = obspos - srfpoint.reshape(3,1)  # srfpoint-observer distance vector

    # If a different reference frame is requested
    if inputframe:
        # Process each time point if t is an array of times
        if np.size(t) > 1:
            # Initialize a rotation matrix list
            rotmat_list = [mat2py_pxform(frame, inputframe, time_point) for time_point in t] # rotation matrix from body-fixed reference frame to the requested one
            # Apply rotation for each time point
            for i, rotmat in enumerate(rotmat_list):
                # Ensure obsvec[:, i] is a 3D vector
                obsvec[:,i] = np.dot(rotmat, obsvec[:,i])
        else:
            # Process a single time point
            rotmat = mat2py_pxform(frame, inputframe, t)
            obsvec = np.dot(rotmat, obsvec)

    # Compute distance
    if np.size(t) >1:
        dist = np.zeros(obsvec.shape[1])
        for i in range(obsvec.shape[1]):
            dist[i] = np.linalg.norm(obsvec[:,i])
    else:
        dist = np.linalg.norm(obsvec)

    return obsvec, dist
