from conversion_functions import *  # Import all conversion functions

def groundtrack(obs, t, target):
    """
    This function returns the spacecraft ground track across the target
    surface, at time t.

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2023

    Usage:        gtlon, gtlat = groundtrack(obs, t, target)

    Inputs:
      > obs:      string SPICE name of the observer body
      > t:        time epoch in TDB seconds past J2000 epoch. It can be
                  either a single point in time or a discretized vector of
                  different time values
      > target:   string SPICE name of the target body

    Outputs:
      > gtlon:    longitude coordinate of the observer ground track, in
                  [deg]
      > gtlat:    latitude coordinate of the observer ground track, in
                  [deg]
    """
    # Pre-allocate variables
    method = 'INTERCEPT/ELLIPSOID'  # target modeling
    abcorr = 'NONE'  # aberration correction
    _, tframe, _ = mat2py_cnmfrm(target)  # body-fixed frame

    # Get ground track
    sctrack, _, _ = mat2py_subpnt(method, target, t, tframe, abcorr, obs)  # sub-spacecraft point
    _, gtlon, gtlat = mat2py_reclat(sctrack)  # convert to latitudinal coordinates

    gtlon = gtlon * mat2py_dpr()  # [rad] to [deg]
    gtlat = gtlat * mat2py_dpr()  # [rad] to [deg]

    return gtlon, gtlat
