import spiceypy as spice
import numpy as np
from PSOA import trgobsvec
#CAMBIAR ANTES DE DORMIR

def trgobsdist(srfpoint, t, target, obs, frame=None):
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
             dist: distance between the observer and the surface point
    """
    _, dist = trgobsvec(srfpoint, t, target, obs, frame = None)
    return dist
