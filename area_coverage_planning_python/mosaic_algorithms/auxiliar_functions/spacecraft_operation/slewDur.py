import numpy as np
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.instpointing import instpointing
def slewDur(p1, p2, t, tobs, inst, target, sc, slew_rate):
    """
    This function determines the time it takes for a spacecraft or an
    instrument mounted on a spacecraft to rotate from one pointing direction
    to another. The calculation is based on the initial and final pointing
    vectors, the slew rate of the spacecraft or instrument, and the current
    observation time.

    Assumptions:
    - We assume that the slew rate is constant
    - We assume that the spacecraft or instrument can slew at this rate in
      any direction (and simultaneously in the three directions)
    - No rotation constraints are considered

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2023

    Usage:        tdur = slewDur(p1, p2, t, tobs, inst, target, sc, slew_rate)

    Inputs:
      > p1:           coordinates [longitude, latitude] of the initial
                      pointing direction, in [deg]
      > p2:           coordinates [longitude, latitude] of the final pointing
                      direction, in [deg]
      > t:            current time in ephemeris seconds past J2000 epoch
      > inst:         string name of the instrument
      > target:       string name of the target body
      > sc:           string name of the spacecraft
      > slew_rate:    slew rate of the spacecraft or instrument, in [rad/s]

    Output:
      > tdur:         duration required to complete the slew, in [sec]
    """

    # Pre-allocate variables
    maxit = 10
    epsilon = 1e-2

    # Rotation matrix corresponding to the initial pointing direction using the
    # instrument pointing information
    _, _, R1, _ = instpointing(inst, target, sc, t, p1[0], p1[1])

    # Initial slew duration
    init_slew = 10.
    for i in range(maxit):
        # Get final pointing matrix
        _, _, R2, _ = instpointing(inst, target, sc, t + tobs + init_slew, p2[0], p2[1])

        # Relative rotation matrix between the two positions
        Rdelta = np.transpose(R1) @ R2

        # Angle of rotation required
        angle = np.arccos((np.trace(Rdelta) - 1) / 2)
        if angle >= np.pi / 2:
            angle = np.pi - angle

        # Duration of the slew
        tdur = angle / slew_rate
        new_slew = tdur

        # Check if the slew duration has converged
        if abs(new_slew - init_slew) < epsilon:
            return tdur

        # Update estimate
        init_slew = new_slew

    return tdur  # Return the final computed duration
