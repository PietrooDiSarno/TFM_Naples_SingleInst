# Import external modules
import numpy as np

# Import the mat2py wrapper functions for SPICE routines
from conversion_functions import *


def instpointing(inst, target, sc, t, *args):
    """
    This function sets the instrument's orientation, provided a target and
    the latitudinal coordinates of the point the instrument should be aiming
    at. It also checks if this point is actually visible from the FOV (could
    be on the dark side of the object as seen from the instrument).

    Usage:
        fovbounds, boresight, rotmat, visible = instpointing(inst, target, sc, t, lon, lat)
        fovbounds, boresight, rotmat, visible, lon, lat = instpointing(inst, target, sc, t)

    Inputs:
        inst:       string name of the instrument
        target:     string name of the target body
        sc:         string name of the spacecraft
        t:          observation time, i.e., the minimum time that the
                    instrument needs to perform an observation, in seconds
        lon:        (optional) longitude coordinate of the target body at which the
                    instrument boresight is pointing, in degrees
        lat:        (optional) latitude coordinate of the target body at which the
                    instrument boresight is pointing, in degrees

    Outputs:
        fovbounds:  FOV bounds in the body-fixed reference frame centered at
                    the spacecraft position at time t
        boresight:  FOV boresight in the body-fixed reference frame
                    centered at the spacecraft position at time t
        rotmat:     rotation matrix from instrument frame to target frame
        visible:    boolean that determines if the point is visible from the
                    instrument's FOV
        lon:        (if not provided) longitude coordinate of the target body at which the
                    instrument boresight is pointing, in degrees
        lat:        (if not provided) latitude coordinate of the target body at which the
                    instrument boresight is pointing, in degrees
    """

    # Pre-allocate variables
    axis3 = False  # Indicates if the spacecraft is 3-axis steerable
    if len(args) >= 2:
        lon = args[0]
        lat = args[1]
        axis3 = True
    elif len(args) == 0:
        lon = None
        lat = None
    else:
        raise ValueError("Both longitude and latitude must be provided if either is provided.")

    method = 'ELLIPSOID'  # Assumption: the target body is modeled as a tri-axial ellipsoid
    _, targetframe, _ = mat2py_cnmfrm(target)  # Retrieve the target frame ID in SPICE
    targetframe = targetframe[0][0]
    abcorr = 'LT'  # One-way light time aberration correction parameter

    # Initialize output variables for optional outputs
    lon_out = None
    lat_out = None

    # Retrieve FOV parameters
    # Equivalent to MATLAB: [shape, instframe, boresight, bounds] = cspice_getfov(cspice_bodn2c(inst), 4)
    inst_code = mat2py_bodn2c(inst)
    shape, instframe, boresight, bounds = mat2py_getfov(inst_code[0], 4)

    if shape == 'CIRCLE' or shape == 'ELLIPSE':
        raise NotImplementedError("Circular and ellipsoidal FOV shapes have not been implemented yet")

    # Pre-allocate arrays for FOV bounds and rotation matrix
    fovbounds = np.zeros((3, bounds.shape[1]))
    rotmat = np.zeros((3, 3))
    visible = False  # Initialize visibility flag

    # Pointing matrix calculation
    if axis3:  # 3-axis steerable spacecraft
        # Convert longitude and latitude to radians
        rpd = mat2py_rpd()  # Radians per degree
        lon_rad = lon * rpd
        lat_rad = lat * rpd

        # Compute the rectangular coordinates of the target point in the body-fixed frame
        target_code = mat2py_bodn2c(target)
        recpoint = mat2py_srfrec(target_code[0], lon_rad, lat_rad)

        # Get the spacecraft position in the body-fixed frame
        instpos, _ = mat2py_spkpos(sc, t, targetframe, abcorr, target)

        # Compute the vector from the instrument to the target point
        v1 = recpoint - instpos

        # Flatten v1 to convert it to shape (3,)
        v1 = v1.flatten()

        # Normalize v1 to get the boresight vector
        boresight = v1 / np.linalg.norm(v1)

        # The z-axis of the instrument is aligned with the boresight
        rotmat[:, 2] = boresight

        # Define a consistent reference vector to avoid singularities
        reference_vector = np.array([0, 0, 1])

        # Check for alignment or anti-alignment with the reference vector
        if abs(np.dot(boresight, reference_vector)) > 0.999:
            # Adjust reference vector if necessary
            reference_vector = np.array([0, 1, 0])

        # Compute the y-axis using the cross product
        yinst = np.cross(boresight, reference_vector)
        yinst /= np.linalg.norm(yinst)  # Normalize

        # Compute the x-axis
        xinst = np.cross(yinst, boresight)
        xinst /= np.linalg.norm(xinst)  # Normalize

        # Complete the rotation matrix
        rotmat[:, 0] = xinst
        rotmat[:, 1] = yinst
        rotmat[:, 2] = boresight

    else:
        # Pointing constrained by kernel (CKernel)
        # Equivalent to MATLAB: rotmat = cspice_pxform(instframe, targetframe, t)
        rotmat = mat2py_pxform(instframe, targetframe, t)

        # Compute the surface intercept point
        # Equivalent to MATLAB: [xpoint, ~, ~, found] = cspice_sincpt(...)
        xpoint, trgepc, srfvec, found = mat2py_sincpt(method, target, t, targetframe, abcorr, sc, instframe, boresight)

        if found:
            # Convert rectangular coordinates to latitude and longitude
            _, lon_rad, lat_rad = mat2py_reclat(xpoint)
            lon_out = lon_rad * mat2py_dpr()  # Degrees per radian
            lat_out = lat_rad * mat2py_dpr()
        else:
            # Point is not visible
            print(f"On {mat2py_et2utc(t, 'C', 0)}, {inst} is not pointing at {target}")
            return fovbounds, boresight, rotmat, visible, lon_out, lat_out

        # The boresight of the instrument must point at the target point
        recpoint = xpoint  # Rectangular coordinates of the target point

        # Get the spacecraft position in the body-fixed frame
        instpos, _ = mat2py_spkpos(sc, t, targetframe, abcorr, target)

        # Compute the vector from the instrument to the target point
        v1 = recpoint - instpos

    # Transform the FOV bounds to the target frame
    for i in range(bounds.shape[1]):
        fovbounds[:, i] = rotmat @ bounds[:, i]

    # Check if the point is visible from the instrument
    if np.dot(v1, recpoint) > 0:
        # Point is not visible; return with visible set to False
        return fovbounds, boresight, rotmat, visible, lon_out, lat_out
    else:
        visible = True  # Point is visible

    # Output values
    return fovbounds, boresight, rotmat, visible, lon_out, lat_out
