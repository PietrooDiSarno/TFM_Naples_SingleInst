import numpy as np

from conversion_functions.mat2py_cnmfrm import mat2py_cnmfrm
from conversion_functions.mat2py_getfov import mat2py_getfov
from conversion_functions.mat2py_bodn2c import mat2py_bodn2c
from conversion_functions.mat2py_rpd    import mat2py_rpd
from conversion_functions.mat2py_srfrec import mat2py_srfrec
from conversion_functions.mat2py_spkpos import mat2py_spkpos
from conversion_functions.mat2py_pxform import mat2py_pxform
from conversion_functions.mat2py_sincpt import mat2py_sincpt
from conversion_functions.mat2py_reclat import mat2py_reclat
from conversion_functions.mat2py_dpr    import mat2py_dpr
from conversion_functions.mat2py_et2utc import mat2py_et2utc

def instpointing(inst, target, sc, t, *args):
    """
    This function sets the instrument's orientation, provided a target and
    the latitudinal coordinates of the point the instrument should be aiming
    at. It also checks if this point is actually visible from the FOV (could
    be on the dark side of the object as seen from the instrument).

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         01/2023

    Usage:        [fovbounds, boresight, rotmat, visible] = instorient(inst,
                      target, sc, t, lon, lat)
                  [fovbounds, boresight, rotmat, visible, lon, lat] = instorient(inst,
                      target, sc, t)

    Inputs:
      > inst:       string name of the instrument
      > target:     string name of the target body
      > sc:         string name of the spacecraft
      > t:          observation time, i.e. the minimum time that the
                    instrument needs to perform an observation, in seconds
      > lon:        longitude coordinate of the target body at which the
                    instrument boresight is pointing, in [deg]
      > lat:        latitude coordinate of the target body at which the
                    instrument boresight is pointing, in [deg]

    Returns:
      > fovbounds:  FOV bounds in the body-fixed reference frame centered at
                    the spacecraft position at time t
      > boresight:  FOV boresight in the body-fixed reference frame
                    centered at the spacecraft position at time t
      > rotmat:     rotation matrix (numpy.ndarray) from instrument frame to target frame
      > visible:    boolean that determines if the point is visible from the
                    instrument's FOV
      > lon:        longitude coordinate of the target body at which the
                    instrument boresight is pointing, in [deg]
      > lat:        latitude coordinate of the target body at which the
                    instrument boresight is pointing, in [deg]
    """
    #Pre-allocate variables
    axis3=False # boolean variable that indicates if the spacecraft is 3-axis steerable
    lon=0
    lat=0
    if len(args) > 0:
        lon = args[0]
        lat = args[1]
        axis3 = True

    method = 'ELLIPSOID' #assumption: ray intercept function is going to
    #model the target body as a tri-axial ellipsoid
    _,targetframe,_ = mat2py_cnmfrm(target) #target frame ID in SPICE
    abcorr = 'LT' #one-way light time aberration correction parameter.

    # Retrieve FOV parameters
    shape,instframe,boresight,bounds = mat2py_getfov((mat2py_bodn2c(inst))[0],4) # instrument FOV's boundary
    # vectors in the instrument frame
    if shape in ["CIRCLE", "ELLIPSE"]:
        raise ValueError("Circular and ellipsoidal FOV shapes have not been implemented yet")

    fovbounds = np.zeros((3, max(np.shape(bounds))))
    rotmat = np.zeros((3, 3))
    visible = False

    # Pointing matrix
    if axis3:  # 3-axis steerable
        lon = lon*mat2py_rpd()
        lat = lat*mat2py_rpd()

        # Boresight of the instrument must point at the target point
        recpoint = mat2py_srfrec(mat2py_bodn2c(target)[0], lon,lat)  # rectangular
        # coordinates of the target point in the body-fixed reference frame

        instpos,_ = mat2py_spkpos(sc, t, targetframe, abcorr,target) # rectangular coordinates
        # of the instrument in the body-fixed reference frame
        v1 = recpoint - instpos  # distance vector to the target point from the instrument in the body-fixed reference frame
        boresight = v1 / np.linalg.norm(v1)  # boresight of the instrument in the
        # body-fixed reference frame

        # The z-axis of the instrument is the boresight
        rotmat[:, 2] = boresight

        # Define a consistent reference vector, e.g., [0, 0, 1] (celestial north)
        reference_vector = np.array([0, 0, 1])

        # Check if boresight is aligned or anti-aligned with reference vector
        if abs(np.dot(boresight, reference_vector)) > 0.999:
            # Adjust reference vector if aligned/anti-aligned to avoid singularity
            reference_vector = np.array([0, 1, 0])

        # Define y-axis using cross product to ensure perpendicularity
        yinst = np.cross(boresight, reference_vector)
        yinst = yinst/np.linalg.norm(yinst)  # normalize yinst

        # Define x-axis using cross product between boresight and yinst
        xinst = np.cross(yinst, boresight)
        xinst = xinst/np.linalg.norm(xinst)  # normalize xinst

        # Assign to rotation matrix
        rotmat[:, 0] = xinst
        rotmat[:, 1] = yinst
        rotmat[:, 2] = boresight

    else:  # Pointing constrained by ckernel
        rotmat = mat2py_pxform(instframe, targetframe, t)
        xpoint, _, _, found = mat2py_sincpt(method, target, t, targetframe, abcorr, sc, instframe, boresight)

        if found:
            _, lon, lat = mat2py_reclat(xpoint)
            lon = lon*mat2py_dpr()
            lat = lat*mat2py_dpr()
        else:
            # not visible
            print(f"On {mat2py_et2utc(t, 'C', 0)}, {inst} is not pointing at {target}")
            if len(args) > 0:
                return fovbounds, boresight, rotmat, visible
            else:
                return fovbounds, boresight, rotmat, visible, lon, lat

        # Boresight of the instrument must point at the target point
        recpoint = xpoint  # rectangular coordinates of the target point in the body-fixed reference frame
        instpos, _ = mat2py_spkpos(sc, t, targetframe, abcorr,
                                  target)  # rectangular coordinates of the instrument in the body-fixed reference frame
        v1 = recpoint - instpos  # distance vector to the target point from the
        #instrument in the body-fixed reference frame

    # Transform coordinates
    fovbounds = np.zeros([3,max(np.shape(bounds))])
    for i in range(max(np.shape(bounds))):
        fovbounds[:, i] = np.dot(rotmat, bounds[:, i])  # instrument FOV's boundary vectors in the target frame

    # Check if the point is visible as seen from the instrument
    if np.dot(v1, recpoint) > 0:  # check if the point is visible as seen from the instrument
        if len(args) > 0:
            return fovbounds, boresight, rotmat, visible
        else:
            return fovbounds, boresight, rotmat, visible, lon, lat
    else:
        visible = True

    # Output values
    if len(args) > 0:
        return fovbounds, boresight, rotmat, visible
    else:
        return fovbounds, boresight, rotmat, visible, lon, lat
