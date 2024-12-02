import numpy as np
from shapely.geometry import Polygon, MultiPolygon

# Import the SPICE wrapper functions
from conversion_functions import *
from mosaic_algorithms.auxiliar_functions.polygon_functions.amsplit import amsplit


def visibleroi(roi, et, target, obs):
    """
    Given a target area (region-of-interest) on a planetary body surface and
    an observer, this function calculates the portion of the former that is
    visible.
    Note: if ROI intercepts the anti-meridian line, this function returns
    the visible ROI polygon divided by this line (see 'amsplit' function).

    Usage: vroi, inter = visibleroi(roi, et, target, obs)

    Inputs:
        roi:    Matrix containing the vertices of the ROI polygon. The
                vertex points are expressed in 2D, in latitudinal
                coordinates [degrees].
            roi[:,0] corresponds to the x values (longitudes) of the vertices
            roi[:,1] corresponds to the y values (latitudes) of the vertices
        et:     Limb projection time, in seconds past J2000 epoch
        target: String name of the target body (SPICE ID)
        obs:    String name of the observer (SPICE ID)

    Outputs:
        vroi:   Matrix containing the vertices of the intersection between
                the input ROI polygon and the limb projection, i.e., the
                visible portion of the ROI on the body surface as seen from
                the observer
        inter:  The intersection polygon (shapely.geometry.Polygon)
    """
    # Previous anti-meridian intersection check
    # Find the discontinuity in the longitude coordinates
    lon_sorted = np.sort(roi[:, 0])
    lon_diff = np.diff(lon_sorted)
    # Find indices where the difference is >= 180 degrees
    discontinuity_indices = np.where(lon_diff >= 180)[0]
    if discontinuity_indices.size > 0:
        # Split the ROI polygon at the anti-meridian
        roi_lon, roi_lat = amsplit(roi[:, 0], roi[:, 1])
        roi = np.column_stack((roi_lon, roi_lat))

    # Parameters for 'cspice_limbpt' function
    method = 'TANGENT/ELLIPSOID'
    # Get the body-fixed frame of the target
    code, targetframe, found = mat2py_cnmfrm(target)  # body-fixed frame
    targetframe = targetframe[0][0]
    abcorr = 'XLT+S'
    corloc = 'CENTER'
    refvec = np.array([0.0, 0.0, 1.0])  # First of the sequence of cutting half-planes
    ncuts = int(1e3)  # Number of cutting half-planes
    # Angular step by which to roll the cutting half-planes about the observer-target vector
    delrol = mat2py_twopi() / ncuts
    schstp = 1.0e-4  # Search angular step size
    soltol = 1.0e-7  # Solution convergence tolerance

    # Limb calculation with 'cspice_limbpt' function
    _, limb_points, _, _ = mat2py_limbpt(
        method, target, et, targetframe, abcorr, corloc, obs,
        refvec, delrol, ncuts, schstp, soltol, ncuts
    )

    # Convert limb points from rectangular to latitudinal coordinates
    _, lblon, lblat = mat2py_reclat(limb_points)
    lblon = lblon * mat2py_dpr()
    lblat = lblat * mat2py_dpr()

    # Check for anti-meridian split in limb points
    lon_sorted = np.sort(lblon)
    lon_diff = np.diff(lon_sorted)
    discontinuity_indices = np.where(lon_diff >= 180)[0]
    if discontinuity_indices.size > 0:
        # Split the limb points at the anti-meridian
        lblon, lblat = amsplit(lblon, lblat)

    # Create polygons and compute the intersection
    # Create limb polygon (poly1) using the limb points
    poly1 = Polygon(np.column_stack((lblon, lblat)))
    # Create ROI polygon (poly2) using the ROI vertices
    poly2 = Polygon(roi)
    # Compute the intersection of the two polygons
    inter = poly1.intersection(poly2)

    # Output visible ROI
    # If the intersection is empty, return empty arrays
    if inter.is_empty:
        vroi = np.array([])
    else:
        # Handle cases where the intersection is a MultiPolygon (multiple disconnected regions)
        if isinstance(inter, MultiPolygon):
            # Combine the coordinates of all polygons
            vroi = []
            for geom in inter.geoms:
                x, y = geom.exterior.coords.xy
                vroi.append(np.column_stack((x, y)))
            # Concatenate all coordinates into a single array
            vroi = np.vstack(vroi)
        else:
            # The intersection is a single polygon
            x, y = inter.exterior.coords.xy
            vroi = np.column_stack((x, y))

    return vroi, inter
