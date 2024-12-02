import numpy as np
from conversion_functions import *
from mosaic_algorithms.auxiliar_functions.observation_geometry.emissionang_gpt import emissionang
from mosaic_algorithms.auxiliar_functions.plot.trgobsvec_gpt import trgobsvec
from mosaic_algorithms.auxiliar_functions.polygon_functions.sortcw import sortcw
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.instpointing_gpt import instpointing


def topo2inst(inputdata, lon, lat, target, sc, inst, et):
    """
    Transforms a set of points from the topographic coordinate system (latitude and longitude on the target body)
    to the instrument frame coordinates.

    Parameters:
        inputdata (list or numpy array): Points in topographic
            coordinates to be transformed. Each point is [longitude, latitude].
            Can be a list of lists (equivalent to MATLAB cell array) or a numpy array.
        lon (float): Longitude of the observation point or area center, in degrees.
        lat (float): Latitude of the observation point or area center, in degrees.
        target (str): Name of the target body.
        sc (str): Name of the spacecraft.
        inst (str): Name of the instrument.
        et (float): Ephemeris time, TDB seconds past J2000 epoch.

    Returns:
        outputData (list or numpy array): The input points transformed
            to the instrument frame coordinates. The format matches the input (list or array).
    """

    # Handle input data in list (cell array) format, ensuring all empty entries are replaced with [NaN, NaN]
    if isinstance(inputdata, list):
        aux = []
        ii = []
        jj = []
        for i in range(len(inputdata)):
            for j in range(len(inputdata[i])):
                if inputdata[i][j] is not None:
                    aux.append(inputdata[i][j])
                else:
                    aux.append([np.nan, np.nan])
                ii.append(i)
                jj.append(j)
        topoPoints = np.array(aux)
    else:
        # Assume inputdata is a numpy array
        topoPoints = inputdata

    # Pre-allocate variables
    # Get target frame name in SPICE
    _, targetframe, _ = mat2py_cnmfrm(target)
    targetframe = targetframe[0][0]

    # Build focal plane
    # Call instpointing to get fovbounds, boresight, rotmat
    fovbounds, boresight, rotmat, visible, lon_out, lat_out = instpointing(inst, target, sc, et, lon, lat)
    # Get spacecraft position relative to target
    vertex = mat2py_spkpos(sc, et, targetframe, 'NONE', target)
    vertex = vertex[0]
    vertex = vertex.flatten()  # Ensure vertex is a 1D array
    # Create a point in the focal plane
    point = vertex + fovbounds[:, 0]
    # Create a plane based on the boresight and a point in the focal plane
    plane = mat2py_nvp2pl(boresight, point)

    # For each topographic point, find its intersection with the focal plane
    spoint = np.zeros((topoPoints.shape[0], 3))
    for i in range(topoPoints.shape[0]):
        if not np.isnan(topoPoints[i, :]).any():
            # Compute the direction vector from target point to spacecraft
            dir_vector, _ = trgobsvec(topoPoints[i, :], et, target, sc)
            dir_vector = -dir_vector
            found, intersection_point = mat2py_inrypl(vertex, dir_vector, plane)
            if found:
                emnang = emissionang(topoPoints[i, :], et, target, sc)
                if emnang >= 90:
                    found = False
            if found:
                spoint[i, :] = intersection_point.flatten()
            else:
                spoint[i, :] = [np.nan, np.nan, np.nan]
        else:
            spoint[i, :] = [np.nan, np.nan, np.nan]

    # Transform coordinates from body-fixed to instrument frame
    tArea = np.zeros_like(spoint)
    for i in range(spoint.shape[0]):
        if not np.isnan(spoint[i, :]).any():
            vpoint = -(vertex - spoint[i, :])
            # Apply inverse rotation to transform to instrument frame
            tArea[i, :] = np.linalg.solve(rotmat, vpoint)
        else:
            tArea[i, :] = [np.nan, np.nan, np.nan]

    # Extract 2D instrument frame coordinates
    instcoord = tArea[:, 0:2]

    # Prepare output data matching the format of the input
    if isinstance(inputdata, list):
        outputData = [[None for _ in row] for row in inputdata]
        idx = 0
        for i, j in zip(ii, jj):
            if not np.isnan(instcoord[idx, :]).any():
                outputData[i][j] = instcoord[idx, :].tolist()
            else:
                outputData[i][j] = None
            idx += 1
    else:
        # Remove rows with NaNs
        valid_indices = ~np.isnan(instcoord[:, 0])
        instcoord = instcoord[valid_indices, :]
        # Sort coordinates clockwise
        aux_x, aux_y = sortcw(instcoord[:, 0], instcoord[:, 1])
        outputData = np.column_stack((aux_x, aux_y))

    return outputData

