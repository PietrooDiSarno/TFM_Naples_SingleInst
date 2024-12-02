import copy
from conversion_functions import *
import numpy as np
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.instpointing import  instpointing
from mosaic_algorithms.auxiliar_functions.plot.trgobsvec import trgobsvec
from mosaic_algorithms.auxiliar_functions.observation_geometry.emissionang import emissionang
from mosaic_algorithms.auxiliar_functions.polygon_functions.sortcw import sortcw

def topo2inst(inputdata_, lon, lat, target, sc, inst, et):
    """
    This function transforms a set of points from the topographic coordinate
    system(latitude and longitude on the target body) to the instrument frame
    coordinates.

    Programmers:  Paula Betriu (UPC/ESEIAAT)
    Date:         09/2022

    Usage:        outputData = topo2inst(inputdata, lon, lat, target, sc, inst, et)

    Inputs:
      > inputdata:     A list of lists  or ndarray of points in topographic
                       coordinates to be transformed. Each point is a row with
                       [longitude, latitude] format
      > lon:           longitude of the observation point or area center, in
                       [deg]
      > lat:           latitude of the observation point or area center, in
                       [deg]
      > target:        string name of the target body
      > sc:            string name of the spacecraft
      > inst:          string name of the instrument
      > et:            ephemeris time, TDB seconds past J2000 epoch

    Outputs:
      > outputData:    A list of lists or ndarray of the input points transformed
                       to the instrument frame coordinates. The format of the
                       output matches the input (cell array or matrix)

    """
    inputdata = copy.deepcopy(inputdata_)
    # Handle input data in list format, ensuring all empty entries are replaced
    # with [NaN, NaN]
    ii, jj,aux = [], [], []

    if isinstance(inputdata, list):
        aux = [[point if point else [np.nan, np.nan] for point in row] for row in inputdata]
        topoPoints = np.vstack([[point for point in sublist] for sublist in aux])
        for i in range(len(inputdata)):
            for j in range(len(inputdata[i])):
                ii.append(i)
                jj.append(j)
    else:
        if np.shape(inputdata) == (2,):
            inputdata = inputdata.reshape(1,2)
        topoPoints = copy.deepcopy(inputdata)

    # Pre-allocate variables
    _, targetframe,_ = mat2py_cnmfrm(target)  # target frame ID in SPICE

    # Build focal plane
    fovbounds, boresight, rotmat ,_ = instpointing(inst, target, sc, et, lon, lat)

    vertex,_ = mat2py_spkpos(sc, et, targetframe, 'NONE', target)
    point = vertex + fovbounds[:, 0]

    # Create a plane based on the boresight and a point in the focal plane
    plane = mat2py_nvp2pl (boresight, point)

    # For each topographic point, find its intersection with the focal plane
    spoint = np.zeros((topoPoints.shape[0], 3))

    for i in range(topoPoints.shape[0]):
        if not np.isnan(topoPoints[i]).any():
            dir = -(trgobsvec(topoPoints[i], et, target, sc))[0]
            found, spoint[i, :] = mat2py_inrypl(vertex, dir, plane)
            if not found:
                print('No intersection')
        else:
            spoint[i, :] = np.full(3, np.nan)

    # Transform coordinates from body-fixed to instrument frame
    tArea = np.zeros([spoint.shape[0],3])
    for i in range(spoint.shape[0]):
        if not np.isnan(spoint[i,:]).any():
            vpoint = -(vertex - spoint[i,:].T) # vector from spacecraft to intersection point
            tArea[i, :] = np.linalg.inv(rotmat).dot(vpoint) # apply inverse rotation to transform
                                                              # to instrument frame
        else:
            tArea[i, :] = np.full(3, np.nan)

    instcoord = tArea[:, :2]  # extract 2D instrument frame coordinates
    # Prepare output data matching the format of the input,i.e., cell array or
    # matrix
    outputData = [[None for _ in row] for row in inputdata]
    if isinstance(inputdata, list):
        for k in range(len(ii)):
            if not np.isnan(instcoord[k]).any():
                outputData[ii[k]][jj[k]] = instcoord[k, :]
    else:
        outputData = copy.deepcopy(instcoord)
    return outputData