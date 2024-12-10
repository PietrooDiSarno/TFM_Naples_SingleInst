import numpy as np
from conversion_functions import *
from shapely.geometry import Polygon
from pySPICElib.kernelFetch import kernelFetch

# Mosaic comparison between the different heuristics

# Relevant paths
kernelpath = '../../../input'

# Case study
# Define mission and spacecraft (SPICE ID of the spacecraft)
mission = 'GALILEO'
sc = 'GALILEO ORBITER'

# Choose instrument (SPICE ID of the instrument)
inst = 'GLL_SSI'

# Define the planetary body (SPICE ID)
target = 'EUROPA'

# Planetary body modelization (DSK or ELLIPSOID)
method = 'ELLIPSOID'

# Clean pool of kernels
mat2py_kclear()

# load kernels
kf = kernelFetch(textFilesPath_=f'{kernelpath}/{mission.lower()}/')
kf.ffFile(metaK='inputkernels.txt')

# Total loaded kernels
print(f"Kernel pool: {mat2py_ktotal('ALL')}")

# Definition of ROIs
# Pre-allocation of variables...
stoptime = mat2py_str2et('1998 MAY 30 00:00:00.000 TDB')  # mosaic end (max)
#stoptime = mat2py_str2et('1998 MAR 29 15:30:00.000 TDB')  # mosaic end (max)
tcadence = 8.5  # [s] between observations
olapx = 20  # [%] of overlap in x direction
olapy = 20  # [%] of overlap in y direction
speedUp = 0
count = 0


roistruct = []  # Pre-allocate roistruct

# Regions of interest

"""
# Pwyll Crater
count += 1
roi = np.array([
    [78, -18],
    [90, -14],
    [109, -13],
    [101, -23],
    [99, -30],
    [85, -33],
    [71, -28]
])

roistruct.append({'vertices': roi.tolist()})
polygon = Polygon(roi)
cx, cy = polygon.centroid.x, polygon.centroid.y
roistruct[count - 1]['cpoint'] = np.array([cx, cy])
roistruct[count - 1]['name'] = "Pwyll Crater"
roistruct[count - 1]['inittime'] = mat2py_str2et('1998 MAR 29 12:38:00.000 TDB')  # closest approach

"""
"""
# Annwn Regio [lon, lat] = [40, 20]º
count += 1
roi = np.array([
    [60, 30],
    [60, 10],
    [40, 10],
    [40, 30]
])
roistruct.append({'vertices': roi.tolist()})
polygon = Polygon(roi)
cx, cy = polygon.centroid.x, polygon.centroid.y
roistruct[count - 1]['cpoint'] = np.array([cx, cy])
roistruct[count - 1]['name'] = "Annwn Regio"
roistruct[count - 1]['inittime'] = mat2py_str2et('1998 MAR 29 12:51:00.000 TDB')  # closest approach
"""

"""
# Niamh
count += 1
roi = np.array([
    [150, 25],
    [150, 15],
    [135, 15],
    [135, 25]
])
roistruct.append({'vertices': roi.tolist()})
polygon = Polygon(roi)
cx, cy = polygon.centroid.x, polygon.centroid.y
roistruct[count - 1]['cpoint'] = np.array([cx, cy])
roistruct[count - 1]['name'] = "Niamh"
roistruct[count - 1]['inittime'] = mat2py_str2et('1998 MAR 29 13:29:00.000 TDB')
"""
#"""
# Cilix crater [lon, lat] = [180, 0]º
count += 1
roi = np.array([
    [-177, 3],
    [-177, -3],
    [177, -3],
    [177, 3]
])
roistruct.append({'vertices': roi.tolist()})
polygon = Polygon(roi)
cx, cy = polygon.centroid.x, polygon.centroid.y
roistruct[count - 1]['cpoint'] = np.array([cx, cy])
roistruct[count - 1]['name'] = "Cilix Crater"
roistruct[count - 1]['inittime'] = mat2py_str2et('1998 MAR 29 13:40:00.000 TDB')  # closest approach
#"""
"""
# Tara Regio
count += 1
roi = np.array([
    [-55, 20],
    [-95, 20],
    [-95, -20],
    [-55, -20]
])
roistruct.append({'vertices': roi.tolist()})
polygon = Polygon(roi)
cx, cy = polygon.centroid.x, polygon.centroid.y
roistruct[count - 1]['cpoint'] = np.array([cx, cy])
roistruct[count - 1]['name'] = "Tara Regio"
roistruct[count - 1]['inittime'] = mat2py_str2et('1998 MAR 29 14:00:00.000 TDB')
"""
"""
# Taliesin
count += 1
roi = np.array([
    [-160, -10],
    [-160, -30],
    [-130, -30],
    [-130, -10]
])
roistruct.append({'vertices': roi.tolist()})
polygon = Polygon(roi)
cx, cy = polygon.centroid.x, polygon.centroid.y
roistruct[count - 1]['cpoint'] = np.array([cx, cy])
roistruct[count - 1]['name'] = "Taliesin"
roistruct[count - 1]['inittime'] = mat2py_str2et('1998 MAR 29 14:21:00.000 TDB')
"""