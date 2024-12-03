import numpy as np
from conversion_functions import *
from shapely.geometry import Polygon
from pySPICElib import *

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
tcadence = 15.  # [s] between observations
olapx = 20.  # [%] of overlap in x direction
olapy = 20.  # [%] of overlap in y direction
speedUp = 0.
count = 0

# Regions of interest
count += 1
roi = np.array([[25., 25.],
                [25., -5.],
                [-5., -5.],
                [-5., 25.]])
roistruct = [{} for _ in range(count)]  # Pre-allocate roistruct
roistruct[count - 1]['vertices'] = roi.tolist()
polygon = Polygon(roi)
cx, cy = polygon.centroid.x, polygon.centroid.y
roistruct[count - 1]['cpoint'] = np.array([cx, cy])
roistruct[count - 1]['name'] = "ROI"
roistruct[count - 1]['inittime'] = mat2py_str2et('1998 MAR 29 12:35:00.000 TDB')  # closest approach
