import numpy as np
from conversion_functions import *
from shapely.geometry import Polygon
from pySPICElib.kernelFetch import kernelFetch

# Mosaic comparison between the different heuristics

# Relevant paths
kernelpath = '../../../input'

# Case study
# Define mission and spacecraft (SPICE ID of the spacecraft)
mission = 'JUICE'
sc = 'JUICE'

# Choose instrument (SPICE ID of the instrument)
inst = 'JUICE_JANUS'

# Define the planetary body (SPICE ID)
target = 'GANYMEDE'

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
stoptime = mat2py_str2et('2034 JUN 07 07:02:00.000 UTC')   # mosaic end (max)
#stoptime = mat2py_str2et('1998 MAR 29 15:30:00.000 TDB')  # mosaic end (max)
tcadence = 8.5  # [s] between observations
olapx = 20  # [%] of overlap in x direction
olapy = 20  # [%] of overlap in y direction
speedUp = 0
count = 0


roistruct = []  # Pre-allocate roistruct

# Regions of interest

# JUICE_ROI_GAN_1_0_07: Galileo Regio
count += 1
roi = np.array([
    [247-360, 9],
    [261-360, 9],
    [261-360, 22],
    [247-360, 22]
])

roistruct.append({'vertices': roi.tolist()})
polygon = Polygon(roi)
cx, cy = polygon.centroid.x, polygon.centroid.y
roistruct[count - 1]['cpoint'] = np.array([cx, cy])
roistruct[count - 1]['name'] = "Selket"
roistruct[count - 1]['inittime'] = mat2py_str2et('2034 JUN 06 06:45:30.000 UTC')