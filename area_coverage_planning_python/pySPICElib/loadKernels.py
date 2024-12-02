from pySPICElib import *
from conversion_functions import *
import os

def loadKernels(metakr):


    # Clean pool of kernels
    endSPICE()

    #Get homespice
    #HOMESPICE = getHomeSpice
    #kernelpath_ = os.path.join(getHomeSpice, 'kernels')

    kernelpath_ = 'C:\\Users\\kekka\\Documents\\SPICE\\kernels'

    # load kernels
    kf=kernelFetch(kernelPath_= kernelpath_)
    kf.ffList(metakr)


    # Total loaded kernels
    print(f"Kernel pool: {mat2py_ktotal('ALL')}")  # in MATLAB this is done in the function loadKernels