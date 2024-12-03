from mosaic_algorithms.paper.figure3.input_data_fig3 import roistruct
from mosaic_algorithms.auxiliar_functions.multiprocess.spawnProcess import *


nroi = len(roistruct) # total number of ROIs to be observed
nroiprocess = 1 # number of ROIs per process
runName = 'main_onlineFrontier_spawn'

sp = spawnProcess(nroi, runName, nroiprocess)
sp.spawn()

