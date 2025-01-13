import copy
import os
import pickle

from mosaic_algorithms.auxiliar_functions.planetary_coverage.roicoverage import roicoverage
from mosaic_algorithms.online_frontier_repair.frontierRepair import frontierRepair
from mosaic_algorithms.paper.video.input_data_juice import *
import matplotlib.pyplot as plt
import importlib
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.computeResMosaic import computeResMosaic

# Revision of grid discretization:
# Grid is going to be built in the camera frame, instead of the body-fixed
# frame. This is somewhat more complicated (it requires more calculations)
# but it could correct the spatial aberration that we presently see
# Re-implementation of sidewinder function with these new feature (and
# other code improvements)


roiname = roistruct[0]['name'].lower().replace(" ", "")

first_mosaic = mat2py_str2et('2034 JUN 06 06:45:30.000 UTC')
last_mosaic = mat2py_str2et('2034 JUN 06 12:00:00.000 UTC')
time_vec = np.linspace(first_mosaic, last_mosaic, 5, endpoint=True)

for i, init_time in enumerate(time_vec):
    # Online Frontier
    A, fpList = frontierRepair(init_time, stoptime, tcadence, inst, sc, target, roi, olapx, olapy, 3 * 1e-3)
    with open (os.path.join('pickles_video', f'pickle_{roiname}_{i}.cfg'), 'wb') as f:
        pickle.dump([A,fpList],f)
        f.close()



