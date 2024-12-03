
import importlib
import sys
from mosaic_algorithms.auxiliar_functions.multiprocess.dataHandling import dataHandling
from mosaic_algorithms.online_frontier_repair.frontierRepair import frontierRepair
from mosaic_algorithms.paper.figure3.input_data_fig3 import *  # Load mission info (kernels, SPICE ids, etc.)


# Revision of grid discretization:
# Grid is going to be built in the camera frame, instead of the body-fixed
# frame. This is somewhat more complicated (it requires more calculations)
# but it could correct the spatial aberration that we presently see
# Re-implementation of sidewinder function with these new feature (and
# other code improvements)

mosaic = 'onlinefrontier'
process = int(sys.argv[1])
nroi = int(sys.argv[2])
rem = int(sys.argv[3])
saveplot = True

if nroi!=rem:
    firstroi = (process-1)*nroi
    lastroi = process*nroi
else:
    firstroi = len(roistruct)-nroi
    lastroi = len(roistruct)

dh = dataHandling()
npoints = 30
for i in range(firstroi,lastroi):
    roiname = roistruct[i]['name'].lower().replace(" ", "")
    roi = np.array(roistruct[i]['vertices'])
    times = np.linspace(roistruct[i]['inittime'], stoptime, npoints)
    makespan = []

    for init_time in times:
        # Online Frontier
        A, fpList = frontierRepair(init_time, stoptime, tcadence, inst, sc, target, roi, olapx, olapy, 3 * 1e-3)
        if not fpList == []:
            makespan.append(fpList[-1]['t'] + tcadence - init_time)
        else:
            makespan.append(None)

    dh.saveValues(mosaic, roiname, times, makespan)
    if saveplot:
           dh.savePlots(mosaic, roiname, times, makespan)

    #a,b,c = dh.getValues(mosaic, roiname)






