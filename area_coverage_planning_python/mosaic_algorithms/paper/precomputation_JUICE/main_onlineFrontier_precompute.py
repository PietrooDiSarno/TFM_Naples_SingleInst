from area_coverage_planning_python.mosaic_algorithms.auxiliar_functions.multiprocess.dataHandling import dataHandling
from area_coverage_planning_python.mosaic_algorithms.online_frontier_repair.frontierRepair import frontierRepair
from area_coverage_planning_python.mosaic_algorithms.auxiliar_functions.spacecraft_operation.computeResMosaic import computeResMosaic
def mosaicOnlineFrontier(timeint, inst, sc, ROIob, instrument, int):

    # Revision of grid discretization:
    # Grid is going to be built in the camera frame, instead of the body-fixed
    # frame. This is somewhat more complicated (it requires more calculations)
    # but it could correct the spatial aberration that we presently see
    # Re-implementation of sidewinder function with these new feature (and
    # other code improvements)

    mosaic = 'onlinefrontier'
    saveplot = True

    dh = dataHandling()
    npoints = len(timeint)

    makespan = []
    nImg = []
    resROI = []


    roiname = ROIob.name
    roi = ROIob.vertices
    target = ROIob.body
    stoptime = timeint[-1] + 3600

    tcadence = 8.5  # [s] between observations
    olapx = 20  # [%] of overlap in x direction
    olapy = 20  # [%] of overlap in y direction


    for init_time in timeint:
        # Online Frontier
        A, fpList = frontierRepair(init_time, stoptime, tcadence, inst, sc, target, roi, olapx, olapy, 3 * 1e-3)
        if not fpList == []:
            makespan.append(fpList[-1]['t'] + tcadence - init_time)
            nImg.append(len(fpList))
            resROI.append(computeResMosaic(fpList, instrument.ifov))
        else:
            makespan.append(None)
            nImg.append(None)
            resROI.append(None)

    #dh.saveValues(mosaic, roiname, timeint, makespan, int)
    if saveplot:
        dh.savePlots(mosaic, roiname, timeint, makespan, nImg, int)

    # a,b,c,d = dh.getValues(mosaic, roiname, int)
    return makespan, nImg, resROI