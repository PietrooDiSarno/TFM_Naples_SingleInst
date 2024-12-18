from area_coverage_planning_python.mosaic_algorithms.auxiliar_functions.multiprocess.dataHandling import dataHandling
from area_coverage_planning_python.mosaic_algorithms.online_frontier_repair.frontierRepair import frontierRepair
from area_coverage_planning_python.mosaic_algorithms.auxiliar_functions.spacecraft_operation.computeResMosaic import computeResMosaic
from mosaic_algorithms.auxiliar_functions.planetary_coverage.roicoverage import roicoverage


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
    cov = []

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
            roi_ = {'vertices': roi}
            cov_ = roicoverage(target, roi_, fpList)[0]
            if cov_ > 95:
                cov.append(100)
            else:
                cov.append(cov_)
        else:
            makespan.append(None)
            nImg.append(None)
            resROI.append(None)
            cov.append(None)


    return makespan, nImg, resROI, cov