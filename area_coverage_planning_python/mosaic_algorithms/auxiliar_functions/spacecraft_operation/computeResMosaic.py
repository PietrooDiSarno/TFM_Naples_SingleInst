from PSOA import pointres
import numpy as np
def computeResMosaic(fpList, ifov):
    r = []
    for fp in fpList:
        srfpoint = [fp['olon'], fp['olat']]
        t = fp['t']
        target = fp['target']
        obs = fp['sc']
        r.append(pointres(ifov, srfpoint, t, target, obs))

    res = np.mean(r)

    return res





