import spiceypy as spice
import matplotlib.pyplot as plt
from pySPICElib.SPICEtools import *

def plotSchedule(bestI, roiL, flybys):
    #col = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
    bestI.stol = bestI.stol[:7]
    roiL = roiL[:7]
    #color = dict(zip(roiL[:].name, col))
    H = 0.5  # just to plot tw, height of each roi tw
    for fb in flybys:
        fbs = spice.str2et(fb[0])
        fbe = spice.str2et(fb[1])
        fig, ax = plt.subplots()
        for i in range(len(bestI.stol)):
            roi = roiL[i]
            for interval in range(spice.wncard(roi.ROI_TW)):
                s,e = spice.wnfetd(roi.ROI_TW, interval)
                if s <= bestI.stol[i] <= e and fbs <= bestI.stol[i] <= fbe:
                    ax.plot(roi.ROI_ObsET[interval], roi.ROI_ObsRes[interval])
                    obs_tw = newTimeWindow(bestI.stol[i], bestI.stol[i] + bestI.obsLength[i])
                    plot_tw(ax, obs_tw, 0 + i * H, 0 + H * (i + 1), 'r')  # plot observation length
    plt.show()