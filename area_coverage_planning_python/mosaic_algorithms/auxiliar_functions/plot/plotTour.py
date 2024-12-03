import matplotlib.pyplot as plt
from shapely.geometry import MultiPolygon, Polygon
import numpy as np
from mosaic_algorithms.auxiliar_functions.plot.groundtrack import groundtrack
from mosaic_algorithms.auxiliar_functions.plot.mapPlot import mapPlot
from mosaic_algorithms.auxiliar_functions.polygon_functions.amsplit import amsplit


def plotTour(tour, fplist, roistruct, sc, target, *args):
    # Pre-allocate variables
    filename = f"../../auxiliar_functions/plot/global-maps/{target.lower()}-map.jpg"
    c1 = (0.72, 0.27, 1.00)
    c2 = (0.35, 0.06, 0.41)
    if len(args) > 0:
        ax = args[0]
        if len(args) > 1:
            video = args[1]
        else:
            video = None
    else:
        fig, ax = mapPlot(filename)
        video = None
    handles = []
    labels = []
    # Plot region(s) of interest
    if video is not None:
        for roi in roistruct:
            x = np.append(np.array(roi['vertices'])[:,0],np.array(roi['vertices'])[0,0])
            y = np.append(np.array(roi['vertices'])[:,1],np.array(roi['vertices'])[0,1])
            ax.plot(x,y,color='k', linewidth=1, linestyle='-', label=roi['name'])
            plt.legend()
            plt.pause(0.5)
        plt.pause(0.5)

    # Check
    if not tour:
        return ax

    # Plot footprint list
    for i in range(len(fplist)):
        # Plot footprint polygon
        if (np.isnan(fplist[i]['bvertices'][:, 0])).any():
            nanindex = np.where(np.isnan(fplist[i]['bvertices'][:, 0]))[0]
            for i in range(len(nanindex)):
                if i == 0:
                    h1 = ax.fill(fplist[i]['bvertices'][:nanindex[0], 0], fplist[i]['bvertices'][:nanindex[0],1],
                                 color=c1, alpha=0.8, edgecolor=c2, linewidth=1, label='Footprint')
                    #plt.pause(0.5)
                    if 'Footprint' not in labels:
                        handles.append(h1[0])
                        labels.append('Footprint')

                else:
                    h1 = ax.fill(fplist[i]['bvertices'][nanindex[i - 1] + 1:nanindex[i], 0],
                                 fplist[i]['bvertices'][nanindex[i - 1] + 1:nanindex[i], 1],
                                 color=c1, alpha=0.8, edgecolor=c2, linewidth=1, label='Footprint')
                    #plt.pause(0.5)
                    if 'Footprint' not in labels:
                        handles.append(h1[0])
                        labels.append('Footprint')
            if ~ np.isnan(fplist[i]['bvertices'][-1, 0]):
                h1 = ax.fill(fplist[i]['bvertices'][nanindex[-1] + 1:, 0],
                             fplist[i]['bvertices'][nanindex[-1] + 1:, 1],
                             color=c1, alpha=0.8, edgecolor=c2, linewidth=1, label='Footprint')
                #plt.pause(0.5)
                if 'Footprint' not in labels:
                    handles.append(h1[0])
                    labels.append('Footprint')
        else:
            h1 = ax.fill(fplist[i]['bvertices'][:, 0], fplist[i]['bvertices'][:, 1], color=c1, alpha=0.8, edgecolor=c2,
                         linewidth=1, label='Footprint')
            #plt.pause(0.5)
            if 'Footprint' not in labels:
                handles.append(h1[0])
                labels.append('Footprint')

        # Plot coverage path
        if i > 0:
            if abs(tour[i-1][0] - tour[i][0]) <= 180:  # no coverage path - a.m. intercept
                h2 = ax.plot([tour[i-1][0], tour[i][0]], [tour[i-1][1], tour[i][1]],
                        color='y', linewidth=1.5, label='Coverage path')
                #plt.pause(0.5)
                if 'Coverage path' not in labels:
                    handles.append(h2[0])
                    labels.append('Coverage path')
        else:
            h2 = None
            h3 = ax.scatter(tour[i][0], tour[i][1], s=100, color='b', marker='^', label='Start point')
            #plt.pause(0.5)
            if 'Start point' not in labels:
                handles.append(h3)
                labels.append('Start point')

        # Plot ground track
        t = fplist[i]['t']
        sclon, sclat = groundtrack(sc, t, target)
        if i > 0:
            h4 = ax.scatter(sclon, sclat, s=8, color='c', marker='o', label='Ground track')
            #plt.pause(0.5)
            if 'Ground track' not in labels:
                handles.append(h4)
                labels.append('Ground track')
        else:
            h4 = None
            ax.scatter(sclon, sclat, s=100, color='b', marker='^')
            #plt.pause(0.5)
        # Save animation
        if video is not None:
            plt.draw()
            video.grab_frame()

    # Re-plot the ROI (for aesthetic purposes)
    for i in range(len(roistruct)):
        roi = np.array(roistruct[i]['vertices'])
        x, y = amsplit(roi[:,0], roi[:,1])

        if (np.isnan(x)).any():
            nanindex = np.where(np.isnan(x))[0]
            polygon_list = []
            for i in range(len(nanindex)):
                if i == 0:
                    h5 = ax.plot(np.append(x[:nanindex[0]],x[0]), np.append(y[:nanindex[0]],y[0]),
                            color='k', linewidth=1, linestyle='-',label=roistruct[i]['name'])
                    #plt.pause(0.5)
                    if roistruct[i]['name'] not in labels:
                        handles.append(h5[0])
                        labels.append(roistruct[i]['name'])
                else:
                    h5 = ax.plot(np.append(x[nanindex[i - 1] + 1:nanindex[i]],x[nanindex[i - 1] + 1]),
                            np.append(y[nanindex[i - 1] + 1:nanindex[i]],y[nanindex[i - 1] + 1]),
                            color='k', linewidth=1, linestyle='-', label=roistruct[i]['name'])
                    #plt.pause(0.5)
                    if roistruct[i]['name'] not in labels:
                        handles.append(h5[0])
                        labels.append(roistruct[i]['name'])
            if ~ np.isnan(x[-1]):
                h5 = ax.plot(np.append(x[nanindex[-1] + 1:],x[nanindex[-1] + 1]),
                        np.append(y[nanindex[-1] + 1:],y[nanindex[-1] + 1]),
                        color='k', linewidth=1, linestyle='-', label=roistruct[i]['name'])
                #plt.pause(0.5)
                if roistruct[i]['name'] not in labels:
                    handles.append(h5[0])
                    labels.append(roistruct[i]['name'])
        else:
            h5 = ax.plot(np.append(x,x[0]), np.append(y,y[0]), color='k', linewidth=1, linestyle='-', label=roistruct[i]['name'])
            #plt.pause(0.5)
            if roistruct[i]['name'] not in labels:
                handles.append(h5[0])
                labels.append(roistruct[i]['name'])
   #plt.pause(0.5)
    # Legend

    unique_labels = dict(zip(labels, handles))
    legend = ax.legend(handles=unique_labels.values(), labels=unique_labels.keys(), loc='best', frameon=True)
    legend.get_frame().set_facecolor('white')
    legend.get_frame().set_alpha(1)
    plt.pause(3)
    if video is not None:
        video.finish()

    return ax