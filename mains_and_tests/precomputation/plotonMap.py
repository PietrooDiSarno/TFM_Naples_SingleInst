from matplotlib.patches import Polygon
import matplotlib.pyplot as plt
import math
import numpy as np
from pySPICElib.SPICEtools import plotGtrack
from PSOA.groundtrack import groundtrack
import spiceypy as spice
from pySPICElib.SPICEtools import newTimeWindow


def plotonMap(matrix, rois, fbs):
    myimage = '../../../data/mosaics/Ganymede_mosaic.jpg'
    minPos = []
    for row in matrix:
        row = [999 if math.isnan(x) else x for x in row]
        minPos.append(row.index(min(row)))
    transposed_matrix = [list(row) for row in zip(*matrix)]
    for fb, row in enumerate(transposed_matrix):
        ggtlon = []
        ggtlat = []
        fig, ax = plt.subplots()
        et0, et1 = fbs[fb]
        et = np.linspace(et0, et1, 500)
        dates = [spice.et2utc(et0,'C',0), spice.et2utc(et1,'C',0)]
        for t in et:
            gtlon, gtlat = groundtrack('JUICE', t, 'GANYMEDE')
            ggtlon.append(gtlon)
            ggtlat.append(gtlat)
        for j, res in enumerate(row):
            if math.isnan(res):
                continue
            else:
                if minPos[j] == fb:
                    color = 'red'
                else:
                    color = 'yellow'
                image = plt.imread(myimage)
                ax.imshow(image, interpolation='none', extent=[-180, 180, -90, 90])
                ax.set_ylim(ymin=-90, ymax=90)
                ax.set_xlim(xmin=-180, xmax=180)
                ax.set_xticks(np.arange(-180, 181, 30))
                ax.set_yticks(np.arange(-90, 91, 45))
                # ax.set_title('Observable ROIS during flyby ' + str(fb+1))

                if rois[j].vertices[0][0] > rois[j].vertices[1][
                    0]:  # if the longitude of the second vertices is lower than the longitude of the second vertice, the roi (assuming squared-shaped) crosses de antimeridian, thus we have an end-of-the-world roi
                    x0 = rois[j].vertices[0][0]
                    x1 = rois[j].vertices[3][0]
                    y0 = rois[j].vertices[0][1]
                    y1 = rois[j].vertices[3][1]
                    x01 = rois[j].vertices[1][0]
                    x11 = rois[j].vertices[2][0]
                    y01 = rois[j].vertices[1][1]
                    y11 = rois[j].vertices[2][1]
                    vertices1 = np.array([[x1, y1], [180., y1], [180., y0], [x0, y0]])
                    vertices2 = np.array([[x11, y11], [-180., y11], [-180., y01], [x01, y01]])

                    p1 = Polygon(vertices1, closed=True, fill=False)
                    p2 = Polygon(vertices2, closed=True, fill=False)
                    ax.add_patch(p1)
                    ax.add_patch(p2)
                    p1.set_color(color)
                    p1.set_linewidth(3)
                    p2.set_linewidth(3)
                    p2.set_color(color)
                    polygon_added = True

                else:
                    p = Polygon(rois[j].vertices, fill=False)
                    p.set_color(color)
                    p.set_linewidth(3)
                    ax.add_patch(p)
                    polygon_added = True
        plotGtrack(ax, ggtlon, ggtlat, label = 'ground track', color='cyan', lw=3)
        plt.title(dates[0] + ' - ' + dates[1] + ' (fb ' + str(fb) + ')')
        plt.legend()
        plt.gcf().set_dpi(200)
    plt.show()
    # if not polygon_added: plt.close(fig)
    # plt.savefig('../obs_plots/mosaic' + str(k))
