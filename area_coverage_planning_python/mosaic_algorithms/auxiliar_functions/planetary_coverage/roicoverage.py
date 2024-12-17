import numpy as np
from shapely.geometry import Polygon, MultiPolygon
#from mosaic_algorithms.auxiliar_functions.planetary_coverage.regionarea import regionarea
from mosaic_algorithms.auxiliar_functions.polygon_functions.amsplit import amsplit
from mosaic_algorithms.auxiliar_functions.polygon_functions.interppolygon import interppolygon


def roicoverage(target, roi, fplist):

    """
     Provided a list of footprints on a body surface and ROIs, this function
     computes the cumulative coverage and overlap on each ROI, in [%]
     Overlap accounts for the % of surface covered at least more than once
     [Note:] This needs a revision, only works for a single roi!!

     Programmers:  Paula Betriu (UPC/ESEIAAT)
     Date:         12/2023
     Revision:     1

     Usage:        [coverage, overlap] = roicoverage(target, roi, fplist)

     Inputs:
       > target:   SPICE ID (int) or name (string) of the body
       > roi:      roi struct that contains, at least, its boundary vertices,
                   in latitudinal coordinates (in [deg])
       > fplist:   list of footprints (struct). See footprint for further
                   information

     Output:
       > coverage: Percentage of surface that the list of provided footprints
                   collectively cover with respect to the ROI surface,
                   in [%]
       > overlap:  Percentage of overlap of the list of provided footprints
                   with respect to the ROI surface, in [%]
    """

    # Pre-allocate variables

    coverage = 0  # initialize coverage
    overlap = 0  # initialize overlap

    if not fplist:
        return coverage, overlap

    # In case roi is not input as a struct but as a matrix (vertices)...
    if not isinstance(roi, dict):
        col1, col2 = amsplit(roi[:, 0], roi[:, 1])
        vertices = np.hstack((col1.reshape(len(col1), 1), col2.reshape(len(col2), 1)))
        del roi
        roi = {'vertices': interppolygon(vertices)}
    else:
        roi['vertices'] = np.array(roi['vertices'])
        col1, col2 = amsplit(roi['vertices'][:, 0], roi['vertices'][:, 1])
        roi['vertices'] = interppolygon(np.hstack((col1.reshape(len(col1), 1), col2.reshape(len(col2), 1))))

    # Cumulative coverage
    for footprint in fplist:
        if 'polyJ' not in roi:
            roi['polyJ'] = Polygon()
            roi['polyO'] = Polygon()

        # Footprint polyshape
        if len(footprint['bvertices']) == 0:
            continue

        x, y = footprint['bvertices'][:, 0], footprint['bvertices'][:, 1]

        if (np.isnan(x)).any():
            nanindex = np.where(np.isnan(x))[0]
            polygon_list = []
            for i in range(len(nanindex)):
                if i == 0:
                    polygon_list.append(Polygon(list(zip(x[:nanindex[0]], y[:nanindex[0]]))))
                else:
                    polygon_list.append(Polygon(
                        list(zip(x[nanindex[i - 1] + 1:nanindex[i]], y[nanindex[i - 1] + 1:nanindex[i]]))))
            if ~ np.isnan(x[-1]):
                polygon_list.append(Polygon(list(zip(x[nanindex[-1] + 1:], y[nanindex[-1] + 1:]))))
            polyFP = MultiPolygon(polygon_list)
        else:
            polyFP = Polygon(list(zip(x,y)))
        polyFP = polyFP.buffer(0)

        # Intersect footprint to cumulative coverage
        polyI = roi['polyJ'].intersection(polyFP)
        roi['polyO'] = roi['polyO'].union(polyI)

        # Join footprint
        roi['polyJ'] = roi['polyJ'].union(polyFP)

    # Calculate ROI specific coverage
    # Get total surface of the ROI
    lon, lat = roi['vertices'][:, 0], roi['vertices'][:, 1]
    # RA, _ = regionarea(target, lon, lat)

    if (np.isnan(lon)).any():
        nanindex = np.where(np.isnan(lon))[0]
        polygon_list = []
        for i in range(len(nanindex)):
            if i == 0:
                polygon_list.append(Polygon(list(zip(lon[:nanindex[0]], lat[:nanindex[0]]))))
            else:
                polygon_list.append(Polygon(
                    list(zip(lon[nanindex[i - 1] + 1:nanindex[i]], lat[nanindex[i - 1] + 1:nanindex[i]]))))
        if ~ np.isnan(lon[-1]):
            polygon_list.append(Polygon(list(zip(lon[nanindex[-1] + 1:], lat[nanindex[-1] + 1:]))))
        polyROI = MultiPolygon(polygon_list)
    else:
        polyROI = Polygon(list(zip(lon, lat)))
    polyROI = polyROI.buffer(0)

    # Intersect ROI with cumulative footprints
    if (isinstance(roi['polyJ'], Polygon) and roi['polyJ'].exterior.coords) or isinstance(roi['polyJ'], MultiPolygon):
        polyI = polyROI.intersection(roi['polyJ'])
        if (isinstance(polyI, Polygon) and polyI.exterior.coords) or isinstance(polyI, MultiPolygon):
            inter = polyROI.difference(polyI)
            coverage = ((polyROI.area - inter.area) * 100) / polyROI.area

    # Intersect ROI with overlap
    if (isinstance(roi['polyO'], Polygon) and roi['polyO'].exterior.coords) or isinstance(roi['polyO'], MultiPolygon):
        polyI = polyROI.intersection(roi['polyO'])
        if (isinstance(polyI, Polygon) and polyI.exterior.coords) or isinstance(polyI, MultiPolygon):
            inter = polyROI.difference(polyI)
            overlap = ((polyROI.area - inter.area) * 100) / polyROI.area

    return coverage, overlap