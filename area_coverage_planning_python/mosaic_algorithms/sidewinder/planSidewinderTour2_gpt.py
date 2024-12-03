# Import necessary modules
from shapely.geometry import Polygon

# Import SPICE wrapper functions with MATLAB-equivalent interfaces
from conversion_functions import *

# Import local modules
from mosaic_algorithms.auxiliar_functions.grid_functions.boustrophedon_gpt import boustrophedon
from mosaic_algorithms.auxiliar_functions.grid_functions.grid2D_gpt import grid2d
from mosaic_algorithms.auxiliar_functions.grid_functions.inst2topo_gpt import inst2topo
from mosaic_algorithms.auxiliar_functions.grid_functions.topo2inst_gpt import topo2inst
from mosaic_algorithms.auxiliar_functions.polygon_functions.closestSide2_gpt import closestSide2
from mosaic_algorithms.auxiliar_functions.polygon_functions.minimumWidthDirection_gpt import minimumWidthDirection
from mosaic_algorithms.auxiliar_functions.plot.groundtrack import groundtrack


# Import or define other required functions:
# topo2inst, groundtrack, closestSide2, minimumWidthDirection,
# grid2D, boustrophedon, inst2topo


def planSidewinderTour2(target, roi, sc, inst, inittime, ovlapx, ovlapy, angle):
    """
    Plans a sidewinder tour for satellite observation activities.

    Parameters:
        target (str): The target body or area.
        roi (list): Region of interest as a list of (x, y) coordinates.
        sc (str): Spacecraft identifier.
        inst (str): Instrument identifier.
        inittime (float): Initial time of observation.
        ovlapx (float): Overlap in the x-direction.
        ovlapy (float): Overlap in the y-direction.
        angle (float): Angle for the coverage path.

    Returns:
        grid (list): Grid points in instrument coordinates.
        origin (list): Origin point in instrument coordinates.
        itour (list): Boustrophedon tour indices.
        grid_topo (list): Grid points in topographical coordinates.
        tour (list): Tour in topographical coordinates.
        dirx (list): Direction vector in x.
        diry (list): Direction vector in y.
        dir1 (list): First direction vector for the coverage path.
        dir2 (list): Second direction vector for the coverage path.
    """

    # Pre-allocate variables
    origin = [0, 0]  # Origin in instrument coordinates

    # Point camera at ROI's centroid
    # MATLAB equivalent: [cx, cy] = centroid(polyshape(roi(:, 1), roi(:, 2)));
    # Conversion notes:
    # - In MATLAB, 'polyshape' creates a polygon, and 'centroid' computes its centroid.
    # - In Python, we use 'shapely.geometry.Polygon' for polygon creation and centroid computation.
    roi_polygon = Polygon(roi)
    centroid = roi_polygon.centroid
    cx, cy = centroid.x, centroid.y  # Coordinates of the centroid

    # Intersect ROI with focal plane
    # MATLAB equivalent: targetArea = topo2inst(roi, cx, cy, target, sc, inst, inittime);
    # Conversion notes:
    # - We assume 'topo2inst' is implemented with the same interface in Python.
    targetArea = topo2inst(roi, cx, cy, target, sc, inst, inittime)

    # Closest polygon side to the spacecraft's ground track position
    # This will determine the coverage path.
    # MATLAB equivalent:
    # gt1 = groundtrack(sc, inittime, target);
    # gt2 = groundtrack(sc, inittime + 500, target);
    gt1 = groundtrack(sc, inittime, target)
    gt2 = groundtrack(sc, inittime + 500, target)

    # Transform ground track positions to instrument coordinates
    # MATLAB equivalent:
    # gt1 = topo2inst(gt1, cx, cy, target, sc, inst, inittime);
    # gt2 = topo2inst(gt2, cx, cy, target, sc, inst, inittime + 500);
    gt1_inst = topo2inst(gt1, cx, cy, target, sc, inst, inittime)
    gt2_inst = topo2inst(gt2, cx, cy, target, sc, inst, inittime + 500)

    # Determine the closest side of the polygon to the ground track positions
    # MATLAB equivalent: [dir1, dir2] = closestSide2(gt1, gt2, targetArea, angle);
    dir1, dir2 = closestSide2(gt1=gt1_inst, gt2=gt2_inst, targetArea=targetArea, angle=angle)

    # Build reference tile
    # MATLAB equivalent:
    # [~, ~, ~, bounds] = cspice_getfov(cspice_bodn2c(inst), 4);
    # Conversion notes:
    # - 'cspice_bodn2c' converts an instrument name to its NAIF ID code.
    # - 'cspice_getfov' retrieves the field-of-view parameters.
    # - Using mat2py wrapper functions for SPICE routines.
    instrument_id = mat2py_bodn2c(inst)
    # Retrieve FOV parameters; we only need 'bounds'
    _, _, _, bounds = mat2py_getfov(instrument_id, 4)

    # Compute minimum width and height of the FOV
    # MATLAB equivalent: [~, width, height, ~] = minimumWidthDirection(bounds(1, :), bounds(2, :));
    # Conversion notes:
    # - Ensure that 'bounds' indexing matches between MATLAB and Python.
    # - Implement 'minimumWidthDirection' with the same logic as in MATLAB.
    _, width, height, _ = minimumWidthDirection(x=bounds[0], y=bounds[1])

    # Define focal plane reference parameters
    fpref = {
        'width': width,
        'height': height,
        'angle': angle
    }

    # Focal plane grid discretization
    # MATLAB equivalent: [grid, dirx, diry] = grid2D(fpref, ovlapx, ovlapy, origin, targetArea);
    # Conversion notes:
    # - Ensure 'grid2D' is implemented to accept parameters in the same way as in MATLAB.
    grid, dirx, diry = grid2d(fpref, ovlapx, ovlapy, origin, targetArea, fpThreshold=0.2)

    # Boustrophedon decomposition
    # MATLAB equivalent: itour = boustrophedon(grid, dir1, dir2);
    itour = boustrophedon(grid, dir1, dir2)

    # Convert grid and tour to topographical coordinates
    # MATLAB equivalent:
    # grid_topo = inst2topo(grid, cx, cy, target, sc, inst, inittime);
    # tour = inst2topo(itour, cx, cy, target, sc, inst, inittime);
    grid_topo = inst2topo(grid, cx, cy, target, sc, inst, inittime)
    tour = inst2topo(itour, cx, cy, target, sc, inst, inittime)

    # Remove empty entries from 'tour'
    # MATLAB equivalent:
    # indel = [];
    # for i=1:numel(tour)
    #     if isempty(tour{i}), indel = [indel i]; end
    # end
    # tour(indel) = [];
    # Conversion notes:
    # - In Python, we iterate over 'tour' and remove empty entries.
    tour = [t for t in tour if t]

    return grid, origin, itour, grid_topo, tour, dirx, diry, dir1, dir2
