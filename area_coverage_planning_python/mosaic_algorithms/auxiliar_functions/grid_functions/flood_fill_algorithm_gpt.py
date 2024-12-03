def flood_fill_algorithm(w, h, olapx, olapy, gamma,
                         target_area, perimeter_area, grid_points, v_points, method,
                         fpThreshold):
    """
    Flood-fill recursive algorithm that discretizes the target area by
    "flooding" the region with 2D rectangular elements. The grid is
    determined by the input width, height, and overlaps in both directions.

    Usage:
        grid_points = flood_fill_algorithm(w, h, olapx, olapy,
                                           gamma, target_area, grid_points, method)

    Inputs:
        - w:            Horizontal resolution. Units are irrelevant as long as
                        they are consistent.
        - h:            Vertical resolution. Units are irrelevant as long as
                        they are consistent.
        - olapx:        Grid footprint overlap in the horizontal direction.
                        Units are in percentage of width.
        - olapy:        Grid footprint overlap in the vertical direction. Units
                        are in percentage of height.
        - gamma:        Grid origin point (seed), a 2D point [x, y].
        - target_area:  N x 2 array containing the vertices of the ROI polygon.
                        target_area[:, 0] corresponds to x values,
                        target_area[:, 1] corresponds to y values.
        - perimeter_area:
                        N x 2 array containing the vertices of the polygon that
                        encloses all the uncovered area. Initially,
                        perimeter_area = target_area, but it changes as observations
                        advance.
        - grid_points:  List of grid points representing the centers of the
                        rectangular elements used to fill the region.
                        When calling this function: grid_points = [].
        - v_points:     List of visited points to prevent gridlock.
                        When calling this function: v_points = [].
        - method:       String name of the method. '4fill' fills the ROI by
                        searching the cardinal directions. '8fill' considers
                        also the diagonal neighbors.
        - fpThreshold:  Threshold for XXX dismissal.

    Outputs:
        - grid_points:  Updated list containing the discretized grid points of the
                        region-of-interest.
        - v_points:     Updated list containing the visited points to prevent gridlock.

    Note:
        This function creates a convex polygon that encloses all the
        uncovered area (even when this is divided into portions) and tours the
        whole area. It is less computationally efficient than the classic
        flood-fill algorithm but helps prevent suboptimal fillings of the
        uncovered area (isolated points).
    """

    # Import necessary libraries
    import numpy as np
    from shapely.geometry import Polygon, Point

    # Variables pre-allocation
    inside = False

    # Convert overlaps from percentage to actual units (consistent with w and h)
    ovlapx = olapx * w / 100
    ovlapy = olapy * h / 100

    # Check if the cell has been previously visited
    for vp in v_points:
        if np.linalg.norm(np.array(vp) - np.array(gamma)) < 1e-5:
            return grid_points, v_points  # Early return if already visited

    # Otherwise, add the current gamma to visited points
    v_points.append(gamma)

    # Rectangular element definition (footprint)
    fpx = [gamma[0] - w / 2, gamma[0] - w / 2, gamma[0] + w / 2, gamma[0] + w / 2]
    fpy = [gamma[1] + h / 2, gamma[1] - h / 2, gamma[1] - h / 2, gamma[1] + h / 2]
    fpshape = Polygon(zip(fpx, fpy))

    # Subtract the allocated cell (footprint) from the perimeter_area
    peripshape = Polygon(perimeter_area)
    inter = peripshape.difference(fpshape)
    areaI = inter.area
    areaP = peripshape.area

    # Check if the footprint is larger than the region of interest
    if areaI == 0:
        grid_points.append(gamma)
        return grid_points, v_points

    # Check if the rectangle at gamma and size [w,h] is contained in the perimeter area
    if peripshape.intersects(Point(gamma)):
        inside = True
    else:
        if abs(areaI - areaP) / fpshape.area > 0.2:
            inside = True

    if inside:
        # Disregard cases where the footprint does not cover a certain minimum of the ROI
        targetpshape = Polygon(target_area)
        areaT = targetpshape.area
        inter = targetpshape.difference(fpshape)
        areaI = inter.area
        areaInter = areaT - areaI
        fpArea = fpshape.area

        # Compare area of footprint to intersection of ROI with footprint
        # If this is below fpThreshold, the footprint is dismissed.
        if areaInter / fpArea > fpThreshold:
            grid_points.append(gamma)
            # Visualization code omitted
        else:
            if not grid_points:
                return grid_points, v_points  # First iteration: footprint doesn't cover minimum required
            # Visualization code omitted

        # Check the neighbors recursively
        # West
        grid_points, v_points = flood_fill_algorithm(
            w, h, olapx, olapy, [gamma[0] - w + ovlapx, gamma[1]],
            target_area, perimeter_area, grid_points, v_points, method, fpThreshold)
        # South
        grid_points, v_points = flood_fill_algorithm(
            w, h, olapx, olapy, [gamma[0], gamma[1] - h + ovlapy],
            target_area, perimeter_area, grid_points, v_points, method, fpThreshold)
        # North
        grid_points, v_points = flood_fill_algorithm(
            w, h, olapx, olapy, [gamma[0], gamma[1] + h - ovlapy],
            target_area, perimeter_area, grid_points, v_points, method, fpThreshold)
        # East
        grid_points, v_points = flood_fill_algorithm(
            w, h, olapx, olapy, [gamma[0] + w - ovlapx, gamma[1]],
            target_area, perimeter_area, grid_points, v_points, method, fpThreshold)

        if method == '8fill':
            # Northwest
            grid_points, v_points = flood_fill_algorithm(
                w, h, olapx, olapy, [gamma[0] - w + ovlapx, gamma[1] + h - ovlapy],
                target_area, perimeter_area, grid_points, v_points, method, fpThreshold)
            # Southwest
            grid_points, v_points = flood_fill_algorithm(
                w, h, olapx, olapy, [gamma[0] - w + ovlapx, gamma[1] - h + ovlapy],
                target_area, perimeter_area, grid_points, v_points, method, fpThreshold)
            # Northeast
            grid_points, v_points = flood_fill_algorithm(
                w, h, olapx, olapy, [gamma[0] + w - ovlapx, gamma[1] + h - ovlapy],
                target_area, perimeter_area, grid_points, v_points, method, fpThreshold)
            # Southeast
            grid_points, v_points = flood_fill_algorithm(
                w, h, olapx, olapy, [gamma[0] + w - ovlapx, gamma[1] - h + ovlapy],
                target_area, perimeter_area, grid_points, v_points, method, fpThreshold)

    return grid_points, v_points
