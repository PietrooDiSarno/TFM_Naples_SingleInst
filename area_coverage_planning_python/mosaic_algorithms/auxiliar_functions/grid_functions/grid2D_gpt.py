def grid2d(fpref, olapx, olapy, gamma, target_area, fpThreshold):
    """
    Grid discretization (using flood-fill algorithm) of a region of interest
    given a reference footprint (unit measure to create the allocatable cells).

    Usage:
        matrix_grid, dirx, diry = grid2d(fpref, olapx, olapy, gamma, target_area)

    Inputs:
        - fpref:        Dictionary containing the parameters that define the footprint.
                        The following keys are needed:
            - width:    Footprint size in the x direction (longitude), in degrees.
            - height:   Footprint size in the y direction (latitude), in degrees.
            - angle:    Footprint orientation angle with respect to the meridian-equator axes, in degrees.
        - olapx:        Grid footprint overlap in the x direction, in percentage.
        - olapy:        Grid footprint overlap in the y direction, in percentage.
        - gamma:        Grid origin point (seed), a 2D point [x, y].
        - target_area:  N x 2 array containing the vertices of the ROI polygon.
                        target_area[:, 0] corresponds to x values (longitude),
                        target_area[:, 1] corresponds to y values (latitude).

    Outputs:
        - matrix_grid:  2D list (matrix) containing the grid discretization of the region-of-interest (ROI).
                        Each point is defined by the instrument boresight projection onto the body surface,
                        in latitudinal coordinates [lon, lat], in degrees.
        - dirx:         Direction vector x after rotation.
        - diry:         Direction vector y after rotation.

    The matrix sorts the discretized points (flood-fill) by latitude and longitude according to the following structure:

                    longitude
                    (-) --------> (+)
        latitude (+) [a11]  [a12] ⋯
                     [a21]
                      ⋮
                      ∨
                    (-)
    """
    import numpy as np
    from shapely.geometry import Polygon
    from scipy.spatial import ConvexHull
    # Ensure that required libraries are imported
    from mosaic_algorithms.auxiliar_functions.grid_functions.flood_fill_algorithm_gpt import flood_fill_algorithm

    # Pre-allocate variables
    matrix_grid = []

    # Get the footprint angle in radians
    angle = np.deg2rad(-fpref['angle'])

    # Create the rotation matrix
    rotmat = np.array([[np.cos(angle), -np.sin(angle)],
                       [np.sin(angle), np.cos(angle)]])

    # MatrixGrid directions x and y
    dirx = rotmat[0, :]
    diry = rotmat[1, :]

    # Calculate the centroid of the target area
    target_polygon = Polygon(target_area)
    cx, cy = target_polygon.centroid.coords[0]

    # Rotate the region-of-interest to align with the footprint
    oriented_area = np.zeros_like(target_area)
    for j in range(len(target_area)):
        delta = target_area[j, :] - np.array([cx, cy])
        rotated_delta = rotmat.dot(delta)
        oriented_area[j, :] = np.array([cx, cy]) + rotated_delta

    # Rotate the seed point gamma
    gamma = np.array([cx, cy]) + rotmat.dot(np.array(gamma) - np.array([cx, cy]))

    # If the area is divided into smaller regions, get the convex polygon that encloses all of them
    # Remove any NaN values from oriented_area
    aux = oriented_area[~np.isnan(oriented_area).any(axis=1)]

    # Compute the convex hull of the oriented area
    hull = ConvexHull(aux)
    peri_area = aux[hull.vertices]

    # Flood-fill algorithm to get the grid points of the oriented ROI
    grid_points, _ = flood_fill_algorithm(
        fpref['width'], fpref['height'], olapx, olapy,
        gamma, oriented_area, peri_area.tolist(),
        [], [], '4fill', fpThreshold)

    if grid_points:
        # Convert grid_points to NumPy array for processing
        grid_points = np.array(grid_points)

        # Sort grid points by latitude in descending order
        sorted_grid = grid_points[np.argsort(-grid_points[:, 1])]

        # Get unique latitude values
        unique_lat = np.unique(sorted_grid[:, 1])

        # Remove "similar" latitude values (differences less than 1e-5)
        if len(unique_lat) > 1:
            ind = np.abs(np.diff(unique_lat)) < 1e-5
            unique_lat = unique_lat[np.append(~ind, True)]  # Append True to match the array length

        # Get unique longitude values
        unique_lon = np.unique(sorted_grid[:, 0])

        # Remove "similar" longitude values (differences less than 1e-5)
        if len(unique_lon) > 1:
            ind = np.abs(np.diff(unique_lon)) < 1e-5
            unique_lon = unique_lon[np.append(~ind, True)]  # Append True to match the array length

        # Initialize the matrix grid as a list of lists
        matrix_grid = [[None for _ in range(len(unique_lon))] for _ in range(len(unique_lat))]

        # Populate the matrix grid
        for i in range(len(unique_lat)):
            # Latitude from highest to lowest
            lat = unique_lat[len(unique_lat) - 1 - i]
            # Find indices where latitude matches
            indlat = np.abs(sorted_grid[:, 1] - lat) < 1e-5
            # Get the corresponding longitudes
            mrow = sorted_grid[indlat, 0]
            mrow = np.sort(mrow)

            for lon in mrow:
                # Find the index of the longitude in unique_lon
                indlon = np.abs(unique_lon - lon) < 1e-5
                index_lon_array = np.where(indlon)[0]
                if index_lon_array.size > 0:
                    index_lon = index_lon_array[0]
                else:
                    continue  # Longitude not found in unique_lon

                # Rotate the point back to the original ROI orientation
                delta = np.array([lon, lat]) - np.array([cx, cy])
                rotated_point = np.array([cx, cy]) + rotmat.T.dot(delta)

                # Assign the rotated point to the matrix grid
                matrix_grid[i][index_lon] = rotated_point.tolist()

    return matrix_grid, dirx.tolist(), diry.tolist()
