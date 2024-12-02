"""
Test Script for grid2d Function

This script tests the `grid2d` function by applying it to a simple polygonal region.
It visualizes the discretized grid over the target area to verify the correctness of the algorithm.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

# Assume grid2d and flood_fill_algorithm are defined in the same script or imported from another module
# from your_module import grid2d, flood_fill_algorithm
from mosaic_algorithms.auxiliar_functions.grid_functions.grid2D_gpt import grid2d


def main():
    """
    Main function to execute the grid2d test.

    - Defines a target area (polygon).
    - Sets the footprint and grid parameters.
    - Calls the grid2d function.
    - Visualizes the results.
    """

    # Define the target area as a simple polygon (e.g., a rectangle)
    target_area = np.array([
        [0, 0],
        [10, 0],
        [10, 5],
        [5, 8],
        [0, 5],
        [0, 0]
    ])

    # Footprint parameters
    fpref = {
        'width': 1.0,  # Footprint width in degrees
        'height': 1.0,  # Footprint height in degrees
        'angle': 30.0  # Footprint orientation angle in degrees
    }

    # Grid parameters
    olapx = 0.0  # Overlap in x-direction (%)
    olapy = 0.0  # Overlap in y-direction (%)
    gamma = [5.0, 2.5]  # Starting point (seed) inside the target area
    fpThreshold = 0.2  # Dismissal threshold

    # Call the grid2d function
    matrix_grid, dirx, diry = grid2d(fpref, olapx, olapy, gamma, target_area,fpThreshold)

    # Visualization
    visualize_results(target_area, matrix_grid, fpref['width'], fpref['height'], fpref['angle'])


def visualize_results(target_area, matrix_grid, w, h, angle):
    """
    Visualizes the target area and the grid discretization.

    Parameters:
    - target_area: N x 2 array of the target area vertices.
    - matrix_grid: 2D list containing the grid points.
    - w: Width of each grid cell.
    - h: Height of each grid cell.
    - angle: Rotation angle of the grid cells.
    """

    # Create a matplotlib figure and axis
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the target area
    target_polygon = MplPolygon(target_area, closed=True, fill=None, edgecolor='black', linewidth=2,
                                label='Target Area')
    ax.add_patch(target_polygon)

    # Flatten the matrix_grid to a list of points
    grid_points = []
    for row in matrix_grid:
        for point in row:
            if point is not None:
                grid_points.append(point)

    # Plot each grid cell
    for gp in grid_points:
        # Calculate the rotated rectangle corners
        rect = get_rotated_rectangle(gp, w, h, angle)
        rect_patch = MplPolygon(rect, closed=True, facecolor='blue', edgecolor='blue', alpha=0.5)
        ax.add_patch(rect_patch)

    # Set plot limits
    margin = max(w, h)
    ax.set_xlim(min(target_area[:, 0]) - margin, max(target_area[:, 0]) + margin)
    ax.set_ylim(min(target_area[:, 1]) - margin, max(target_area[:, 1]) + margin)
    ax.set_aspect('equal')

    # Add legend and titles
    ax.legend()
    ax.set_title('Grid2D Algorithm Test')
    ax.set_xlabel('Longitude')
    ax.set_ylabel('Latitude')

    # Show grid
    ax.grid(True, linestyle='--', alpha=0.7)

    # Display the plot
    plt.show()


def get_rotated_rectangle(center, width, height, angle):
    """
    Calculates the corners of a rotated rectangle.

    Parameters:
    - center: [x, y] coordinates of the rectangle center.
    - width: Width of the rectangle.
    - height: Height of the rectangle.
    - angle: Rotation angle in degrees.

    Returns:
    - corners: List of [x, y] coordinates of the rectangle corners.
    """

    angle_rad = np.deg2rad(angle)
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)

    # Define rectangle corners relative to center
    dx = width / 2
    dy = height / 2
    corners = np.array([
        [-dx, -dy],
        [-dx, dy],
        [dx, dy],
        [dx, -dy]
    ])

    # Rotate corners
    rotation_matrix = np.array([[cos_a, -sin_a],
                                [sin_a, cos_a]])
    rotated_corners = np.dot(corners, rotation_matrix.T)

    # Translate corners to center
    rotated_corners += center

    return rotated_corners


if __name__ == "__main__":
    main()
