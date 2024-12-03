"""
Test Script for flood_fill_algorithm Function

This script tests the `flood_fill_algorithm` function by applying it to a simple polygonal region.
It visualizes the flood-filled grid points over the target area to verify the correctness of the algorithm.
"""

# Import external packages
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon as MplPolygon

# Import local packages
from mosaic_algorithms.auxiliar_functions.grid_functions.floodFillAlgorithm import floodFillAlgorithm


# Define the main function
def main():
    """
    Main function to execute the flood fill test.

    - Defines a target area (polygon).
    - Sets the grid parameters.
    - Calls the flood_fill_algorithm function.
    - Visualizes the results.
    """

    # Define the target area as a simple polygon (e.g., a rectangle with a triangular cutout)
    target_area = np.array([
        [0, 0],
        [10, 0],
        [10, 5],
        [7, 5],
        [5, 8],
        [3, 5],
        [0, 5],
        [0, 0]
    ])

    # Initial perimeter area is the same as the target area
    perimeter_area = target_area.copy()

    # Grid parameters
    w = 1.0  # Width of each grid cell
    h = 1.0  # Height of each grid cell
    olapx = 0.0  # Overlap in x-direction (%)
    olapy = 0.0  # Overlap in y-direction (%)
    gamma = [5.0, 2.5]  # Starting point (seed) inside the target area
    method = '4fill'  # Flood fill method ('4fill' or '8fill')
    fpThreshold = 0.2  # Dismissal threshold

    # Initialize grid points and visited points lists
    grid_points = []
    v_points = []

    # Call the flood fill algorithm
    grid_points, v_points = floodFillAlgorithm(
        w, h, olapx, olapy, gamma,
        target_area, perimeter_area, grid_points, v_points, method
    )

    # Convert grid points to NumPy array for plotting
    grid_points = np.array(grid_points)

    # Visualization
    visualize_results(target_area, grid_points, w, h)


def visualize_results(target_area, grid_points, w, h):
    """
    Visualizes the target area and the flood-filled grid points.

    Parameters:
    - target_area: N x 2 array of the target area vertices.
    - grid_points: N x 2 array of grid point centers.
    - w: Width of each grid cell.
    - h: Height of each grid cell.
    """

    # Create a matplotlib figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the target area
    target_polygon = MplPolygon(target_area, closed=True, fill=None, edgecolor='black', linewidth=2,
                                label='Target Area')
    ax.add_patch(target_polygon)

    # Plot each grid cell
    for gp in grid_points:
        # Define the rectangle corners
        rect = MplPolygon([
            [gp[0] - w / 2, gp[1] - h / 2],
            [gp[0] - w / 2, gp[1] + h / 2],
            [gp[0] + w / 2, gp[1] + h / 2],
            [gp[0] + w / 2, gp[1] - h / 2]
        ], closed=True, facecolor='blue', edgecolor='blue', alpha=0.5)
        ax.add_patch(rect)

    # Set plot limits
    ax.set_xlim(min(target_area[:, 0]) - w, max(target_area[:, 0]) + w)
    ax.set_ylim(min(target_area[:, 1]) - h, max(target_area[:, 1]) + h)
    ax.set_aspect('equal')

    # Add legend and titles
    ax.legend()
    ax.set_title('Flood Fill Algorithm Test')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')

    # Show grid
    ax.grid(True, linestyle='--', alpha=0.7)

    # Display the plot
    plt.show()


if __name__ == "__main__":
    main()
