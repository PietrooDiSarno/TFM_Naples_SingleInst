"""
Test Script for boustrophedon Function

This script tests the `boustrophedon` function by applying it to a grid of points.
It visualizes the planned tour over the grid to verify the correctness of the boustrophedon path planning.
"""

import numpy as np
import matplotlib.pyplot as plt

# Assume boustrophedon is defined in the same script or imported from another module
# from your_module import boustrophedon
from mosaic_algorithms.auxiliar_functions.grid_functions.boustrophedon import boustrophedon


def main():
    """
    Main function to execute the boustrophedon path planning test.

    - Defines a grid with observation points.
    - Sets the sweeping directions.
    - Calls the boustrophedon function.
    - Visualizes the planned tour.
    """

    # Define grid dimensions
    num_rows = 5
    num_cols = 5
    grid_spacing = 1.0  # Distance between grid points

    # Generate grid points
    grid = [[None for _ in range(num_cols)] for _ in range(num_rows)]
    for i in range(num_rows):
        for j in range(num_cols):
            x = j * grid_spacing
            y = i * grid_spacing
            grid[i][j] = (x, y)

    # Introduce some empty cells (None) to simulate obstacles or missing points
    grid[2][2] = None  # Remove center point
    grid[0][4] = None  # Remove a corner point
    grid[4][0] = None  # Remove another corner point

    # Set sweeping directions
    dir1 = 'north'  # Primary direction ('north', 'south', 'east', 'west')
    dir2 = 'east'  # Secondary direction ('north', 'south', 'east', 'west')

    # Call the boustrophedon function
    tour = boustrophedon(grid, dir1, dir2)

    # Visualization
    visualize_tour(grid, tour, grid_spacing, num_rows, num_cols)


def visualize_tour(grid, tour, grid_spacing, num_rows, num_cols):
    """
    Visualizes the grid and the planned boustrophedon tour.

    Parameters:
    - grid: 2D list of grid points.
    - tour: Ordered list of points representing the tour.
    - grid_spacing: Distance between grid points.
    """

    # Extract all grid points for plotting
    all_points = [point for row in grid for point in row if point is not None]
    all_x = [p[0] for p in all_points]
    all_y = [p[1] for p in all_points]

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot grid points
    ax.scatter(all_x, all_y, c='lightgray', s=100, label='Grid Points')

    # Plot the tour path
    tour_x = [point[0] for point in tour]
    tour_y = [point[1] for point in tour]
    ax.plot(tour_x, tour_y, c='blue', linewidth=2, label='Planned Tour')
    ax.scatter(tour_x, tour_y, c='red', s=50, zorder=5)

    # Annotate points with their sequence number in the tour
    for idx, point in enumerate(tour):
        ax.text(point[0], point[1] + 0.1, str(idx + 1), ha='center', va='bottom', fontsize=9, color='black')

    # Set plot limits
    ax.set_xlim(-grid_spacing, num_cols * grid_spacing)
    ax.set_ylim(-grid_spacing, num_rows * grid_spacing)
    ax.set_aspect('equal')

    # Add grid lines
    ax.set_xticks(np.arange(0, num_cols * grid_spacing, grid_spacing))
    ax.set_yticks(np.arange(0, num_rows * grid_spacing, grid_spacing))
    ax.grid(True, linestyle='--', alpha=0.5)

    # Add titles and labels
    ax.set_title('Boustrophedon Path Planning')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.legend()

    # Display the plot
    plt.show()


if __name__ == "__main__":
    main()
