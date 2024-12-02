"""
Test Script for planSidewinderTour2 Function

This script tests the `planSidewinderTour2` function by planning a sidewinder
tour for satellite observation activities over a specified region of interest (ROI).
It provides a toy example with mocked dependencies to facilitate verification and
validation of the function's behavior. The results, including grid points and
the planned tour, are visualized for intuitive understanding.

"""

# Import external modules
import matplotlib.pyplot as plt
import numpy as np
from pySPICElib import kernelFetch
from shapely.geometry import Polygon

# Helper functions (amsplit, grid2D, boustrophedon, etc.)
# Assuming these functions are defined in the same module or are imported
from mosaic_algorithms.sidewinder.planSidewinderTour import planSidewinderTour


# Define the main function to execute the test
def main():
    """
    Main function to execute the planSidewinderTour2 test.

    - Defines a target area (ROI).
    - Sets planning parameters.
    - Calls the planSidewinderTour2 function.
    - Visualizes the grid points and the planned tour.
    """

    # Load SPICE kernels
    kf = kernelFetch(textFilesPath_='../')
    kf.ffFile(metaK='input/galileo/inputkernels.txt', forceDownload=False)

    # Define the target area as a simple rectangular ROI
    roi =np.array([
        [0, 0],
        [20, 0],
        [20, 10],
        [0, 10],
        [0, 0]
    ])  # Rectangle from (0,0) to (20,10)

    # Define planning parameters
    target = 'EUROPA'  # Target body
    sc = 'GALILEO ORBITER'  # Spacecraft identifier
    inst = 'GLL_SSI'  # Instrument identifier
    inittime = 0.0  # Initial time of observation (arbitrary units)
    ovlapx = 0.2  # 20% overlap in x-direction
    ovlapy = 0.2  # 20% overlap in y-direction
    angle = 45.0  # Coverage path angle in degrees

    # Call the planSidewinderTour2 function
    grid, origin, itour, grid_topo, tour, dirx, diry, dir1, dir2 = planSidewinderTour(
        target=target,
        roi=roi,
        sc=sc,
        inst=inst,
        inittime=inittime,
        olapx=ovlapx,
        olapy=ovlapy
        #angle=angle
    )


    # Convert lists to NumPy arrays for easier manipulation
    grid = np.array(grid)
    grid_topo = np.array(grid_topo)
    tour = np.array(tour)

    # Visualization
    visualize_results(roi, grid, grid_topo, tour, itour, dir1, dir2)


def visualize_results(roi, grid, grid_topo, tour, itour, dir1, dir2):
    """
    Visualizes the ROI, grid points, and the planned sidewinder tour.

    Parameters:
    - roi (list): Original region of interest coordinates.
    - grid (np.ndarray): Grid points in instrument coordinates.
    - grid_topo (np.ndarray): Grid points in topographical coordinates.
    - tour (np.ndarray): Ordered list of grid points representing the tour.
    - itour (list): Indices of the tour in the grid.
    - dir1 (list): First direction vector for coverage path.
    - dir2 (list): Second direction vector for coverage path.
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    # Plot the ROI polygon
    roi_polygon = Polygon(roi)
    x_roi, y_roi = roi_polygon.exterior.xy
    ax.plot(x_roi, y_roi, 'k-', linewidth=2, label='ROI')

    # Plot the grid points
    ax.scatter(grid_topo[:, 0], grid_topo[:, 1], c='blue', label='Grid Points')

    # Plot the tour path
    if len(tour) > 0:
        ax.plot(tour[:, 0], tour[:, 1], 'r--', linewidth=1.5, label='Planned Tour')
        ax.scatter(tour[:, 0], tour[:, 1], c='red')

    # Plot direction vectors
    origin = [0, 0]
    scale = 5  # Scale for visualization
    ax.arrow(origin[0], origin[1], dir1[0] * scale, dir1[1] * scale,
             head_width=0.5, head_length=0.5, fc='green', ec='green', label='Direction 1')
    ax.arrow(origin[0], origin[1], dir2[0] * scale, dir2[1] * scale,
             head_width=0.5, head_length=0.5, fc='purple', ec='purple', label='Direction 2')

    # Configure plot aesthetics
    ax.set_title('Sidewinder Tour Planning Test')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.set_aspect('equal', 'box')
    ax.grid(True, linestyle='--', alpha=0.7)

    # Create a custom legend to handle arrows
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color='k', lw=2, label='ROI'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='blue', markersize=8, label='Grid Points'),
        Line2D([0], [0], color='r', linestyle='--', lw=2, label='Planned Tour'),
        Line2D([0], [0], marker=r'$\rightarrow$', color='w', markerfacecolor='green', markersize=10,
               label='Direction 1'),
        Line2D([0], [0], marker=r'$\rightarrow$', color='w', markerfacecolor='purple', markersize=10,
               label='Direction 2')
    ]
    ax.legend(handles=custom_lines, loc='upper right')

    # Show the plot
    plt.show()


if __name__ == "__main__":
    main()
