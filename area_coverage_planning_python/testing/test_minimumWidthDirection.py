"""
Test Script for minimumWidthDirection Function

This script tests the `minimumWidthDirection` function by applying it to a sample polygon.
It visualizes the polygon, as well as the minimum width direction and height,
to verify the function's outputs in a visually intuitive manner.

"""

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from matplotlib.patches import Polygon as MplPolygon
from matplotlib.lines import Line2D

# Assume minimumWidthDirection is defined in the same script or imported from another module
# from your_module import minimumWidthDirection
from mosaic_algorithms.auxiliar_functions.polygon_functions.minimumWidthDirection import minimumWidthDirection


def main():
    """
    Main function to execute the minimum width direction test.

    - Defines a sample polygon.
    - Calls the minimumWidthDirection function to calculate the minimum width direction.
    - Visualizes the polygon with its minimum width and height direction.
    """

    # Define the polygon vertices as arrays of x and y coordinates
    x = [0, 10, 15, 8, 3, 0]
    y = [0, 0, 7, 12, 10, 5]

    # Call the function to calculate minimum width direction and other outputs
    thetamin, minwidth, height, axes = minimumWidthDirection(x, y)

    # Print results to verify calculations
    print(f"Minimum Width Angle: {thetamin} degrees")
    print(f"Minimum Width: {minwidth}")
    print(f"Height: {height}")
    print(f"Direction Vector (Axes): {axes}")

    # Visualization
    visualize_results(x, y, thetamin, minwidth, height, axes)


def visualize_results(x, y, thetamin, minwidth, height, axes):
    """
    Visualizes the polygon, minimum width direction, and height based on the function's outputs.

    Parameters:
    - x, y: Arrays of x and y coordinates of the polygon vertices.
    - thetamin: Angle (in degrees) at which the polygon's width is minimized.
    - minwidth: Minimum width of the polygon.
    - height: Height of the polygon, at the orientation specified by thetamin.
    - axes: Unit vector representing the direction of the minimum width.
    """

    # Create a matplotlib figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the polygon
    vertices = np.column_stack((x, y))
    polygon_shape = MplPolygon(vertices, closed=True, fill=None, edgecolor='black', linewidth=2, label='Polygon')
    ax.add_patch(polygon_shape)

    # Calculate and plot the minimum width and height lines
    centroid = Polygon(vertices).centroid
    cx, cy = centroid.x, centroid.y

    # Plot the centroid
    ax.plot(cx, cy, 'ro', label="Centroid")

    # Minimum width line endpoints
    direction_vector = np.array(axes) * minwidth / 2
    point1 = [cx - direction_vector[0], cy - direction_vector[1]]
    point2 = [cx + direction_vector[0], cy + direction_vector[1]]
    min_width_line = Line2D([point1[0], point2[0]], [point1[1], point2[1]], color='blue', linewidth=2,
                            label='Min Width')

    # Height line endpoints, perpendicular to min width
    height_vector = np.array([-axes[1], axes[0]]) * height / 2
    height_point1 = [cx - height_vector[0], cy - height_vector[1]]
    height_point2 = [cx + height_vector[0], cy + height_vector[1]]
    height_line = Line2D([height_point1[0], height_point2[0]], [height_point1[1], height_point2[1]], color='green',
                         linewidth=2, label='Height')

    # Add the lines to the plot
    ax.add_line(min_width_line)
    ax.add_line(height_line)

    # Set plot limits
    ax.set_xlim(min(x) - 2, max(x) + 2)
    ax.set_ylim(min(y) - 2, max(y) + 2)
    ax.set_aspect('equal')

    # Add legend, titles, and labels
    ax.legend()
    ax.set_title('Minimum Width Direction Test')
    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')

    # Show grid
    ax.grid(True, linestyle='--', alpha=0.7)

    # Display the plot
    plt.show()


if __name__ == "__main__":
    main()
