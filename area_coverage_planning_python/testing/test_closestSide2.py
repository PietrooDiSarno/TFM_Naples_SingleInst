"""
Test Script for closestSide2 Function

This script tests the `closestSide2` function by analyzing a spacecraft's ground track in relation to a defined
region-of-interest (ROI). It computes the closest side of the ROI and the spacecraft’s movement direction.
The results are visualized to verify correctness.

"""

# Import external modules
import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, LineString, Point

# Import local modules
from mosaic_algorithms.auxiliar_functions.polygon_functions.closestSide2_gpt import closestSide2


# Main function
def main():
    """
    Main function to execute the closestSide2 test.

    - Defines a sample target area (ROI).
    - Sets initial and later spacecraft ground track positions.
    - Specifies footprint angle.
    - Calls closestSide2 function.
    - Visualizes the results.
    """
    # Define a sample ROI as a rectangle with vertices in latitude and longitude
    targetArea = np.array([
        [10, 10],
        [20, 10],
        [20, 20],
        [10, 20]
    ])  # Square ROI from (10, 10) to (20, 20) degrees

    # Initial and later positions of the spacecraft ground track
    gt1 = [15, 25]  # Initial position [longitude, latitude]
    gt2 = [18, 12]  # Later position [longitude, latitude]

    # Footprint angle (rotation of ROI)
    angle = 30  # degrees

    # Call the closestSide2 function
    dir1, dir2 = closestSide2(gt1, gt2, targetArea, angle)

    # Print the results
    print(f"Closest side of ROI to the spacecraft's initial position: {dir1}")
    print(f"Direction of spacecraft movement with respect to ROI: {dir2}")

    # Visualize the results
    visualize_results(targetArea, gt1, gt2, angle, dir1)


def visualize_results(targetArea, gt1, gt2, angle, dir1):
    """
    Visualizes the target area (ROI), the spacecraft ground track positions, and the closest side of the ROI.

    Parameters:
    - targetArea: Nx2 array of ROI vertices.
    - gt1: Initial spacecraft ground track position.
    - gt2: Later spacecraft ground track position.
    - angle: Footprint angle in degrees.
    - dir1: Closest side of the ROI to the spacecraft.
    """
    # Create a figure
    fig, ax = plt.subplots(figsize=(8, 6))

    # Rotate the ROI according to the footprint angle
    angle_rad = -np.radians(angle)
    rotmat = np.array([[np.cos(angle_rad), -np.sin(angle_rad)],
                       [np.sin(angle_rad), np.cos(angle_rad)]])

    # Compute centroid of the target area
    target_polygon = Polygon(targetArea)
    centroid = target_polygon.centroid
    cx, cy = centroid.x, centroid.y

    # Rotate the ROI around its centroid
    roi_rotated = np.array([cx, cy]) + (targetArea - np.array([cx, cy])) @ rotmat.T

    # Plot the rotated ROI
    roi_polygon = Polygon(roi_rotated)
    x_roi, y_roi = roi_polygon.exterior.xy
    ax.plot(x_roi, y_roi, 'b-', label='Rotated ROI')

    # Plot the spacecraft ground track positions
    ax.plot(gt1[0], gt1[1], 'ro', label='Initial Position (gt1)')
    ax.plot(gt2[0], gt2[1], 'go', label='Later Position (gt2)')

    # Highlight the closest side direction
    bbox = target_polygon.bounds
    if dir1 == 'north':
        ax.hlines(bbox[3], bbox[0], bbox[2], colors='orange', linestyle='--', label='Closest Side (North)')
    elif dir1 == 'south':
        ax.hlines(bbox[1], bbox[0], bbox[2], colors='orange', linestyle='--', label='Closest Side (South)')
    elif dir1 == 'east':
        ax.vlines(bbox[2], bbox[1], bbox[3], colors='orange', linestyle='--', label='Closest Side (East)')
    elif dir1 == 'west':
        ax.vlines(bbox[0], bbox[1], bbox[3], colors='orange', linestyle='--', label='Closest Side (West)')

    # Add plot details
    ax.set_title('Spacecraft Position and Closest Side to ROI')
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.legend()
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.set_aspect('equal')

    # Show plot
    plt.show()


if __name__ == "__main__":
    main()
