"""
Test Script for closestSide Function

This script tests the `closestSide` function, which calculates the spacecraft's ground track position relative
to the edges of a target region of interest (ROI) on a planetary body's surface.
It identifies the closest side of the ROI to the spacecraft's sub-observer point and the movement direction.
"""

import numpy as np
from pySPICElib import kernelFetch
from shapely.geometry import Polygon
import matplotlib.pyplot as plt
from spiceypy import utc2et

# Load local modules
from mosaic_algorithms.auxiliar_functions.polygon_functions.closestSide_gpt import closestSide


def main():
    """
    Main function to execute the closestSide test.

    - Loads SPICE kernels.
    - Defines a target area (ROI).
    - Sets observation parameters.
    - Calls the closestSide function.
    - Visualizes the ROI, the spacecraft's sub-observer position, and the direction results.
    """
    # Load SPICE kernels
    # Example kernel paths (modify these paths to point to your SPICE kernels)
    # furnsh('naif0012.tls')  # Leapseconds kernel
    # furnsh('de430.bsp')  # Planetary ephemeris
    # furnsh('pck00010.tpc')  # Planetary constants kernel
    # furnsh('spacecraft.bsp')  # Spacecraft SPK kernel
    # furnsh('spacecraft_frame.tf')  # Spacecraft frames kernel

    # Load SPICE kernels
    kf = kernelFetch()
    kf.ffFile(metaK='input/galileo/inputkernels.txt', forceDownload=False)

    # Define a sample target area (ROI) as a polygon in lat/lon
    targetArea = np.array([
        [10, -5],
        [20, -5],
        [20, 5],
        [10, 5],
        [10, -5]
    ])  # Rectangle from 10° to 20° longitude and -5° to 5° latitude

    # Observation parameters
    t = utc2et('2023-01-01T12:00:00')  # Observation time
    target = 'EUROPA'  # Target body
    sc = 'GALILEO ORBITER'  # Spacecraft identifier
    angle = 45  # Rotation angle in degrees

    # Call the closestSide function
    dir1, dir2 = closestSide(target, sc, t, targetArea, angle)

    # Display the results
    print(f"Closest side of ROI to the spacecraft's ground track: {dir1}")
    print(f"Spacecraft movement direction: {dir2}")

    # Visualize the results
    visualize_results(targetArea, dir1, dir2)


def visualize_results(targetArea, dir1, dir2):
    """
    Visualizes the ROI and labels the closest side and movement direction.

    Parameters:
    - targetArea: N x 2 array of ROI vertices in [longitude, latitude] degrees.
    - dir1: Closest side of the ROI to the spacecraft's ground track.
    - dir2: Spacecraft movement direction.
    """
    # Create a figure and axis for plotting
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the ROI
    roi_polygon = Polygon(targetArea)
    x_roi, y_roi = roi_polygon.exterior.xy
    ax.plot(x_roi, y_roi, 'k-', label='ROI')
    ax.fill(x_roi, y_roi, alpha=0.3, color='blue')

    # Label the closest side and movement direction
    ax.text(float(np.mean(x_roi)), float(np.mean(y_roi)), f"Closest Side: {dir1}\nMovement: {dir2}",
            ha='center', va='center', fontsize=12, bbox=dict(facecolor='white', alpha=0.8))

    # Plot settings
    ax.set_title('Closest Side of ROI Relative to Spacecraft Ground Track')
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_xlim(min(x_roi) - 5, max(x_roi) + 5)
    ax.set_ylim(min(y_roi) - 5, max(y_roi) + 5)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()

    plt.show()


if __name__ == "__main__":
    main()
