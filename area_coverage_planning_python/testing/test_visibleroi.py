"""
Test Script for visibleroi Function

This script tests the `visibleroi` function by computing the visible portion of a specified region of interest (ROI)
on a planetary body's surface from a given observer's viewpoint at a specific time. It visualizes both the ROI and
the visible portion to verify the correctness of the algorithm.

"""

import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from spiceypy import str2et
from pySPICElib.kernelFetch import kernelFetch

# Import modules
from mosaic_algorithms.auxiliar_functions.polygon_functions.visibleroi import visibleroi


# Main function
def main():
    """
    Main function to execute the visibleroi test.

    - Loads SPICE kernels.
    - Defines a target area (ROI).
    - Sets observation parameters.
    - Calls the visibleroi function.
    - Visualizes the results.
    """

    # Load SPICE kernels
    kf = kernelFetch(textFilesPath_='../')
    kf.ffFile(metaK='input/galileo/inputkernels.txt', forceDownload=False)

    # Define the target area (ROI) as a polygon in lat/lon
    roi = np.array([
        [0, -10],
        [30, -10],
        [30, 10],
        [0, 10],
        [0, -10]
    ])  # Simple rectangle from 0° to 30° longitude and -10° to 10° latitude

    # Observation parameters
    et = str2et('1998 MAY 30 00:00:00.000 TDB')  # Observation time
    target = 'EUROPA'  # Target body
    obs = 'GALILEO ORBITER'  # Observer

    # Call the visibleroi function
    vroi, inter = visibleroi(roi, et, target, obs)

    # Visualize the results
    visualize_results(roi, vroi, target, obs)


def visualize_results(roi, vroi, target, obs):
    """
    Visualizes the ROI and the visible portion from the observer's viewpoint.

    Parameters:
    - roi: N x 2 array of ROI vertices.
    - vroi: N x 2 array of visible ROI vertices.
    - target: Target body name (string).
    - obs: Observer name (string).
    """
    fig, ax = plt.subplots(figsize=(8, 6))

    # Plot the ROI
    roi_polygon = Polygon(roi)
    x_roi, y_roi = roi_polygon.exterior.xy
    ax.plot(x_roi, y_roi, 'k-', label='ROI')

    # Plot the visible ROI
    if vroi.size > 0:
        if vroi.ndim == 1:
            vroi = [vroi]
        else:
            vroi = [vroi]
        for v in vroi:
            ax.fill(v[:, 0], v[:, 1], facecolor='green', alpha=0.5, label='Visible ROI')

    # Configure the plot
    ax.set_title(f'Visible ROI on {target} from {obs}')
    ax.set_xlabel('Longitude (degrees)')
    ax.set_ylabel('Latitude (degrees)')
    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.grid(True, linestyle='--', alpha=0.7)
    ax.legend()

    plt.show()


if __name__ == "__main__":
    main()
