"""
Test Script for inst2topo and topo2inst Functions

This script tests the `inst2topo` and `topo2inst` functions by simulating a spacecraft instrument pointing scenario.
It converts a grid of points from instrument frame coordinates to topographical coordinates and back,
and visualizes the results to verify the correctness of the transformations.
"""

import numpy as np
import matplotlib.pyplot as plt
import spiceypy as spice
from pySPICElib.kernelFetch import kernelFetch

# Assume inst2topo and topo2inst are defined in the same script or imported from another module
# from your_module import inst2topo, topo2inst
from mosaic_algorithms.auxiliar_functions.grid_functions.inst2topo import inst2topo
from mosaic_algorithms.auxiliar_functions.grid_functions.topo2inst import topo2inst


def main():
    """
    Main function to execute the test of inst2topo and topo2inst.

    - Defines a grid in instrument frame coordinates.
    - Sets the spacecraft and instrument parameters.
    - Converts the grid to topographical coordinates using inst2topo.
    - Converts the topographical coordinates back to instrument frame using topo2inst.
    - Visualizes the original and reconstructed grids to verify the transformations.
    """

    # Load SPICE kernels
    kf = kernelFetch(textFilesPath_='../')
    kf.ffFile(metaK='input/galileo/inputkernels.txt', forceDownload=False)

    # Define instrument frame grid points (e.g., a simple 2D grid)
    grid_size = 5
    grid_range = np.linspace(-0.0001, 0.0001, grid_size)  # Small angles in radians
    grid_x, grid_y = np.meshgrid(grid_range, grid_range)
    grid = []
    for i in range(grid_size):
        grid_row = []
        for j in range(grid_size):
            grid_row.append([grid_x[i, j], grid_y[i, j]])
        grid.append(grid_row)

    # Spacecraft and instrument parameters
    lon = 0.0  # Instrument pointing longitude in degrees
    lat = 0.0  # Instrument pointing latitude in degrees
    target = 'EUROPA'  # Target body
    sc = 'GALILEO ORBITER'  # Spacecraft identifier
    inst = 'GLL_SSI'  # Instrument identifier
    et = spice.str2et('1998 MAY 30 00:00:00.000 TDB')  # Ephemeris time

    # Call inst2topo to convert instrument frame grid to topographical coordinates
    grid_topo = inst2topo(grid, lon, lat, target, sc, inst, et)

    # Call topo2inst to convert topographical coordinates back to instrument frame
    grid_reconstructed = topo2inst(grid_topo, lon, lat, target, sc, inst, et)

    # Visualize the original and reconstructed grids
    visualize_results(grid, grid_reconstructed)


def visualize_results(grid_original, grid_reconstructed):
    """
    Visualizes the original and reconstructed grids in instrument frame coordinates.

    Parameters:
    -----------
    grid_original : list of lists
        Original grid points in instrument frame coordinates.
    grid_reconstructed : list of lists
        Reconstructed grid points after transformations.
    """

    # Extract x and y coordinates from the grids
    x_original = np.array([[point[0] for point in row] for row in grid_original])
    y_original = np.array([[point[1] for point in row] for row in grid_original])

    x_reconstructed = np.array(
        [[point[0] if point is not None else np.nan for point in row] for row in grid_reconstructed])
    y_reconstructed = np.array(
        [[point[1] if point is not None else np.nan for point in row] for row in grid_reconstructed])

    # Plot the original grid
    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    plt.title('Original Instrument Frame Grid')
    plt.xlabel('X (radians)')
    plt.ylabel('Y (radians)')
    plt.scatter(x_original.flatten(), y_original.flatten(), color='blue', label='Original Grid')
    plt.legend()
    plt.grid(True)

    # Plot the reconstructed grid
    plt.subplot(1, 2, 2)
    plt.title('Reconstructed Instrument Frame Grid')
    plt.xlabel('X (radians)')
    plt.ylabel('Y (radians)')
    plt.scatter(x_reconstructed.flatten(), y_reconstructed.flatten(), color='red', label='Reconstructed Grid')
    plt.legend()
    plt.grid(True)

    # Display the plots
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
