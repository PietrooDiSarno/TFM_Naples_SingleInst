"""
Test Script for emissionang Function

This script tests the `emissionang` function by calculating the emission angle
between a target surface normal vector and the observer's position vector at
a given surface point and time(s). Results are visualized to verify the correctness
of the function.

"""

# Import external modules
import numpy as np
import matplotlib.pyplot as plt
from pySPICElib import kernelFetch

# Assume emissionang is defined in the same script or imported from another module
# from your_module import emissionang
from mosaic_algorithms.auxiliar_functions.observation_geometry.emissionang_gpt import emissionang


def main():
    """
    Main function to execute the emission angle test.

    - Defines a sample surface point and observation parameters.
    - Calls the emissionang function with various time values.
    - Visualizes the emission angle over time.
    """

    # Load SPICE kernels
    kf = kernelFetch(textFilesPath_='../')
    kf.ffFile(metaK='input/galileo/inputkernels.txt', forceDownload=False)

    # Define a sample surface point on the target (e.g., latitudinal coordinates for Mars)
    srfpoint = [0.0, 0.0]  # [lon, lat] in degrees, pointing to the equator at prime meridian

    # Define time values: a single time or an array of times to evaluate the emission angle
    t = np.linspace(0, 86400, 100)  # Times in TDB seconds past J2000 (1 day sampled at 100 points)

    # Specify the target and observer bodies
    target = 'EUROPA'  # Target body
    obs = 'GALILEO ORBITER'  # Spacecraft identifier

    # Calculate the emission angle
    angle = emissionang(srfpoint, t, target, obs)

    # Visualization
    visualize_results(t, angle)


def visualize_results(t, angle):
    """
    Visualizes the emission angle over time.

    Parameters:
    - t: Time values (array of floats).
    - angle: Emission angle at each time value (array of floats).
    """

    # Flatten inputs
    t = t.flatten()
    angle = angle.flatten()

    # Create a matplotlib figure and axis
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot the emission angle over time
    ax.plot(t, angle, label="Emission Angle", linewidth=2)
    ax.set_title("Emission Angle Over Time")
    ax.set_xlabel("Time (TDB seconds past J2000)")
    ax.set_ylabel("Emission Angle (degrees)")
    ax.legend()

    # Show grid
    ax.grid(True, linestyle='--', alpha=0.7)

    # Display the plot
    plt.show()


if __name__ == "__main__":
    main()
