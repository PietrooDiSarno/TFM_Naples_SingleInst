"""
Test Script for instpointing Function

This script tests the `instpointing` function by setting up a mock scenario involving a spacecraft,
a target body, and an instrument. It visualizes the instrument's field of view (FOV) and boresight
in a simple 3D plot to verify the orientation and visibility of the target point.

"""

# Import external modules
import matplotlib.pyplot as plt
from pySPICElib import kernelFetch

# Assume instpointing is defined in the same script or imported from another module
# from your_module import instpointing
from mosaic_algorithms.auxiliar_functions.spacecraft_operation.instpointing import instpointing


def main():
    """
    Main function to execute the instpointing test.

    - Sets up parameters for a sample instrument, target, and spacecraft.
    - Calls the instpointing function with mock data.
    - Visualizes the FOV bounds, boresight, and target visibility.
    """

    # Load SPICE kernels
    kf = kernelFetch(textFilesPath_='../')
    kf.ffFile(metaK='input/galileo/inputkernels.txt', forceDownload=False)

    # Instrument, target, spacecraft parameters
    inst = 'GLL_SSI'  # Instrument identifier
    target = 'EUROPA'  # Target body
    sc = 'GALILEO ORBITER'  # Spacecraft identifier
    t = 1000  # Observation time in seconds

    # Define optional longitude and latitude (in degrees)
    lon, lat = 45, -30  # Mock values for testing

    # Call the instpointing function
    #fovbounds, boresight, rotmat, visible, lon_out, lat_out = instpointing(inst, target, sc, t, lon, lat)
    fovbounds, boresight, rotmat, visible = instpointing(inst, target, sc, t, lon, lat)
    lon_out = lon
    lat_out = lat
    # Display the results
    print("FOV Bounds:", fovbounds)
    print("Boresight:", boresight)
    print("Rotation Matrix:\n", rotmat)
    print("Is Target Visible:", visible)
    print("Longitude:", lon_out, "Latitude:", lat_out)

    # Visualization of the FOV and boresight
    visualize_instpointing(fovbounds, boresight, visible, lon_out, lat_out)


def visualize_instpointing(fovbounds, boresight, visible, lon_out, lat_out):
    """
    Visualizes the instrument's field of view (FOV) and boresight in a 3D plot.

    Parameters:
    - fovbounds: FOV bounds in the body-fixed reference frame.
    - boresight: Boresight vector in the body-fixed reference frame.
    - visible: Boolean indicating if the target point is visible.
    - lon_out: Longitude of the target point.
    - lat_out: Latitude of the target point.
    """

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Plot FOV bounds
    for i in range(fovbounds.shape[1]):
        ax.plot([0, fovbounds[0, i]], [0, fovbounds[1, i]], [0, fovbounds[2, i]], 'b-')

    # Plot boresight vector
    ax.quiver(0, 0, 0, boresight[0], boresight[1], boresight[2], color='red', length=1.5, normalize=True,
              label='Boresight')

    # Plot the target visibility status
    target_color = 'green' if visible else 'gray'
    ax.scatter([0], [0], [0], color=target_color, s=100, label=f"Target ({'Visible' if visible else 'Not Visible'})")

    # Add labels and legend
    ax.set_xlabel("X-axis")
    ax.set_ylabel("Y-axis")
    ax.set_zlabel("Z-axis")
    ax.legend()
    ax.set_title(f"Instrument FOV and Target Visibility\nLon: {lon_out}°, Lat: {lat_out}°")

    # Set aspect ratio and display
    ax.set_box_aspect([1, 1, 1])
    plt.show()


if __name__ == "__main__":
    main()
