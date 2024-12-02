"""
Test Script for trgobsvec Function

This script tests the `trgobsvec` function by calculating the observer vector and distance
from a specified target surface point to an observer at a given epoch.
The results are visualized to validate the spatial relationship.

"""

import numpy as np
import matplotlib.pyplot as plt
from pySPICElib.kernelFetch import kernelFetch


# Assume trgobsvec is defined in the same script or imported from another module
# from your_module import trgobsvec
from mosaic_algorithms.auxiliar_functions.plot.trgobsvec import trgobsvec


def main():
    """
    Main function to execute the trgobsvec test.

    - Defines a surface point and observer-target configuration.
    - Calls trgobsvec to compute the observer vector and distance.
    - Visualizes the observer vector and distance.
    """

    # Load SPICE kernels
    kf = kernelFetch(textFilesPath_='../')
    kf.ffFile(metaK='input/galileo/inputkernels.txt', forceDownload=False)

    # Define a sample surface point on the target body
    srfpoint = [10.0, -45.0]  # Latitudinal coordinates in [longitude, latitude] degrees
    srfpoint3d = [10.0, -45.0, 1560.0]  # For EUROPA, mean radius ~1560.8 km

    # Define the observer and target configuration
    t = 0.0  # Time in TDB seconds past J2000 epoch
    target = 'EUROPA'  # Target body
    obs = 'GALILEO ORBITER'  # Spacecraft identifier
    frame = None  # Reference frame (if None: use the target frame)

    # Call the trgobsvec function
    obsvec, dist = trgobsvec(srfpoint, t, target, obs, frame)

    # Visualize results
    visualize_results(srfpoint3d, obsvec, dist)


def visualize_results(srfpoint, obsvec, dist, coord_type='latitudinal'):
    """
    Visualizes the target surface point and the observer vector in 3D.

    Parameters:
    - srfpoint: Target surface point. If `coord_type` is 'latitudinal', it should be a tuple/list
                (longitude in degrees, latitude in degrees, radius in km). If 'cartesian', (x, y, z).
    - obsvec: Observer position vector as seen from the target surface point (dx, dy, dz) in kilometers.
    - dist: Distance between the observer and the target surface point in kilometers.
    - coord_type: Coordinate system of `srfpoint`. Either 'cartesian' or 'latitudinal'.
    """

    # Convert latitudinal to Cartesian if necessary
    if coord_type == 'latitudinal':
        if len(srfpoint) != 3:
            raise ValueError("For 'latitudinal' coord_type, srfpoint must have three elements (longitude, latitude, radius).")
        lon_deg, lat_deg, radius = srfpoint
        lon_rad = np.radians(lon_deg)
        lat_rad = np.radians(lat_deg)
        x0 = radius * np.cos(lat_rad) * np.cos(lon_rad)
        y0 = radius * np.cos(lat_rad) * np.sin(lon_rad)
        z0 = radius * np.sin(lat_rad)
    elif coord_type == 'cartesian':
        if len(srfpoint) != 3:
            raise ValueError("For 'cartesian' coord_type, srfpoint must have three elements (x, y, z).")
        x0, y0, z0 = srfpoint
    else:
        raise ValueError("coord_type must be either 'cartesian' or 'latitudinal'.")

    dx, dy, dz = obsvec
    observer_pos = (x0 + dx, y0 + dy, z0 + dz)

    # Ensure 'dist' is a scalar float
    if isinstance(dist, np.ndarray):
        if dist.size == 1:
            dist_scalar = dist.item()
        else:
            raise ValueError(f"Expected 'dist' to be a single value, but got an array with shape {dist.shape}.")
    else:
        dist_scalar = dist

    # Verify that the provided distance matches the observer vector's magnitude
    computed_dist = np.linalg.norm(obsvec)
    if not np.isclose(computed_dist, dist_scalar, atol=1e-2):
        raise ValueError(f"The provided distance {dist_scalar} km "
                         f"does not match the magnitude of obsvec {computed_dist:.2f} km.")

    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.view_init(elev=20, azim=30)

    # Plot the target surface point
    ax.scatter(x0, y0, z0, color='red', s=100, label="Target Surface Point")

    # Plot the observer vector from the surface point
    ax.quiver(x0, y0, z0, dx, dy, dz, color='blue', arrow_length_ratio=0.1, label="Observer Vector")

    # Mark the observer's position
    ax.scatter(observer_pos[0], observer_pos[1], observer_pos[2], color='blue', s=50, label="Observer Position")

    # Display distance as an annotation at the observer's position
    # ax.text(observer_pos[0], observer_pos[1], observer_pos[2], f"Dist: {dist_scalar:.2f} km", color='black')

    # Determine plot limits to ensure both points are visible with some padding
    all_x = [x0, observer_pos[0]]
    all_y = [y0, observer_pos[1]]
    all_z = [z0, observer_pos[2]]
    max_range = max(
        abs(max(all_x) - min(all_x)),
        abs(max(all_y) - min(all_y)),
        abs(max(all_z) - min(all_z))
    ) * 1.5  # Add 50% padding

    mid_x = (max(all_x) + min(all_x)) / 2
    mid_y = (max(all_y) + min(all_y)) / 2
    mid_z = (max(all_z) + min(all_z)) / 2

    ax.set_xlim(mid_x - max_range / 2, mid_x + max_range / 2)
    ax.set_ylim(mid_y - max_range / 2, mid_y + max_range / 2)
    ax.set_zlim(mid_z - max_range / 2, mid_z + max_range / 2)

    # Set labels and title
    ax.set_xlabel("X [km]")
    ax.set_ylabel("Y [km]")
    ax.set_zlabel("Z [km]")
    ax.set_title("Observer Vector and Distance from Target Surface Point")

    # Add legend
    ax.legend()

    # Optional: Enhance the view angle for better visualization
    ax.view_init(elev=20, azim=30)

    plt.show()




if __name__ == "__main__":
    main()
