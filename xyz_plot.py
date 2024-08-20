from simplekml import Kml
import numpy as np
import os
from pyproj import Transformer

# File path (replace this with your actual path)
filename = '/home/mathjo/catkin_ws/src/GILO-EKF/output/m2dgr/street_05/odo.txt'

# Set up the UTM to WGS84 transformer (replace with your correct UTM zone)
transformer = Transformer.from_crs("epsg:32633", "epsg:4326")  # Replace EPSG:32633 with your UTM zone EPSG code

try:
    # Load the data from the file, ignoring comment lines
    data = np.loadtxt(filename, comments='#')

    # Extract XYZ columns (x, y, z in meters)
    x_coords = data[:, 1]
    y_coords = data[:, 2]
    z_coords = data[:, 3]

    # Convert Cartesian (x, y) to geographic coordinates (longitude, latitude)
    longitudes = []
    latitudes = []
    for x, y in zip(x_coords, y_coords):
        lon, lat = transformer.transform(x, y)
        longitudes.append(lon)
        latitudes.append(lat)

    # Create a KML object
    kml = Kml()

    # Add the trajectory to the KML file as a line
    linestring = kml.newlinestring(name="3D Trajectory")
    linestring.coords = [(lon, lat, z) for lon, lat, z in zip(longitudes, latitudes, z_coords)]
    linestring.altitudemode = "absolute"  # Use absolute altitude mode for true 3D positioning
    linestring.extrude = 0  # Disable extrusion (no vertical line to the ground)
    linestring.style.linestyle.color = "ff0000ff"  # Red color (in KML format)
    linestring.style.linestyle.width = 3  # Line width

    # KML file path
    kml_file_path = '/home/mathjo/catkin_ws/src/GILO-EKF/output/m2dgr/street_05/trajectory_xyz.kml'

    # Ensure the directory exists
    os.makedirs(os.path.dirname(kml_file_path), exist_ok=True)

    # Save the KML file
    kml.save(kml_file_path)
    print(f"KML file saved as '{kml_file_path}'.")

except FileNotFoundError:
    print(f"File '{filename}' not found.")
except Exception as e:
    print(f"An error occurred: {e}")
