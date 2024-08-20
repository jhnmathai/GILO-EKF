from simplekml import Kml
import numpy as np

# File path (replace this with your actual path)
filename = '/home/mathjo/catkin_ws/src/GILO-EKF/output/m2dgr/street_05/odo.txt'

# Load the data from the file, ignoring comment lines
data = np.loadtxt(filename, comments='#')

# Extract BLH columns (latitude, longitude in radians and height)
latitudes_rad = data[:, 1]
longitudes_rad = data[:, 2]
heights = data[:, 3]

# Convert radians to degrees for latitude and longitude
latitudes_deg = np.degrees(latitudes_rad)
longitudes_deg = np.degrees(longitudes_rad)

# Create a KML object
kml = Kml()

# Add the trajectory to the KML file as a line
linestring = kml.newlinestring(name="Trajectory")
linestring.coords = [(lon, lat, height) for lon, lat, height in zip(longitudes_deg, latitudes_deg, heights)]
linestring.altitudemode = "clampToGround"  # Clamp the line to the terrain
linestring.extrude = 0  # Disable extrusion (no vertical line to the ground)
linestring.style.linestyle.color = "ff0000ff"  # Red color (in KML format)
linestring.style.linestyle.width = 3  # Line width


# Save the KML file
kml_file_path = '/home/mathjo/catkin_ws/src/GILO-EKF/output/m2dgr/street_05/trajectory.kml'
kml.save(kml_file_path)

print(f"KML file saved as 1'{kml_file_path}'.")
