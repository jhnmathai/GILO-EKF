import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyproj import Transformer

# Function to convert ECEF to UTM
def ecef_to_utm(x, y, z):
    transformer_to_latlong = Transformer.from_crs("epsg:4978", "epsg:4326", always_xy=True)  # ECEF to lat/lon
    transformer_to_utm = Transformer.from_crs("epsg:4326", "epsg:32633", always_xy=True)  # lat/lon to UTM (zone 33)

    lon, lat, alt = transformer_to_latlong.transform(x, y, z)
    easting, northing = transformer_to_utm.transform(lon, lat)

    return easting, northing, alt

# Function to calculate heading
def calculate_heading(easting, northing):
    delta_e = np.diff(easting)
    delta_n = np.diff(northing)
    heading = np.arctan2(delta_n, delta_e)
    heading = np.degrees(heading) % 360
    return np.concatenate(([heading[0]], heading))  # Repeat the first heading value

# Read data
ground_truth = pd.read_csv('gt_street_04.txt', header=None, delim_whitespace=True)
ekf_1 = pd.read_csv('gilo_street_04.txt', header=None, delim_whitespace=True)
ekf_2 = pd.read_csv('lio_street_04.txt', header=None, delim_whitespace=True)

# Select relevant columns
ground_truth = ground_truth.iloc[:, [0, 1, 2, 3]]
ekf_1 = ekf_1.iloc[:, [0, 1, 2, 3]]
ekf_2 = ekf_2.iloc[:, [0, 1, 2, 3]]

# Now assign the column names
ground_truth.columns = ['timestamp', 'x', 'y', 'z']
ekf_1.columns = ['timestamp', 'x', 'y', 'z']
ekf_2.columns = ['timestamp', 'x', 'y', 'z']

# Convert ECEF to UTM
gt_easting, gt_northing, _ = ecef_to_utm(ground_truth['x'], ground_truth['y'], ground_truth['z'])
ekf1_easting, ekf1_northing, _ = ecef_to_utm(ekf_1['x'], ekf_1['y'], ekf_1['z'])
ekf2_easting, ekf2_northing, _ = ecef_to_utm(ekf_2['x'], ekf_2['y'], ekf_2['z'])

# Interpolate EKF results to match ground truth timestamps
ekf1_easting_interp = np.interp(ground_truth['timestamp'], ekf_1['timestamp'], ekf1_easting)
ekf1_northing_interp = np.interp(ground_truth['timestamp'], ekf_1['timestamp'], ekf1_northing)
ekf2_easting_interp = np.interp(ground_truth['timestamp'], ekf_2['timestamp'], ekf2_easting)
ekf2_northing_interp = np.interp(ground_truth['timestamp'], ekf_2['timestamp'], ekf2_northing)

# Calculate heading
gt_heading = calculate_heading(gt_easting, gt_northing)
ekf1_heading = calculate_heading(ekf1_easting_interp, ekf1_northing_interp)
ekf2_heading = calculate_heading(ekf2_easting_interp, ekf2_northing_interp)

# Plotting the results
plt.figure(figsize=(15, 10))

# Plot Easting vs Northing with equal scaling
plt.subplot(2, 2, 1)
plt.plot(gt_easting, gt_northing, label='Ground Truth', color='black', linestyle='--')
plt.plot(ekf1_easting_interp, ekf1_northing_interp, label='GILO-EKF', color='blue')
plt.plot(ekf2_easting_interp, ekf2_northing_interp, label='LIO-EKF', color='red')
plt.xlabel('Easting (m)')
plt.ylabel('Northing (m)')
plt.title('Trajectory (Easting vs Northing)')
plt.legend()
plt.axis('equal')  # Set equal scaling

# Convert Pandas Series to Numpy array before plotting
gt_timestamp = ground_truth['timestamp'].to_numpy()

# Plot Easting Error
plt.subplot(2, 2, 2)
plt.plot(gt_timestamp, gt_easting - ekf1_easting_interp, label='GILO-EKF Error', color='blue')
plt.plot(gt_timestamp, gt_easting - ekf2_easting_interp, label='LIO-EKF Error', color='red')
plt.xlabel('Timestamp')
plt.ylabel('Easting Error (m)')
plt.title('Easting Error')
plt.legend()

# Plot Northing Error
plt.subplot(2, 2, 3)
plt.plot(gt_timestamp, gt_northing - ekf1_northing_interp, label='GILO-EKF Error', color='blue')
plt.plot(gt_timestamp, gt_northing - ekf2_northing_interp, label='LIO-EKF Error', color='red')
plt.xlabel('Timestamp')
plt.ylabel('Northing Error (m)')
plt.title('Northing Error')
plt.legend()

# Plot Heading Error
plt.subplot(2, 2, 4)
plt.plot(gt_timestamp, gt_heading - ekf1_heading, label='GILO-EKF Error', color='blue')
plt.plot(gt_timestamp, gt_heading - ekf2_heading, label='LIO-EKF Error', color='red')
plt.xlabel('Timestamp')
plt.ylabel('Heading Error (degrees)')
plt.title('Heading Error')
plt.legend()

plt.tight_layout()
plt.show()
