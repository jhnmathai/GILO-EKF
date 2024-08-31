import pyproj
import xml.etree.ElementTree as ET

# Function to convert ECEF to latitude, longitude, and height
def ecef_to_llh(ecef_coords):
    ecef = pyproj.Proj(proj='geocent', ellps='WGS84', datum='WGS84')
    lla = pyproj.Proj(proj='latlong', ellps='WGS84', datum='WGS84')
    transformer = pyproj.Transformer.from_proj(ecef, lla)
    lon, lat, alt = transformer.transform(ecef_coords[0], ecef_coords[1], ecef_coords[2], radians=False)
    return lat, lon, alt

# Reading data from the file
def read_data(file_name):
    with open(file_name, 'r') as file:
        lines = file.readlines()
    
    ecef_data = []
    timestamps = []
    
    for line in lines:
        parts = line.split()
        timestamp = float(parts[0])
        ecef_x = float(parts[1])
        ecef_y = float(parts[2])
        ecef_z = float(parts[3])
        
        timestamps.append(timestamp)
        ecef_data.append((ecef_x, ecef_y, ecef_z))
    
    return timestamps, ecef_data

# Create KML file
def create_kml(timestamps, lat_lon_alt):
    kml = ET.Element('kml')
    kml.set('xmlns', 'http://www.opengis.net/kml/2.2')
    document = ET.SubElement(kml, 'Document')
    
    for i, (timestamp, (lat, lon, alt)) in enumerate(zip(timestamps, lat_lon_alt)):
        placemark = ET.SubElement(document, 'Placemark')
        name = ET.SubElement(placemark, 'name')
        name.text = str(i)
        description = ET.SubElement(placemark, 'description')
        description.text = f'Timestamp: {timestamp}'
        point = ET.SubElement(placemark, 'Point')
        coordinates = ET.SubElement(point, 'coordinates')
        coordinates.text = f'{lon},{lat},{alt}'
    
    tree = ET.ElementTree(kml)
    tree.write('output.kml', xml_declaration=True, encoding='utf-8')

# Main process
def main():
    file_name = '/home/mathjo/msr_ws/src/GILO-EKF/output/m2dgr/street_03/odo.txt'
    timestamps, ecef_data = read_data(file_name)
    
    lat_lon_alt = [ecef_to_llh(ecef) for ecef in ecef_data]
    
    create_kml(timestamps, lat_lon_alt)
    print('KML file generated as output.kml')

if __name__ == "__main__":
    main()
