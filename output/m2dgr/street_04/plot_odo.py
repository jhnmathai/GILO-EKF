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

# Create KML file with a path
def create_kml_with_path(timestamps, lat_lon_alt):
    kml = ET.Element('kml')
    kml.set('xmlns', 'http://www.opengis.net/kml/2.2')
    document = ET.SubElement(kml, 'Document')
    
    # Define a style for the line
    style = ET.SubElement(document, 'Style', id='lineStyle')
    linestyle = ET.SubElement(style, 'LineStyle')
    color = ET.SubElement(linestyle, 'color')
    color.text = 'ff0000ff'  # Blue color (in AABBGGRR format)
    width = ET.SubElement(linestyle, 'width')
    width.text = '4'  # Width of the line
    
    # Create a Placemark for the path
    placemark = ET.SubElement(document, 'Placemark')
    name = ET.SubElement(placemark, 'name')
    name.text = 'Trajectory Path'
    styleUrl = ET.SubElement(placemark, 'styleUrl')
    styleUrl.text = '#lineStyle'
    lineString = ET.SubElement(placemark, 'LineString')
    tessellate = ET.SubElement(lineString, 'tessellate')
    tessellate.text = '1'
    coordinates = ET.SubElement(lineString, 'coordinates')
    
    # Add the coordinates for the path
    coord_string = ""
    for (lat, lon, alt) in lat_lon_alt:
        coord_string += f'{lon},{lat},{alt} '
    
    coordinates.text = coord_string.strip()
    
    tree = ET.ElementTree(kml)
    tree.write('output_path.kml', xml_declaration=True, encoding='utf-8')

# Main process
def main():
    file_name = '/home/mathjo/msr_ws/src/GILO-EKF/output/m2dgr/street_04/odo.txt'
    timestamps, ecef_data = read_data(file_name)
    
    lat_lon_alt = [ecef_to_llh(ecef) for ecef in ecef_data]
    
    create_kml_with_path(timestamps, lat_lon_alt)
    print('KML file generated as output_path.kml')

if __name__ == "__main__":
    main()
