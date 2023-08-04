#!/usr/bin/env python3

import utm, copy, shapely, re, sys
from shapely.geometry import Polygon,MultiPolygon,GeometryCollection,LineString
from shapely import ops
import airspace_config
import matplotlib.pyplot as plt

IGNORE_LIMIT = 0.005  # If an area airsport area intersects an a TMA with less than this
                      # ratio we will ignore it to not create too many tiny areas
IGNORE_SIZE = 20 * 1000 * 1000  # The intersecting area also needs to be smaller than this (meters squared)

FILENAME = 'NorwayAirspace 20230425 revA-fixed.txt'
#FILENAME = 'luftrom-2023.fl.txt'
SECTORS_FILENAME = 'acc-sectors.txt'
OUTPUT='Norway2023-modified.txt'

PLOT_SUBTRACTIONS=False

def to_dms(dd):
    mnt,sec = divmod(dd*3600, 60)
    deg,mnt = divmod(mnt, 60)
    return map(int, (deg,mnt,sec))

class Coordinate:
    def __init__(self, lat=0.0, lon=0.0):
        if lat < 0.0:
            raise Exception()
        self.lat = lat
        self.lon = lon

    def to_utm(self):
        """Assumes zone 32V which should work ok for Norway"""
        e,n,_,_ = utm.from_latlon(self.lat, self.lon, force_zone_number=32, force_zone_letter='V')
        return e,n
    
    def to_dms(self):
        ew = 'E' if self.lon >= 0 else 'W'
        return '{0:02d}:{1:02d}:{2:02d} N  {3:02d}:{4:02d}:{5:02d} {ew}'.format(*to_dms(self.lat), *to_dms(abs(self.lon)), ew=ew)

    def tuple(self):
        return self.lat, self.lon
    
    @staticmethod
    def from_string(dms_string):
        lat,lon = dms_string.split('N')
        lon = lon.replace('E', '')

        d,m,s = lat.split(':')
        d,m = map(int, (d,m))
        s = float(s)
        lat = d + m/60.0 + s / 3600.0
        
        d,m,s = lon.split(':')
        d,m = map(int, (d,m))
        s = float(s)
        lon = d + m/60.0 + s / 3600.0

        return Coordinate(lat,lon)

    @staticmethod
    def from_utm(e, n):
        lat,lon = utm.to_latlon(e, n, 32, 'V', strict=False)
        return Coordinate(lat,lon)

    @staticmethod
    def from_aip(coord):
        lat,lon = [c.strip() for c in coord.split()]
        d,m,s = [int(c) for c in (lat[:2],lat[2:4],lat[4:6])]
        lat = d + m/60.0 + s / 3600.0
        
        d,m,s = [int(c) for c in (lon[:3],lon[3:5],lon[5:7])]
        lon = d + m/60.0 + s / 3600.0       
        return Coordinate(lat, lon)

    def __repr__(self):
        return '{0:.4f}N {1:.4f}E'.format(self.lat, self.lon)

# To check that all frequencies in config file are used
# (could catch error where an airspace has changed name from config)
found_frequencies = set()

class Airspace:
    FREQUENCY_RE = re.compile(r'(.*)(\d\d\d\.\d{1,3})$')
    
    def __init__(self):
        # Fields from OpenAir format
        self.comment = None
        self.cls = None
        self.name = None
        self.key = None
        self.limit_low = None
        self.limit_high = None
        self.center = None
        self.radius = None
        self.arc_radius = None
        self.arc = None
        self.frequency = None
        self.controller = None
        self.coordinates = []

        # Other fields
        self.airsport_names = []  # Airsport areas which are part of this TMA
        self.area = None          # Shapely area representing this airspace
        self.containing_tma = None  # For airsport areas, the controlled airspaces they are part of
        self.split_areas = []      # Sometimes areas need to be split, and then this will link to the split
                                   # airspaces which replace self

    def add_frequency(self):
        # Look up frequency in config file
        frequency,airspace_name = airspace_config.lookup_frequency(self.name)

        # To check that all frequencies in config file are used
        # (could catch error where an airspace has changed name from config)
        if frequency:
            found_frequencies.add(airspace_name)

        if not self.frequency:
            # If freqency is written as last part of airspace name, them prefer that
            frequency_match = self.FREQUENCY_RE.match(self.name)
            if frequency_match:
                self.frequency = frequency_match.group(2)
            elif frequency:
                # use looked up frequency from config file
                self.frequency = frequency
                # Add frequency to end of name of airspace
                self.name += f" {self.frequency}"
                                   
    def process_area(self, *airspaces):
        """Creates shapely areas for those airspaces we need to work on.
           Converts from latitude/longitude to UTM coordinates in the process, so
           shapely can be used."""
        for space in airspaces:
            if self.key in space and not self.area:
                self.area = Polygon([c.to_utm() for c in self.coordinates])
        return self.area

    def intersect_area(self, splitting_airspace):
        """Returns a new Airspace which is the intersection with splitting_airspace"""

        if splitting_airspace.area.disjoint(self.area) or splitting_airspace.area.touches(self.area):
            # No point in subtracting if the areas don't overlap
            return None
    
        try:
            intersection = self.area.intersection(splitting_airspace.area)
        except shapely.errors.TopologicalError as e:
            # If there is an error, make a graphical plot so it's easier to see what is wrong
            plt.title(self.name + ' intersecting ' + splitting_airspace.name)
            self.plot(show_points=True)
            splitting_airspace.plot('r-')
#            c = Coordinate.from_utm(437819.00011702196, 6826375.3107311446)
#            print(c.to_dms())
#            plt.plot(c.lon, c.lat, 'ro')
            plt.show()

        
        # If the intersection area is tiny, or covers almost all of the original area,
        # then don't bother with intersecting.
        # This avoids tiny areas which are just a nuisance.  The airspace definitions are
        # unfortunately not accurate enough so that edges for airsport areas completely
        # match the containing TMA airspace
        diff = intersection.area / self.area.area
        if diff < IGNORE_LIMIT and intersection.area < 20.0 * 1000 * 1000: # or diff > (1.0 - IGNORE_LIMIT):
            return None

        # Create a new airspace object for the intersecting part of the airspace.
        # Only the intersections will be used in final file format
        new = copy.deepcopy(self)
        new.name = splitting_airspace.name
        new.comment = '* Part of {} which is intersecting {} at {}'.format(new.name, splitting_airspace.name, splitting_airspace.limit_low)
        new.area = intersection
        new.split_areas = []
        new.frequency = splitting_airspace.frequency

        return new
            
    def plot(self, color='b-', area=None, show_points=False):
        """Make a plot of an area.  Use plt.show() to show the plot later"""
        if not area:
            area = self.area

        # Get all the coordinates of all outlines we need to plot
        coords = []
        if isinstance(area, MultiPolygon) or isinstance(area, GeometryCollection):
            for a in area.geoms:
                if isinstance(a, LineString):
                    continue
                #                if a.interiors:
                #                    print('Has holes:', self.name)
                coords.append(a.exterior.coords)
        else:
            coords.append(area.exterior.coords)

        # Make a plot for each outline
        for c in coords:
            # Transform coordinates to latitude and logitude
            coord = [Coordinate.from_utm(e,n).tuple() for e,n in c]
            # Plot the outline
            plt.plot(*reversed(tuple(zip(*coord))), color)
            # Plot coordinate
            if show_points:
                for cor in coord:
                    plt.plot(*reversed(cor), color[0] + 'o')

    def sectorize(self, sector):
        """Intersects an airspace with a sectors, so that it can be divided up into sectors"""
        sector_area = self.intersect_area(sector)
        if not sector_area:
            return

        self.split_areas.append(sector_area)
                    
    def subtract(self, airspace):
        """Subtracts the @airspace from this airspace's area"""

        if self.area.area == 0.0:
            # No area left
            return

        splitoff_airspace = self.intersect_area(airspace)
        if not splitoff_airspace:
            # Areas don't intersect, nothing more to do here
            return

        # The splitting airspace (airsport area) has this as the containing TMA
        airspace.containing_tma = self

        if PLOT_SUBTRACTIONS:
        #if self.name.startswith('Farris TMA 6'):
            plt.title(self.name + ' minus ' + airspace.name)
            self.plot()
            airspace.plot('r-')

        orig_area_size = self.area.area

        # Uncomment to convert UTM co-ordinates from shaply when there are geometries with errors
        #    c = Coordinate.from_utm(655753.67995712569, 7005045.6455644593)
        #    print(airspace.name, c.to_dms())

        # The intersecting area is split off from this airspace
        self.area = self.area.difference(splitoff_airspace.area)

        if isinstance(self.area, MultiPolygon):
            # Remove tiny nuisance areas created by the subtraction
            filtered = []
            for a in self.area.geoms:
                if a.area / orig_area_size > IGNORE_LIMIT:
                    filtered.append(a)
            if len(filtered) == 1:
                self.area = filtered[0]
            else:
                self.area = MultiPolygon(filtered)

        if PLOT_SUBTRACTIONS:
#        if self.name.startswith('Farris TMA 6'):
            plt.title(self.name + ' minus ' + airspace.name)
            plt.show()

                
    def remove_holes(self):
        if not self.area or not self.area.interiors:
            return
        
        nearest = ops.nearest_points(self.area.exterior, self.area.interiors)
        e_coords = [x for x in self.area.exterior.coords]
        i_coords = [x for x in self.area.interiors[0].coords]

        print(self.area.exterior.distance(self.area.interiors))


    def __str__(self):
        """A string in OpenAir format representing this airspace"""
        if self.split_areas:
            # One entry for each split off sub area
            return ''.join([str(a) for a in self.split_areas])

        if self.area and self.area.area == 0.0:
            return ''

        lines = []

        if self.containing_tma:
            assert not self.containing_tma.split_areas, "Overlying TMA should not itself have split off areas, it is itself a split off area from the main TMA"

            # Create an entry representing the TMA over the airsport area
            overlying_tma = copy.copy(self.containing_tma)
            overlying_tma.comment = "Part of {} lying over {}".format(overlying_tma.name, self.name)
            overlying_tma.area = self.area                # Same area as air sport area
            overlying_tma.limit_low = self.limit_high
            overlying_tma.limit_high = self.containing_tma.limit_high


            # Set lower limit of airsport area to be lower limit of TMA, not ground
            self.limit_low = self.containing_tma.limit_low
        else:
            overlying_tma = None
        
        if not (self.cls and self.name and self.limit_low and self.limit_high) or\
           not (self.coordinates or self.center or self.arc_radius):
            print ('Missing information! ' + self.name)
            return ''

        if self.area is not None:
            if isinstance(self.area, MultiPolygon) or isinstance(self.area, GeometryCollection):
                # Create one entry for each area in multipolygon
                coordinates = []
                for a in self.area.geoms:
                    if not isinstance(a, Polygon):
                        continue
                    if a.interiors:
                        print('Has holes:', self.name)
        
                    coordinates.append([Coordinate.from_utm(e,n) for e,n in a.exterior.coords])
            else:
                # We don't handle holes, so warn about those
                if self.area.interiors:
                    print('Has holes:', self.name)

                coordinates = [[Coordinate.from_utm(e,n) for e,n in self.area.exterior.coords]]
        else:
            # If we didn't change the area, then use the original coordinates read in
            coordinates = [self.coordinates]

        # Create OpenAir entry for each set of coordinates
        for coords in coordinates:
            if self.comment:
                lines.append(self.comment)
            lines.append('AC ' + self.cls)
            lines.append('AN ' + self.name)
            if self.frequency:
                lines.append('AF ' + self.frequency)
            if self.controller:
                lines.append('AG ' + self.controller)
            lines.append('AL ' + self.limit_low)
            lines.append('AH ' + self.limit_high)
            if self.center:
                lines.append('V ' + self.center)
            if self.arc_radius:
                lines.append('DA ' + self.arc_radius)
            if self.arc:
                lines.append('DB ' + self.arc)
            if self.radius:
                lines.append('DC ' + self.radius)

            for c in coords:
                lines.append('DP ' + c.to_dms())
            lines.append('\n*')

        res = '\n'.join(lines) + '\n'

        # Create another entry for the overlying TMA, if there is one
        if overlying_tma:
            res += str(overlying_tma)
        
        return res
        
def parse(filename):
    airspace = None
    airspaces = []
    content = []     # Output content, contains either strings, or Airspace objects which can be printed
    
    polaris_airspaces = {}     # Polaris airspace
    tma_airspaces = {}         # Each TMA which has aerial sport areas
    airsport_airspaces = {}    # Link to TMA airspace from each aerial sport areas

    for tma,areas in airspace_config.airsport_areas.items():
        tma_airspaces[tma.upper()] = []
        for a in areas:
            airsport_airspaces[a] = []
    

    for line in open(filename, encoding='iso8859-1'):
#    for line in open(filename, encoding='utf-8'):        
        line = line.strip()
        if len(line) < 2 or line.startswith('*'):
            content.append(line)
            continue

        field,rest = line[:2],line[2:].strip()

        if field == 'AC':
            airspace = Airspace()
            airspace.cls = rest
            airspaces.append(airspace)
            content.append(airspace)
        elif not airspace:
            continue
        elif field == 'AN':
            airspace.name = rest

            if airspace.name.upper().startswith('POLARIS'):
                airspace.key = airspace.name
                polaris_airspaces[airspace.name] = airspace

            # Check if this TMA has any airsport areas
            for tma,airsport in airspace_config.airsport_areas.items():
                tma = tma.upper()
                if airspace.name.upper().startswith(tma):
                    airspace.airsport_names = airsport
                    airspace.key = tma
                    tma_airspaces[tma].append(airspace)
                    
                for name in airsport:
                    if airspace.name.upper().startswith(name.upper()):
                        airspace.key = name
                        airsport_airspaces[name].append(airspace)
        elif field == 'AL':
            airspace.limit_low = rest
        elif field == 'AH':
            airspace.limit_high = rest
        elif field == 'AF':
            airspace.frequency = rest
        elif field == 'AG':
            airspace.controller = rest
        elif field == 'DP':
            coord = Coordinate.from_string(rest)
            airspace.coordinates.append(coord)
        elif field == 'V ':
            airspace.center = rest
        elif field == 'DA':
            airspace.arc_radius = rest
        elif field == 'DB':
            airspace.arc = rest
        elif field == 'DC':
            airspace.radius = rest
        else:
            print('Unhandled: ', line)

    r = re.compile(r'((.*)(TMA|CTA|TIA|TIZ|CTR)\s+(\d*))')
    freqs = []
    for a in airspaces:
        a.process_area(airsport_airspaces, tma_airspaces, polaris_airspaces)

        if False and a.frequency:
            name = a.name
            if m := a.FREQUENCY_RE.match(name):
                name = m.group(1).strip()
            
            if m := r.match(name):
                name = m.group(0)
            freqs.append((f"    '{name}' : '{a.frequency}',"))
        a.add_frequency()
        
    # Warn about airspaces in config file which we didn't not find in airspace file
    for freq_name in airspace_config.airspace_frequencies.keys():
        if freq_name not in found_frequencies:
            continue
            print(f'Warning: Airspace from frequencies config file not found in airspace file: {freq_name}')

    return content, tma_airspaces, polaris_airspaces, airsport_airspaces

def parse_acc_sectors(filename):
    """Parses a simple file with ACC sectors.
    This is unexplainably not part of AIP, so we have a special file just with sector coordinates"""
    
    name_re = re.compile(r'ACC SECTOR (\d+)')
    a = Airspace()
    spaces = []
    
    for line in open(filename):
        if match := name_re.match(line):
            a = Airspace()
            a.name = f'Polaris S{match.group(1)}'
            a.key = a.name
            spaces.append(a)
        elif line.strip():  # ignore blank lines
            coord = Coordinate.from_string(line)
            a.coordinates.append(coord)

    for a in spaces:
        a.process_area([a.key])
        a.add_frequency()
            
    return spaces

def subtract_airsport_airspaces():
    # Go through all TMAs with airsport areas
    for tma in tma_airspaces.values():
        # Go through each sub area in TMA
        for airspace in tma:
            # Go through and subtract all airsport airspaces
            for airsport in airspace.airsport_names:
                airsports = airsport_airspaces[airsport]
                # Subtract each area in airsport airspace
                for area in airsports:
                    area.cls = airspace.cls
                    airspace.subtract(area)


def sectorize_polaris(areas, sectors):
    # Go through all Polaris areas and divide into sectors
    for area in areas:
        for sector in sectors:
            area.sectorize(sector)

test_data = \
"""592953N 0093542E - 592426N 0093524E - 592112N 0092731E - 592255N 0091923E - 592953N 0093542E"""

def parse_aip_coordinates(data):
    fields = [c.strip() for c in data.split('-')]
    return [Coordinate.from_aip(f) for f in fields]

def convert_aip_to_openair(data):
    """Used to convert coordinates copied directly from AIP into OpenAir format"""
    coords = parse_aip_coordinates(data)
    for c in coords:
        print('DP ' + c.to_dms())

#print(convert_aip_to_openair(test_data))

if __name__ == '__main__':
    #convert_aip_to_openair(test_data)
    polaris_sectors = parse_acc_sectors(SECTORS_FILENAME)

    content, tma_airspaces, polaris_airspaces, airsport_airspaces = parse(FILENAME)
    sectorize_polaris(polaris_airspaces.values(), polaris_sectors)
    
    subtract_airsport_airspaces()

    output = open(OUTPUT, 'w')
    for line in content:
        output.write(str(line) + '\n')
    output.close()
    print('Result written to', OUTPUT)
