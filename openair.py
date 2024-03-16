#!/usr/bin/env python3

import utm, copy, shapely, re, sys, requests, time
from shapely.geometry import Polygon,MultiPolygon,GeometryCollection,LineString
from shapely import ops, set_precision
import airspace_config
import matplotlib.pyplot as plt

IGNORE_LIMIT = 0.0004  # If an area airsport area intersects an a TMA with less than this
                      # ratio we will ignore it to not create too many tiny areas
IGNORE_SIZE = 0.5 * 1000 * 1000  # The intersecting area also needs to be smaller than this (meters squared)
PRECISION = 5.0   # Precision used in calculations in m.  If this is too small it will create
                   # almost empty "lines" when doing the area subtractions due to rounding errors

URL = 'https://raw.githubusercontent.com/relet/pg-xc/master/openair/luftrom.fl.txt'
LOCAL_ADDITIONS = 'static/local-additions.txt'
#FILENAME="polaris.txt"
SECTORS_FILENAME = 'static/acc-sectors.txt'
CHANGELOG_FILENAME = 'static/changelog.txt'
OUTPUT='Norway2023-modified.txt'

PLOT_SUBTRACTIONS=False
USE_EXTENDED_OPENAIR=False

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

    def to_openair(self):
        return 'DP ' + self.to_dms()

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

class ArcSegment:
    def __init__(self, coord1, coord2):
        self.coord1 = coord1
        self.coord2 = coord2

    def to_openair(self):
        return 'DB ' + self.coord1.to_dms() + ', ' + self.coord2.to_dms()

class Circle:
    def __init__(self, radius):
        self.radius = radius

    def to_openair(self):
        return 'DC ' + self.radius

class Parameter:
    def __init__(self, direction=None, center=None, width=None, zoom=None):
        self.direction = direction
        self.center = center
        self.width = width
        self.zoom = zoom
    
    def to_openair(self):
        if self.direction:
            return f'V D={self.direction}'
        if self.center:
            return f'V X={self.center.to_dms()}'
        if self.width:
            return f'V W={self.width}'
        if self.zoom:
            return f'V Z={self.zoom}'

# To check that all frequencies in config file are used
# (could catch error where an airspace has changed name from config)
found_frequencies = set()

class Airspace:
    FREQUENCY_RE = re.compile(r'(.*)(\d\d\d\.\d{1,3})$')
    
    def __init__(self):
        # Fields from OpenAir format
        self.comment = None
        self.cls = None
        self.type = None
        self.name = None
        self.identifier = None
        self.key = None
        self.limit_low = None
        self.limit_high = None
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
        frequency,controller,airspace_name = airspace_config.lookup_frequency(self.name)

        # To check that all frequencies in config file are used
        # (could catch error where an airspace has changed name from config)
        if frequency:
            found_frequencies.add(airspace_name)

            if not self.frequency:
                # use looked up frequency from config file
                self.frequency = frequency
                # Add frequency to end of name of airspace
                self.name += f" {self.frequency}"
                self.controller = controller
                                   
    def process_area(self, *airspaces):
        """Creates shapely areas for those airspaces we need to work on.
           Converts from latitude/longitude to UTM coordinates in the process, so
           shapely can be used."""
        for space in airspaces:
            if self.key in space and not self.area:
                self.area = set_precision(Polygon([c.to_utm() for c in self.coordinates]), PRECISION)
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
        if diff < IGNORE_LIMIT:
            #print(f"Not splitting {self.name} and {splitting_airspace.name} because the diff is too low: {diff}")
            return None
        if intersection.area < IGNORE_SIZE: # or diff > (1.0 - IGNORE_LIMIT):
            #print(f"Not splitting {self.name} and {splitting_airspace.name} because size of new area too low {intersection.area}")
            return None

        # Create a new airspace object for the intersecting part of the airspace.
        # Only the intersections will be used in final file format
        new_containing = copy.deepcopy(self)
        new_containing.name = splitting_airspace.name
        new_containing.comment = f'Part of {self.name} which is intersecting {splitting_airspace.name}'
        if splitting_airspace.limit_low:
            new_containing.comment += f'at {splitting_airspace.limit_low}'
        new_containing.area = intersection
        new_containing.split_areas = []
        new_containing.frequency = splitting_airspace.frequency
        new_containing.controller = splitting_airspace.controller

        # If the splitting airsport area does is not entirely lie within self's airspace,
        # then we need to split the airsport area too
        new_airsport = copy.deepcopy(splitting_airspace)
        new_airsport.name = splitting_airspace.name
        new_airsport.comment = f'Part of {splitting_airspace.name} which lies inside {self.name} {self.limit_low}'
        new_airsport.area = intersection
        new_airsport.limit_low = self.limit_low
        new_airsport.containing_tma = None  # Assumes that the same containing sector lies over the entire airsport area
        new_airsport.split_areas = []
        splitting_airspace.split_areas.append(new_airsport)

        return new_containing

    def sectorize(self, sector):
        """Intersects an airspace with sectors, so that it can be divided up into sectors"""
        sector_area = self.intersect_area(sector)
        if not sector_area:
            return

        self.split_areas.append(sector_area)
                    
    def subtract(self, airspace):
        """Subtracts the @airspace from this airspace's area"""

        if self.area.area == 0.0:
            # No area left
            return

        # If the area we are subtracting from already consist of several sub areas, then
        # we need to subtract from each of them.
        # (typically Polaris has been split into sectors, which we then have to subtract wave
        #  sectors from)
        if self.split_areas:
            subareas = self.split_areas
        else:
            subareas = [self]

        for subarea in subareas:
            splitoff_airspace = subarea.intersect_area(airspace)

            if not splitoff_airspace:
                # Areas don't intersect, nothing more to do here
                continue

            if PLOT_SUBTRACTIONS:
                plt.title(subarea.name + ' minus ' + airspace.name)
                subarea.plot()
                airspace.plot('r-')
                plt.title(subarea.name + ' minus ' + airspace.name)
                plt.show()

            # The splitting airspace (airsport area) has this as the containing TMA
            airspace.containing_tma = subarea

            orig_area_size = subarea.area.area

            # Uncomment to convert UTM co-ordinates from shaply when there are geometries with errors
            #    c = Coordinate.from_utm(655753.67995712569, 7005045.6455644593)
            #    print(airspace.name, c.to_dms())

            # The intersecting area is split off from this airspace
            subarea.area = subarea.area.difference(splitoff_airspace.area)

            if isinstance(subarea.area, MultiPolygon):
                # Remove tiny nuisance areas created by the subtraction
                filtered = []
                for a in subarea.area.geoms:
                    if a.area / orig_area_size > IGNORE_LIMIT:
                        filtered.append(a)
                if len(filtered) == 1:
                    subarea.area = filtered[0]
                else:
                    subarea.area = MultiPolygon(filtered)

    def remove_holes(self):
        if not self.area or not self.area.interiors:
            return
        
        nearest = ops.nearest_points(self.area.exterior, self.area.interiors)
        e_coords = [x for x in self.area.exterior.coords]
        i_coords = [x for x in self.area.interiors[0].coords]

        print(self.area.exterior.distance(self.area.interiors))


    def __str__(self):
        """A string in OpenAir format representing this airspace"""

        if self.area is not None and self.area.is_empty:
            assert not self.split_areas
            return ''

        if self.containing_tma:
            #assert not self.containing_tma.split_areas, ("Overlying TMA (%s) should not itself have split off areas, it is itself a split off area from the main TMA" % {self.containing_tma.name})

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

        if self.split_areas:
            # One entry for each split off sub area
            res = ''.join([str(a) for a in self.split_areas])

            # Add one entry for overlying airspace
            # (assumes the same airspace lies over all of the the airsport area)
            if overlying_tma:
                return str(overlying_tma) + res
            else:
                return res

        lines = []

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
                lines.append('* ' + self.comment)
            lines.append('AC ' + self.cls)
            if USE_EXTENDED_OPENAIR:
                lines.append('AY ' + self.type)
            lines.append('AN ' + self.name)
            if USE_EXTENDED_OPENAIR:
                lines.append('AI ' + self.identifier)
            lines.append('AL ' + self.limit_low)
            lines.append('AH ' + self.limit_high)
            if self.frequency:
                lines.append('AF ' + self.frequency)
                if self.controller:
                    lines.append('AG ' + self.controller)

            for c in coords:
                lines.append(c.to_openair())

            lines.append('\n*')

        res = '\n'.join(lines) + '\n'

        # Create another entry for the overlying TMA, if there is one
        if overlying_tma:
            res += str(overlying_tma)
        
        return res

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

def parse(*filenames):
    airspace = None
    airspaces = []
    content = []     # Output content, contains either strings, or Airspace objects which can be printed
    
    polaris_airspaces = {}     # Polaris airspace
    tma_airspaces = {}         # Each TMA which has aerial sport areas
    airsport_airspaces = {}    # Link to TMA airspace from each aerial sport areas

    numbered_airspace_re = re.compile(r'(.*(TMA|TIA|CTA)) \d+')

    # Prepare airspace config first
    for tma,areas in airspace_config.airsport_areas.items():
        tma_airspaces[tma.upper()] = []
        for a in areas:
            airsport_airspaces[a] = []
    
    for filename in filenames:
        for line in open(filename, encoding='utf-8'):        
            line = line.strip()
            if len(line) < 2 or line.startswith('*'):
                content.append(line)
                continue

            field,rest = line[:2],line[2:].strip()

            if field == 'AC':
                airspace = Airspace()
                cls = rest
                if USE_EXTENDED_OPENAIR:
                    # In extended OpenAir class should purely be ICAO class, so
                    # change to class G what is not A,B,C,D,E (concerns warning areas, restrictions areas)
                    if cls not in 'ABCDEG':
                        cls = 'G'
                airspace.cls = cls
            elif not airspace:
                continue
            elif field == 'AY':
                airspace.type = rest
            elif field == 'AN':
                airspace.name = rest

                if airspace.name in airspace_config.REMOVE_AIRSPACES:
                    # ignore this airspace
                    airspace = None
                    continue

                # Add to our airspaces
                airspaces.append(airspace)
                content.append(airspace)

                if not airspace.type and USE_EXTENDED_OPENAIR:
                    # Figure out type from name
                    if airspace.name.startswith('EN R'):
                        airspace.type = 'R'
                    elif airspace.name.startswith('EN D') or airspace.name.startswith('END'):
                        airspace.type = 'Q'
                    elif 'CTA' in airspace.name:
                        airspace.type = 'CTA'
                    elif 'CTR' in airspace.name:
                        airspace.type = 'CTR'
                    elif 'TMA' in airspace.name:
                        airspace.type = 'TMA'
                    elif 'TIA' in airspace.name or 'TIZ' in airspace.name:
                        airspace.type = 'RMZ'
                    elif airspace.cls == 'G':
                        airspace.type = 'Q'
                    else:
                        print(f'Unable to find airspace type for {airspace.name}')
                
                if airspace.name.upper().startswith('POLARIS'):
                    airspace.key = airspace.name
                    polaris_airspaces[airspace.name] = airspace

                # Create an ID from the name
                m = numbered_airspace_re.match(airspace.name) 
                if m:
                    # Drop the number
                    airspace.identifier = m.group(1)
                else:
                    airspace.identifier = airspace.name

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
                parameter,value = rest.split('=', 1)
                if parameter == 'X':
                    coord = Coordinate.from_string(value)
                    airspace.coordinates.append(Parameter(center=coord))
                elif parameter == 'D':
                    airspace.coordinates.append(Parameter(direction=value))
                elif parameter == 'W':
                    airspace.coordinates.append(Parameter(width=value))
                elif parameter == 'Z':
                    airspace.coordinates.append(Parameter(zoom=value))
            elif field == 'DA':
                raise Exception()
                airspace.arc_radius = rest
            elif field == 'DB':
                coord1,coord2 = rest.split(',', 1)
                airspace.coordinates.append(ArcSegment(
                                            Coordinate.from_string(coord1),
                                            Coordinate.from_string(coord2)))
            elif field == 'DC':
                airspace.coordinates.append(Circle(rest))
            else:
                print('Unhandled: ', line)

    # Add radio frequencies to airspaces
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
            a.identifier = a.name
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
        print(c.to_openair())

#print(convert_aip_to_openair(test_data))

def download_luftrom_info():
    r = requests.get(URL, allow_redirects=True)
    filename = 'Downloads/luftrom.fl.txt'
    open(filename, 'wb').write(r.content)
    return filename

if __name__ == '__main__':
    #convert_aip_to_openair(test_data)
    polaris_sectors = parse_acc_sectors(SECTORS_FILENAME)

    luftrom = download_luftrom_info()
    content, tma_airspaces, polaris_airspaces, airsport_airspaces = parse(
        luftrom, LOCAL_ADDITIONS)
    sectorize_polaris(polaris_airspaces.values(), polaris_sectors)
    
    subtract_airsport_airspaces()

    output = open(OUTPUT, 'w')
    output.write(open(CHANGELOG_FILENAME).read())

    empty_line = False
    for line in content:
        content = str(line).strip()
        if not content:
            # Skip multiple empty lines
            if empty_line:
                continue
            empty_line = True
        else:
            empty_line = False

        output.write(content + '\n')
    output.close()
    print('Result written to', OUTPUT)
