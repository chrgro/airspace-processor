#!/usr/bin/env python3

import sys
import pprint
import shapely
import utm
from shapely.geometry import LineString

# Read in NorwayAirspace....txt or similar OpenAir format
airspace_file = sys.argv[1]

# Converting to DMS format from just a float degree value
def to_dms(dd):
    mnt,sec = divmod(dd*3600, 60)
    deg,mnt = divmod(mnt, 60)
    return map(int, (deg,mnt,sec))

# Coordinate class, stolen from Arne Martins code
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


##
## Parsing OpenAir format document from here
##

# Generic common element in the file, didn't end up using this after all
class FileElement:
    def __init__(self):
        pass
    
# An airspace definition in the file
class AirspaceDef(FileElement):
    def __init__(self):
        self.definitionlines = []

        self.ac = None
        self.an = None
        self.al = None
        self.ah = None
        self.dp = []

        self.v_or_db = []

    # Parse the lines of an airspace definition
    def parse(self, lines, i):
        while i < len(lines):
            if lines[i].startswith("*"):
                break
            else:
                self.definitionlines.append(lines[i])
                if lines[i].strip() == "":
                    pass
                if lines[i].startswith("AC"):
                    assert(self.ac == None)
                    self.ac = lines[i].replace("AC ", "").strip()
                elif lines[i].startswith("AN"):
                    assert(self.an == None)
                    self.an = lines[i].replace("AN ", "").strip()
                elif lines[i].startswith("AL"):
                    assert(self.al == None)
                    self.al = lines[i].replace("AL ", "").strip()
                elif lines[i].startswith("AH"):
                    assert(self.ah == None)
                    self.ah = lines[i].replace("AH ", "").strip()
                elif lines[i].startswith("DP"):
                    self.dp.append(lines[i].replace("DP ", "").strip())
                elif lines[i].startswith("V") or lines[i].startswith("DB"):
                    self.v_or_db.append(lines[i].strip())
                else:
                    print("Error: unknown line type: ", i, lines[i])
                    exit(1)
            i+=1
        return i

    # Reconstruct the openair file based on the definition
    def reconstruct(self):
        ret = f"AC {self.ac}\nAN {self.an}\nAL {self.al}\nAH {self.ah}\n"
        if "Wave" in self.an:
            dp_ordered = reversed(self.dp)
        else:
            dp_ordered = self.dp
        for dp in dp_ordered:
            ret += "DP "+dp+"\n"
        for v_or_db in self.v_or_db:
            ret += v_or_db+"\n"
        return ret
        
        self.v_or_db = []
        return "".join(definitionline for definitionline in self.definitionlines)

# A comment line, to be able to 1:1 reconstruct an existing file
class Comment(FileElement):
    def __init__(self):
        self.commentlines = []

    def parse(self, lines, i):
        while i < len(lines) and lines[i].startswith("*"):
            self.commentlines.append(lines[i])
            i+=1
        return i

    def reconstruct(self):
        return "".join(commentline for commentline in self.commentlines)

# Empty lines, to be able to reconstruct 1:1 original file
class EmptyLine(FileElement):
    def __init__(self):
        self.empty_lines = []

    def parse(self, lines, i):
        self.empty_lines.append(lines[i])
        return i+1

    def reconstruct(self):
        return "".join(empty_line for empty_line in self.empty_lines)

# Top element when parsing an openair .txt file
class Document(FileElement):
    def __init__(self):
        self.elements = []

    def parse(self, lines, i):
        while i < len(lines):
            if lines[i].startswith("*"):
                print("Parsing comment from line", i)
                comment = Comment()
                i = comment.parse(lines, i)
                print("Comment done at line", i)
                self.elements.append(comment)
                continue
            elif lines[i].strip() == "":
                print("empty line", i, len(lines))
                emptyline = EmptyLine()
                i = emptyline.parse(lines, i)
                self.elements.append(emptyline)
                continue
            else:
                print("Parsing airspace from line", i)
                airspace = AirspaceDef()
                i = airspace.parse(lines, i)
                print("Airspace done at line", i)
                self.elements.append(airspace)
                continue

    def __str__(self):
        return str(self.elements)

    def reconstruct(self):
        return "".join(element.reconstruct() for element in self.elements)

# Read the airspace file
# Encoding seems to be latin1 or something, its at least not UTF8
with open(airspace_file, encoding='ISO-8859-1') as fp:
    lines = fp.readlines()

# Parse the document
document = Document()
document.parse(lines, 0)

# File encoding is CRLF
with open("original_reconstructed.txt", "w", encoding='ISO-8859-1', newline='\r\n') as fp:
    fp.write(document.reconstruct())

# Find the airspace we want work on for Bjorli wavecamp usage
# Looking for 3 CTA airspaces and 3 waveboxes
cta_spaces = []
wave_spaces = []
for item in document.elements:
    if type(item) == AirspaceDef:
        if item.an in ["Polaris CTA 12", "Polaris CTA 13", "Polaris CTA 14"]:
            print("Found this airspace:", item.an)
            #pprint.pprint(vars(item))
            cta_spaces.append(item)

        elif item.an in ["Bjorli Wave 125.70", "Lesja Wave 125.70", "Dovre Wave 125.70"]:
            print("Found this wavebox airspace:", item.an)
            #pprint.pprint(vars(item))
            wave_spaces.append(item)

assert(len(cta_spaces) == 3)
assert(len(wave_spaces) == 3)

# Create a dummy airspace def so we can later visualize all crossings
all_intersects = AirspaceDef()
all_intersects.ac = "D"
all_intersects.an = "All Intersect Points"
all_intersects.al = "FL 50"
all_intersects.ah = "FL 100"

# Compute all intersecting points between a CTA and a wavebox
# We will need these points later when constructing the new airspace locations
for wave in wave_spaces:
    for cta in cta_spaces:
        #print("Processing", cta.an)
        #print("Processing", wave.an)
        print("***")
        for i in range(len(cta.dp)-1):
            cta_a = Coordinate.from_string(cta.dp[i])
            cta_b = Coordinate.from_string(cta.dp[i+1])
            cta_line = LineString([cta_a.to_utm(), cta_b.to_utm()])
            for j in range(len(wave.dp)-1):
                wave_a = Coordinate.from_string(wave.dp[j])
                wave_b = Coordinate.from_string(wave.dp[j+1])
                wave_line = LineString([wave_a.to_utm(), wave_b.to_utm()])

                cross = cta_line.intersection(wave_line)
                #print(type(cross))
                if not cross.is_empty:
                    print("Intersection between",cta.an,"and",wave.an,":")
                    cross_coord = Coordinate.from_utm(cross.x, cross.y)
                    print(cross_coord.to_dms())
                    all_intersects.dp.append(cross_coord.to_dms())
                    print("Between Wave", wave.dp[j], "--", wave.dp[j+1], "and", cta.dp[i],"--", cta.dp[i+1])
                    print()
                    
# Printout the intersections
print(all_intersects.reconstruct())
# CTA: Counter clockwise
# Wave: Clockwise

# The rest is manual work, to construct new airspace boxes!
# Use a GUI tool for this to click-and-create airspace definitions

