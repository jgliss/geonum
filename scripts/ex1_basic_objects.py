"""
geonum example script 1
"""

from geonum.base import GeoPoint

### Create 2 GeoPoint objects 

# Summit region of Mt. Etna
p1 = GeoPoint(37.751005,  14.993435, altitude = 3264.0, name = "Etna")
# Position of volcanological observatory of Mt. Etna
#p2 = GeoPoint(37.765755,  15.016696, altitude = 2820.0, name = "Observatory")
p2 = GeoPoint(37.765755,  15.016696, name = "Observatory")
# Print info (string represenation of both points")
print "Point1: %s, point 2: %s" %(p1, p2)

# Get the connection vector of both points
connection_vector = p2 - p1

# The connection vector of both objects
print connection_vector

# Output:
# GeoVector3D Etna->Observatory
# Azimuth: 51.3786772497, Elevation: -9.60641740857, Magnitude: 2.66060747963 m


