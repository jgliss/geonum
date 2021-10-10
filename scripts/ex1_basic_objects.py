"""geonum example script 1

Introduction into the geonum base classes GeoPoint and GeoVector3D
"""

from geonum import GeoPoint
from SETTINGS import OPTPARSE
from numpy import testing as npt
from matplotlib.pyplot import show
### Create 2 GeoPoint objects 

# Summit region of Mt. Etna
p1 = GeoPoint(latitude=37.751005, longitude=14.993435,
              altitude=3264.0, name="Etna")

# Position of volcanological observatory of Mt. Etna
p2 = GeoPoint(latitude=37.765755, longitude=15.016696,
              auto_topo_access=True,
              name="Observatory")

# Print info (string represenation of both points")
print(("Point1: %s, point 2: %s" %(p1, p2)))

# Get and print the connection vector of both points
connection_vector = p2 - p1

print(connection_vector)

# Import script options
(options, args) = OPTPARSE.parse_args()

# If applicable, do some tests. This is done only if TESTMODE is active: 
# testmode can be activated globally (see SETTINGS.py) or can also be 
# activated from the command line when executing the script using the 
# option --test 1
if int(options.test):
    from os.path import basename
    
    actual = [connection_vector.azimuth,
              connection_vector.elevation, 
              connection_vector.magnitude]
    assert connection_vector.anchor is p1
    npt.assert_allclose(actual=actual,
                        desired=[51.378677249653983,
                                 -9.6064174085658465,
                                 2.6606074796318557],
                        rtol=1e-7)
    
    

    print(("All tests passed in script: %s" %basename(__file__))) 
try:
    if int(options.show) == 1:
        show()
except:
    print("Use option --show 1 if you want the plots to be displayed")



# Output:
# GeoVector3D Etna->Observatory
# Azimuth: 51.38°, Elevation: -9.6064°, Magnitude: 2.66 km (hor: 2.62 km)


