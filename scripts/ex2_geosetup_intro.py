# -*- coding: utf-8 -*-
"""
geonum example script 1

@author: Jonas Gliß
@email: jg@nilu.no
@Copyright: Jonas Gliß
"""


import geonum
from matplotlib.pyplot import show

s = geonum.GeoSetup() 
s.create_test_data() #Etna test data

s.plot_2d()
s.plot_3d()

show()

