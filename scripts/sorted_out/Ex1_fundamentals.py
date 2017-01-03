# -*- coding: utf-8 -*-
"""
geonum example script 1

@author: Jonas Gliß
@email: jg@nilu.no
@Copyright: Jonas Gliß
"""
#==============================================================================
# import sys
# import PIL.Image
# sys.modules['Image'] = PIL.Image
#==============================================================================

import geonum 



p1=gnp.GeoPoint(-18.093200, -69.215787, name = "Observer")      
p2=gnp.GeoPoint(-18.163532, -69.14264, name="Parinacota", topoPath=topoPath)                  

vc=p2-p1

print vc            

#:DEFINE CORNERS FOR RETRIEVAL OF TOPOGRAPHY DATA
#Top left point of map
pTL=p1.offset(azimuth=-45, dist_hor=3.0)
#Lower righ point of map
pLR=p2.offset(azimuth=135, dist_hor=3.0)

#Now load topodata for the area spanned by the just defined corners
topo=pTL.get_topo_data(pLR)[0]
print topo.resolution
#topo.increase_resolution(0.03) #increase topo resolution grid to 30m

p1.set_topo_data(topo)
p2.set_topo_data(topo)

dz=abs(np.nanmax(topo.data)-np.nanmin(topo.data))
m=topo.plot_3d()
p1.plot(m, addName=True, dzText=dz*0.1)
p2.plot(m, addName=True, dzText=dz*0.1)

#Now create a setup from this stuff
stp=gnp.base.GeoSetup("Altiplano", localTopoPath=topoPath,points=[p1,p2],vectors=[vc])

#==============================================================================
# 
# ep1=p1.get_elevation_profile(azimuth=vc.azimuth, dist_hor=vc.magnitude*1.4)
# fig, ax=plt.subplots(2,1, figsize=(18,8))
# ep1.plot(ax=ax[0])
# 
#==============================================================================


