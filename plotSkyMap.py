from matplotlib.font_manager import FontProperties
from SkyMap import *
from matplotlib import cm

s = SkyMap()


mapfile ="Results/mapIC86_ndof2_2DegCoarseBins_RealData.root"
#mapfile ="ScrambledSkyMaps/mapIC86_989_fixedsigma_ndof2_2DegCoarseBins_ns0_width0.root"

#my_cmap = matplotlib.cm.binary

proj = "hammer"
s.AddHistogramFile(mapfile)
#s.AddGalPlane("resources/GalPlane3600.coords", "blue", "--", 2)
#s.evxcoord = [260,275]

#s.evycoord = [50,60]
s.hotspot = True
s.set_latitude_grid(4)
s.SetGridColor("blue")
#s.cmap = my_cmap
s.Plot(proj)


savefig("RealMapFine_BigSpot")

print "It's Showtime!"

#show()



