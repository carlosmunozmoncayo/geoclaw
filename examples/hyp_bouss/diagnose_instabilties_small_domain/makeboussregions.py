from pylab import *
import os,sys
from clawpack.amrclaw import region_tools
from clawpack.geoclaw import topotools, dtopotools, kmltools, fgmax_tools, marching_front

#Read topo file
CLAW = os.environ['CLAW']
scratch_dir = os.path.join(CLAW, 'geoclaw', 'scratch')
topo_fname = 'etopo10min120W60W60S0S.asc'
topo = topotools.Topography(os.path.join(scratch_dir,topo_fname), topo_type=2)

# Specify a RuledRectangle the our flagregion should lie in:
rrect = region_tools.RuledRectangle()
rrect.ixy = 'y'  # s = latitude
rrect.s = linspace(topo.extent[-2],topo.extent[-1],4)
rrect.lower =  -85*ones(rrect.s.shape)
rrect.upper = topo.extent[1]*ones(rrect.s.shape)
rrect.method = 1
xr,yr = rrect.vertices()

# Start with a mask defined by the ruled rectangle `rrect` defined above:
mask_out = rrect.mask_outside(topo.X, topo.Y)

# select onshore points within 2 grip points of shore:
pts_chosen_Zabove0 = marching_front.select_by_flooding(topo.Z, mask=mask_out, 
                                                       prev_pts_chosen=None, 
                                                       Z1=0, Z2=1e6, max_iters=None)
# select offshore points down to 1000 m depth:
pts_chosen_Zbelow0 = marching_front.select_by_flooding(topo.Z, mask=None, 
                                                       prev_pts_chosen=None, 
                                                       Z1=0, Z2=-100., max_iters=None)
# buffer offshore points with another 10 grid cells:
pts_chosen_Zbelow0 = marching_front.select_by_flooding(topo.Z, mask=None, 
                                                       prev_pts_chosen=pts_chosen_Zbelow0, 
                                                       Z1=0, Z2=-5000., max_iters=5)

# Take the intersection of the two sets of points selected above:
nearshore_pts = where(pts_chosen_Zabove0+pts_chosen_Zbelow0 == 2, 1, 0)

rr2 = region_tools.ruledrectangle_covering_selected_points(topo.X, topo.Y,
                                                          nearshore_pts, ixy='y', method=1,
                                                          verbose=True)
#Save rr2 to a file
rr2.write('rr_nearshore_continental_SouthAmerica.data')


# Specify a RuledRectangle the our flagregion should lie in:
rrect = region_tools.RuledRectangle()
rrect.ixy = 'y'  # s = latitude
rrect.s = linspace(topo.extent[-2],topo.extent[-1],4)
rrect.lower =  topo.extent[0]*ones(rrect.s.shape)
rrect.upper = -85*ones(rrect.s.shape)
rrect.method = 1

# Start with a mask defined by the ruled rectangle `rrect` defined above:
mask_out = rrect.mask_outside(topo.X, topo.Y)

# select onshore points within 2 grip points of shore:
pts_chosen_Zabove0 = marching_front.select_by_flooding(topo.Z, mask=mask_out, 
                                                       prev_pts_chosen=None, 
                                                       Z1=0, Z2=1e6, max_iters=None)
# select offshore points down to 1000 m depth:
pts_chosen_Zbelow0 = marching_front.select_by_flooding(topo.Z, mask=None, 
                                                       prev_pts_chosen=None, 
                                                       Z1=0, Z2=-100., max_iters=None)
# buffer offshore points with another 10 grid cells:
pts_chosen_Zbelow0 = marching_front.select_by_flooding(topo.Z, mask=None, 
                                                       prev_pts_chosen=pts_chosen_Zbelow0, 
                                                       Z1=0, Z2=-5000., max_iters=5)

# Take the intersection of the two sets of points selected above:
nearshore_pts = where(pts_chosen_Zabove0+pts_chosen_Zbelow0 == 2, 1, 0)

rr3 = region_tools.ruledrectangle_covering_selected_points(topo.X, topo.Y,
                                                          nearshore_pts, ixy='y', method=1,
                                                          verbose=True)

#Save rr3 to a file
rr3.write('rr_nearshore_ocean_SouthAmerica.data')