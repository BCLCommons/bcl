from pymol.cgo import *
from pymol import cmd
momentssse_a=[]
cur_mom=[ ALPHA,1,CYLINDER,34.591,-17.827,11.134,33.3761,-19.4017,11.7637,0.1,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,0.5,CYLINDER,34.591,-17.827,11.134,33.3761,-19.4017,11.7637,0.2,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,1,CYLINDER,36.949,-17.244,14.022,35.7341,-18.8187,14.6517,0.1,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,0.5,CYLINDER,36.949,-17.244,14.022,35.7341,-18.8187,14.6517,0.2,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,1,CYLINDER,39.75,-18.933,12.138,38.5351,-20.5077,12.7677,0.1,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,0.5,CYLINDER,39.75,-18.933,12.138,38.5351,-20.5077,12.7677,0.2,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,1,CYLINDER,39.089,-16.764,9.12,37.8741,-18.3387,9.74975,0.1,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,0.5,CYLINDER,39.089,-16.764,9.12,37.8741,-18.3387,9.74975,0.2,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,1,CYLINDER,39.189,-13.663,11.268,37.9741,-15.2377,11.8977,0.1,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,0.5,CYLINDER,39.189,-13.663,11.268,37.9741,-15.2377,11.8977,0.2,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,1,CYLINDER,42.486,-14.743,12.756,41.2711,-16.3177,13.3857,0.1,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,0.5,CYLINDER,42.486,-14.743,12.756,41.2711,-16.3177,13.3857,0.2,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,1,CYLINDER,43.912,-15.304,9.306,42.6971,-16.8787,9.93575,0.1,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,0.5,CYLINDER,43.912,-15.304,9.306,42.6971,-16.8787,9.93575,0.2,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,1,CYLINDER,42.791,-11.864,8.232,41.5761,-13.4387,8.86175,0.1,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,0.5,CYLINDER,42.791,-11.864,8.232,41.5761,-13.4387,8.86175,0.2,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,1,CYLINDER,44.425,-10.332,11.271,43.2101,-11.9067,11.9007,0.1,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cur_mom=[ ALPHA,0.5,CYLINDER,44.425,-10.332,11.271,43.2101,-11.9067,11.9007,0.2,0,1,1,0,1,1]
momentssse_a.extend(cur_mom)

cmd.load_cgo( momentssse_a, 'momentssse_a')
calculated_exposure_centeredsse_a=[]
cur_mom=[ ALPHA,1,CYLINDER,40.6195,-15.3889,10.9666,39.4045,-16.9636,11.5963,0.1,0,1,1,0,1,1]
calculated_exposure_centeredsse_a.extend(cur_mom)

cmd.load_cgo( calculated_exposure_centeredsse_a, 'calculated_exposure_centeredsse_a')
experimental_access_centeredsse_a=[]
cur_mom=[ ALPHA,0.5,CYLINDER,40.6195,-15.3889,10.9666,39.4045,-16.9636,11.5963,0.2,0,1,1,0,1,1]
experimental_access_centeredsse_a.extend(cur_mom)

cmd.load_cgo( experimental_access_centeredsse_a, 'experimental_access_centeredsse_a')
calculated_exposure_xysse_a=[]
cur_mom=[ ALPHA,1,CYLINDER,40.6195,-15.3889,10.9666,41.2683,-16.0579,11.3293,0.1,0,1,1,0,1,1]
calculated_exposure_xysse_a.extend(cur_mom)

cmd.load_cgo( calculated_exposure_xysse_a, 'calculated_exposure_xysse_a')
experimental_access_xysse_a=[]
cur_mom=[ ALPHA,0.5,CYLINDER,40.6195,-15.3889,10.9666,41.2683,-16.0579,11.3293,0.2,0,1,1,0,1,1]
experimental_access_xysse_a.extend(cur_mom)

cmd.load_cgo( experimental_access_xysse_a, 'experimental_access_xysse_a')
xaxissse_a=[CYLINDER,40.6195,-15.3889,10.9666,39.996,-14.8105,10.4404,0.1,1,1,0,1,1,0]

cmd.load_cgo( xaxissse_a, 'x-axissse_a')
yaxissse_a=[CYLINDER,40.6195,-15.3889,10.9666,40.425,-14.8518,11.7874,0.1,1,1,0,1,1,0]

cmd.load_cgo( yaxissse_a, 'y-axissse_a')
zaxissse_a=[CYLINDER,40.6195,-15.3889,10.9666,41.3768,-14.7749,10.7442,0.1,1,1,0,1,1,0]

cmd.load_cgo( zaxissse_a, 'z-axissse_a')
