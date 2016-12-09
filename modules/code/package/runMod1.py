
import FRF2Shape
from numpy import *

# -- To be completed at page load --
## Read in FRF data
dat_fname = 'mod1_FRF.csv'
FRF, w = FRF2Shape.importFRF(dat_fname)
## Read in dof location coordinates
coord_fname = 'DOF_coord.csv'
DOF_coord = FRF2Shape.importCSV(coord_fname)
# boundary coordinates
boundx = array([0,576])
boundy = array([0,72,144,216,288])
boundx, boundy = meshgrid(boundx,boundy)
boundx = boundx.flatten()
boundy = boundy.flatten()
bcoords = hstack((boundx[...,None],boundy[...,None],zeros(boundx.shape[0])[...,None]))
## Compile input parameters
# input dof locations
ins = arange(6,11)
# output dof locations
outs = arange(6,22)
# other future capabilities


## Compute Shapes from FRF subset
# use predetermined poles
poles = array([5.54, 7.84, 11.67, 15.34, 19.47])
# Run getCMIF
CMIF, shapes, FRFsub = FRF2Shape.getCMIF(FRF,w,pole=poles,inDOF=ins,outDOF=outs)

## Plot results
# Plot FRF

# Plot CMIF
# Plot curve fit of ShapeArray (with DOF coords)
for eachshape in range(0,poles.shape[0]):
    FRF2Shape.modeInterp(DOF_coord[outs,:],shapes[:,eachshape],bcoords,75)
