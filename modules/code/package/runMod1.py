
import FRF2Shape
from numpy import *

# -- To be completed at page load --
## Read in FRF data
dat_fname = 'mod1_FRF.csv'
FRF, w = FRF2Shape.importFRF(dat_fname)
## Read in dof location coordinates
## Compile input parameters
# input dof locations
ins = arange(5,11)
# output dof locations
outs = arange(5,21)
# other future capabilities


## Compute Shapes from FRF subset
# use predetermined poles
poles = array([5.54, 7.84, 11.67, 15.34, 19.47])
# Run getCMIF
CMIF, shapes = FRF2Shape.getCMIF(FRF,w,pole=poles,inDOF=ins,outDOF=outs)

print(CMIF.shape)
print(shapes.shape)

## Plot results
# Plot FRF
# Plot CMIF
# Plot curve fit of ShapeArray (with DOF coords)
