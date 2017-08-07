"""
Computes mode shapes from FRF matrix
 Input:
   FRF - nxnxm matrix, where n = number of DOF, m = frequency bins
   iDOF - array of integers that identify which DOF are impacted
   oDOF - array of integers that identify which DOF are sampled (i.e. sensor locations)
   fn - Natural frequencies for which to compute shapes

 Output:
   shapes = nxm array for which each element is the relative amplitude of each DOF at specified of m modes
   FRFsub = nxm array that is a subset of input FRF where n = number of iDOF, and m = number of oDOF
   CMIF = singular values at every frequency line
"""

from numpy import *
from scipy import linalg, interpolate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def importCSV(filename):
    # Import from text file (ignore header row)
    dat = loadtxt(filename,delimiter=',',skiprows=1)
    return dat


def importFRF(filename):
    # Import from text file (ignore header row)
    filedat = loadtxt(filename,delimiter=',',skiprows=1)
    # First column in spectral lines
    w = filedat[:,0]
    # Recombine real and imaginary components
    numcolumns = filedat.shape[1]
    IOsize = rint((numcolumns-1)/2)
    realdat = filedat[:,1:IOsize.astype(int)+1]
    imagdat = filedat[:,IOsize.astype(int)+1:]*1j
    frf = realdat+imagdat
    # Reshape such that noxnixns
    FRF = frf[..., newaxis]
    FRF = transpose(FRF,(2,1,0)).reshape(25,25,-1)
    return FRF, w

def getCMIF(FRF,w,pole=None,outDOF=None,inDOF=None):
    """
    Performs a singular value decomposition of the chosen subset of the FRF at every frequency line
    If no input or output degrees of freedom are specified, the svd is performed on the full FRF
    If poles are specified the shapes at those pole are returned, otherwise all the left singular vectors are returned
    """

    # Resample FRF to form subset
    if outDOF is None:
        no = FRF.shape[0]
        outs = arange(0,no-1)
    else:
        outs = outDOF-1
        no = outDOF.shape[0]

    if inDOF is None:
        ni = FRF.shape[1]
        ins = arange(0,ni-1)
    else:
        ins = inDOF-1
        ni = inDOF.shape[0]

    FRFsub = FRF[outs,:,:][:,ins,:]

    if FRFsub.shape[1]>=FRFsub.shape[0]:
        inswitchout = 1
        ni, no = no, ni
    else:
        inswitchout=0

    # Compute CMIF
    # preallocate
    ns = FRFsub.shape[2]
    uu = zeros((no,ni,ns),dtype=complex64)
    vv = zeros((ni,ni,ns),dtype=complex64)
    ss = zeros((ns,ni),dtype=complex64)
    for ii in range(0,ns-1):
        if inswitchout==1:
            frf = asmatrix(FRFsub[:,:,ii]).H
        else:
            frf = asmatrix(FRFsub[:,:,ii])

        u, s, v = linalg.svd(frf,full_matrices=False, compute_uv=True)
        uu[:,:,ii] = u
        vv[:,:,ii] = v
        ss[ii,:] = s

    CMIF = ss

    if pole is not None:

        if inswitchout==1:
            allshapes = transpose(real(vv),(1,0,2))
        else:
            allshapes = imag(uu)
        # use predefined pole locations
        ns_ind = zeros(pole.shape[0])
        shapes = zeros((allshapes.shape[0],pole.shape[0]))
        for ii in range(0,pole.shape[0]):
            # find corresponding spectral line
            ns_ind[ii] = abs(w-pole[ii]*2*pi).argmin()
            # pull shapes from singular vector
            shapes[:,ii] = allshapes[:,ss[ns_ind[ii].astype(int),:].argmax(),ns_ind[ii].astype(int)]

        return CMIF, shapes, FRFsub
    else:
        return CMIF, FRFsub

#def plotFRF(FRF):

def modeInterp(coords,z,bounds,scale, saveloc, savename):
    # concat geometry
    xTot = r_[coords[:,0],bounds[:,0]]
    yTot = r_[coords[:,1],bounds[:,1]]
    zTot = r_[z,bounds[:,2]]

    # define grid resolution
    xres = 49
    yres = 49

    # interpolate
    xv = linspace(xTot.min(), xTot.max(), xres)
    yv = linspace(yTot.min(), yTot.max(), yres)
    xInter, yInter = meshgrid(xv, yv)
    tck = interpolate.bisplrep(xTot,yTot,zTot)
    zInter = interpolate.bisplev(xInter[0,:],yInter[:,0],tck)

    # plot mesh surface
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    rainbow = cm = plt.get_cmap('gist_rainbow')
    surf = ax.plot_surface(xInter, yInter, zInter*scale, rstride=1,
                cstride=1, cmap=rainbow, linewidth=0, antialiased=False)

    plt.tight_layout()
    fig.subplots_adjust(bottom=0.1)

    plt.savefig(saveloc + '/' + savename) # bbox_inches='tight')
    #plt.show()

    plt.close(fig)



def main():
    print('{} called'.format(__name__))

if __name__ == "__main__":
    main()
