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
import scipy.linalg

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
    FRFsub = FRF

    if outDOF is None:
        no = FRF.shape[0]
    else:
        FRFsub = FRFsub[outDOF-1,:,:]
        no = outDOF.shape[0]

    if inDOF is None:
        ni = FRF.shape[1]
    else:
        FRFsub = FRFsub[:,inDOF-1,:]
        ni = inDOF.shape[0]

    if FRFsub.shape[1]>FRFsub.shape[0]:
        FRFsub = transpose(FRFsub,(1,0,2))
        no = ni
        ni = FRFsub.shape[1]

    # Compute CMIF
    freq_lines = FRFsub.shape[2]
    # preallocate
    ns = FRFsub.shape[2]
    uu = zeros((no,ni,ns),dtype=complex64)
    vv = zeros((ni,ni,ns),dtype=complex64)
    ss = zeros((ns,ni),dtype=complex64)
    for ii in range(0,freq_lines-1):
        u, s, v = linalg.svd(asmatrix(FRFsub[:,:,ii]),full_matrices=False, compute_uv=True)
        uu[:,:,ii] = u
        vv[:,:,ii] = v.T
        ss[ii,:] = s

    CMIF = ss

    if pole is not None:
        # use predefined pole locations
        ns_ind = zeros(pole.shape[0])
        for ii in range(0,pole.shape[0]):
            ns_ind[ii] = abs(w-pole[ii]).argmin()
        ##Compute mode shapes at pole locations indices
        print(w[ns_ind.astype(int)])
        shapes = uu[:,:,ns_ind.astype(int)]

        return CMIF, shapes
    else:
        return CMIF


def main():
    print('{} called'.format(__name__))

if __name__ == "__main__":
    main()
