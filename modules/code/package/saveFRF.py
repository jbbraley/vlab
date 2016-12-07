# resave txt file with complex numbers

from numpy import *
filedat = loadtxt('mod1_FRF.csv',delimiter=',',skiprows=1)
w = filedat[:,0]
numcolumns = filedat.shape[1]
IOsize = rint((numcolumns-1)/2)
realdat = filedat[:,1:IOsize.astype(int)+1]
imagdat = filedat[:,IOsize.astype(int)+1:]*1j
FRF = realdat+imagdat

FRF.reshape()
