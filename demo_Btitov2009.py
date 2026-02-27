import numpy as np
from slipq import slipq

nx=151
ny=151
nz=301

xa=np.linspace(-15., 15., nx)
ya=np.linspace(-15., 15., ny)
za=np.linspace(  0., 60., nz)

by=np.zeros((nz,ny,nx),'f4')+0.2
bz=np.zeros((nz,ny,nx),'f4')
for i in range(nx): bz[:,:,i]=xa[i]

Ly=5.
Lz=5.
dummy=np.zeros((nz,ny,nx),'f4')
for k in range(nz):
    for j in range(ny):
        dummy[k,j,:]=(1.-(za[k]/Lz)**2)/((1+(ya[j]/Ly)**2)*(1+(za[k]/Lz)**2)**2)

bx8 =-1.- 8*dummy
bx13=-1.-13*dummy

qsf, qsl8 =slipq(bx8, by, bz, bx13, by, bz, xa=xa, ya=ya, za=za, delta=0.05, \
xreg=[-10,10], yreg=[-10,10], fname='Btitov2009_qsf', preview=True)

qsb, qsl13=slipq(bx13, by, bz, bx8, by, bz, xa=xa, ya=ya, za=za, delta=0.05, \
xreg=[-10,10], yreg=[-10,10], fname='Btitov2009_qsb', preview=True)