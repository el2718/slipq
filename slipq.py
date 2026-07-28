import numpy as np
from fastqsl import fastqsl
import matplotlib.pyplot as plt

def slipq(bx0, by0, bz0, bx1, by1, bz1, xa=None, ya=None, za=None, spherical=False, \
xreg=None, yreg=None, factor=4, delta=None, lon_delta=None, lat_delta=None, preview=False, fname='slipq'):
# ------------------------------------------------------------
    qsl0=fastqsl(bx0, by0, bz0, xa=xa, ya=ya, za=za, spherical=spherical, \
    xreg=xreg, yreg=yreg, factor=factor, delta=delta, lon_delta=lon_delta, lat_delta=lat_delta, \
    rF_out=True, seed=True, B_out=True, silent=True, fname=fname+'_t0', preview=preview)

    nq2, nq1=qsl0['dim']

    mapt0=np.zeros((nq2, nq1, 3),'f4')
    for j in range(nq2):
        for i in range(nq1):
            if   qsl0['sign2d'][j,i] == -1: mapt0[j,i,:]=qsl0['rFs'][j,i,:]
            elif qsl0['sign2d'][j,i] ==  0: mapt0[j,i,:]=qsl0['seed'][j,i,:]
            elif qsl0['sign2d'][j,i] ==  1: mapt0[j,i,:]=qsl0['rFe'][j,i,:]
# ------------------------------------------------------------
    if by0 is not None and bz0 is None : bx1=by0
    qsl1=fastqsl(bx1, by1, bz1, xa=xa, ya=ya, za=za, spherical=spherical, seed=mapt0, \
             rF_out=True, B_out=True, targetB_out=True, silent=True)
    
    mapt0mapt1=np.zeros((nq2, nq1, 3),'f4')
    Bn_target_t1=np.zeros((nq2, nq1),'f4')
    for j in range(nq2):
        for i in range(nq1):
            if   qsl1['B'][j,i,2] < 0.:
                mapt0mapt1[j,i,:]=qsl1['rFs'][j,i,:]
                Bn_target_t1[j,i]=qsl1['Bs'][j,i,2]
            elif qsl1['B'][j,i,2] == 0.:
                mapt0mapt1[j,i,:]=qsl1['seed'][j,i,:]
                Bn_target_t1[j,i]=qsl1['B'][j,i,2]
            elif qsl1['B'][j,i,2] > 0:
                mapt0mapt1[j,i,:]=qsl1['rFe'][j,i,:]
                Bn_target_t1[j,i]=qsl1['Be'][j,i,2]
# ------------------------------------------------------------
    slipq01=np.zeros((nq2, nq1),'f4')

    if spherical:
        g_launch_t0=np.cos(qsl0.seed[:,0,1])
        g_target_t1=np.cos(mapt0mapt1[:,:,1])
    else:
        g_launch_t0=np.zeros((nq2),'f4')+1.
        g_target_t1=np.zeros((nq2, nq1),'f4')+1.

    if spherical:
        area=4*qsl0['lon_delta']*qsl0['lat_delta']
    else:
        area=4*qsl0['delta']**2.

    for j in range(1,nq2-1):
        for i in range(1,nq1-1):
            if qsl0['rboundary'][j-1,i] == 11 and \
               qsl0['rboundary'][j,i-1] == 11 and \
               qsl0['rboundary'][j,i]   == 11 and \
               qsl0['rboundary'][j,i+1] == 11 and \
               qsl0['rboundary'][j+1,i] == 11 and \
               qsl1['rboundary'][j-1,i] == 11 and \
               qsl1['rboundary'][j,i-1] == 11 and \
               qsl1['rboundary'][j,i]   == 11 and \
               qsl1['rboundary'][j,i+1] == 11 and \
               qsl1['rboundary'][j+1,i] == 11 :
                slipq01[j,i]=(((mapt0mapt1[j,i+1,0]-mapt0mapt1[j,i-1,0])*g_target_t1[j,i]/g_launch_t0[j])**2+ \
                              ((mapt0mapt1[j,i+1,1]-mapt0mapt1[j,i-1,1])/g_launch_t0[j])**2                 + \
                              ((mapt0mapt1[j+1,i,0]-mapt0mapt1[j-1,i,0])*g_target_t1[j,i])**2               + \
                               (mapt0mapt1[j+1,i,1]-mapt0mapt1[j-1,i,1])**2) \
                             /area*np.abs(Bn_target_t1[j,i]/qsl0['B'][j,i,2])
            else: slipq01[j,i]=np.nan
    

    slipq01[slipq01 <=2.]=2.

    slipq_tmp=slipq01.copy()
    slipq_tmp[np.isnan(slipq_tmp)]=2.
    slipq_tmp[np.isinf(slipq_tmp)]=2.

    if preview: 
        odir='fastqsl/'
        plt.imsave(odir+fname+'.png', np.log10(slipq_tmp), vmin=1., vmax=5., origin='lower', cmap='gray')
        print(odir+fname+'_magnetogram.png')

    return slipq01, qsl0
