dx=0.2
dy=0.2
dz=0.2
x0=-15.
y0=-15.
z0=0.
nx=long((-x0)/dx)*2+1
ny=long((-y0)/dy)*2+1
nz=long(30/dx)*2+1

xa=findgen(nx)*dx+x0
ya=findgen(ny)*dy+y0
za=findgen(nz)*dz+z0

By=fltarr(nx,ny,nz)+0.2

Bz=fltarr(nx,ny,nz)
for i=0, nx-1 do Bz[i,*,*]=xa[i]

Ly=5.
Lz=5.
dummy=fltarr(nx,ny,nz)
for k=0, nz-1 do begin
    z=k*dz+z0
for j=0, ny-1 do begin
    y=j*dy+y0
    dummy[*,j,k]=(1.-(z/Lz)^2.)/((1+(y/Ly)^2.)*(1+(z/Lz)^2.)^2.)
endfor
endfor

Bx8 =-1.- 8*dummy
Bx13=-1.-13*dummy

qsf=slipq(bx8, by, bz, bx13, by, bz, xa=xa, ya=ya, za=za, delta=0.05, $
xreg=[-10,10], yreg=[-10,10], fname='Btitov2009_qsf', /preview)
qsb=slipq(bx13, by, bz, bx8, by, bz, xa=xa, ya=ya, za=za, delta=0.05, $
xreg=[-10,10], yreg=[-10,10], fname='Btitov2009_qsb', /preview)

end
