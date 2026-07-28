function slipq, bx0, by0, bz0, bx1, by1, bz1, xa=xa, ya=ya, za=za, spherical=spherical, $
xreg=xreg, yreg=yreg, factor=factor, delta=delta, lon_delta=lon_delta, lat_delta=lat_delta, $
qsl0=qsl0, preview=preview, fname=fname
;------------------------------------------------------------
fastqsl, bx0, by0, bz0, xa=xa, ya=ya, za=za, spherical=spherical, $
xreg=xreg, yreg=yreg, factor=factor, delta=delta, lon_delta=lon_delta, lat_delta=lat_delta, $
/rf, /seed, /b, qsl=qsl0, odir=odir, preview=preview, fname=fname+'_t0', /silent

nq1=qsl0.dim[0]
nq2=qsl0.dim[1]

mapt0=fltarr(3, nq1, nq2)

for j=0, nq2-1 do begin
for i=0, nq1-1 do begin
    case qsl0.sign2d[i,j] of
	-1: mapt0[*,i,j]=qsl0.rfs[*,i,j]
     0: mapt0[*,i,j]=qsl0.seed[*,i,j]
     1: mapt0[*,i,j]=qsl0.rfe[*,i,j]
    endcase
endfor
endfor
;------------------------------------------------------------
if N_PARAMS() eq 2 then bx1=by0
fastqsl, bx1, by1, bz1, xa=xa, ya=ya, za=za, spherical=spherical, seed=mapt0, $
/rf, /b, /targetB, qsl=qsl1, /silent

mapt0mapt1=fltarr(3, nq1, nq2)
Bn_target_t1=fltarr(nq1, nq2)

for j=0, nq2-1 do begin
for i=0, nq1-1 do begin
    if (qsl1.b[2,i,j] lt 0.) then begin
        mapt0mapt1[*,i,j]=qsl1.rfs[*,i,j]
        Bn_target_t1[i,j]=qsl1.Bs[2,i,j]
    endif else if (qsl1.b[2,i,j] gt 0.) then begin
        mapt0mapt1[*,i,j]=qsl1.rfe[*,i,j]
        Bn_target_t1[i,j]=qsl1.Be[2,i,j]
    endif else if (qsl1.b[2,i,j] eq 0.) then begin
        mapt0mapt1[*,i,j]=qsl1.seed[*,i,j]
        Bn_target_t1[i,j]=qsl1.B[2,i,j]
    endif
endfor
endfor
;------------------------------------------------------------
slipq01=fltarr(nq1, nq2) 

if (spherical) then begin
    g_launch_t0=cos(reform(qsl0.seed[1,0,*]))
    g_target_t1=cos(reform(mapt0mapt1[1,*,*]))
endif else begin
    g_launch_t0=fltarr(nq2)+1.
    g_target_t1=fltarr(nq1, nq2)+1.
endelse

if spherical then area=4*qsl0.lon_delta*qsl0.lat_delta else area=4*qsl0.delta^2

for j=1, nq2-2 do begin
for i=1, nq1-2 do begin
if qsl0.rboundary[i,j-1] eq 11 and $
   qsl0.rboundary[i-1,j] eq 11 and $
   qsl0.rboundary[i,j]   eq 11 and $
   qsl0.rboundary[i+1,j] eq 11 and $
   qsl0.rboundary[i,j+1] eq 11 and $
   qsl1.rboundary[i,j-1] eq 11 and $
   qsl1.rboundary[i-1,j] eq 11 and $
   qsl1.rboundary[i,j]   eq 11 and $
   qsl1.rboundary[i+1,j] eq 11 and $
   qsl1.rboundary[i,j+1] eq 11 then begin
    slipq01[i,j]=(((mapt0mapt1[0,i+1,j]-mapt0mapt1[0,i-1,j])*g_target_t1[i,j]/g_launch_t0[j])^2.+ $
                  ((mapt0mapt1[1,i+1,j]-mapt0mapt1[1,i-1,j])/g_launch_t0[j])^2.                 + $
                  ((mapt0mapt1[0,i,j+1]-mapt0mapt1[0,i,j-1])*g_target_t1[i,j])^2.               + $
                   (mapt0mapt1[1,i,j+1]-mapt0mapt1[1,i,j-1])^2.) $
                 /area*abs(Bn_target_t1[i,j]/qsl0.b[2,i,j])
endif else slipq01[i,j]= !values.f_nan
endfor
endfor

abnormal= where(slipq01 lt 2.)
if abnormal[0] ne -1 then slipq01[abnormal]=2.

if preview then begin
    write_png, odir+fname+'.png', bytscl(alog10(slipq01), min=1., max=5.,/nan)
    print, odir+fname+'.png'
endif

return, slipq01

end
