function func, x
common params, hsmall, qsmall, adia, mm, kx, molec 

ii = dcomplex(0d0, 1d0)

zeta = x(0)

mu  = molec
eps = 0.022*zeta^(2./7.) 
smallq= qsmall 
gmma = adia 
khat = kx
mmode = mm 

bcool = (4.4d5/mu/(gmma-1d0))*zeta^(-57./14.)
bcool*= 8.3d-9*zeta^(33./7.) + 1d0/khat^2 
sig   = x(1) 

ksq  = khat^2

bigA = 1d0 + ii*eps*smallq*khat

c0 = mmode*(mmode+1d0)*bigA^2

c1 = ii*bcool*( (1d0-gmma) - ksq*(1d0+2d0*mmode)^2*(gmma-1d0) + 4d0*ii*eps*smallq*khat*mmode*(mmode+1d0)*(gmma-1d0) - 2d0*bigA^2*gmma*mmode*(mmode+1d0) )

c2 = (ksq+1d0)*bigA + bcool^2*( (1d0-gmma) + gmma*(1d0-gmma)*( ksq*(1d0+2d0*mmode)^2 - 4d0*ii*eps*smallq*khat*mmode*(mmode+1d0) ) - gmma^2*bigA^2*mmode*(mmode+1d0) )

c3 = bcool*( -3d0*ii - 2d0*ii*ksq + eps*khat*smallq + gmma*(ii + smallq*eps*khat*(1d0+2d0*ksq)) )

c4 = -(1d0+ksq)^2 + bcool^2*(1d0+gmma*ksq)*(-2d0 + gmma - ii*eps*smallq*khat*gmma)

c5 = 2d0*ii*bcool*(1d0+ksq)*(1d0+gmma*ksq)

c6 = bcool^2*(1d0+gmma*ksq)^2

coeff = [c0, c1, c2, c3, c4, c5, c6]

res = dcomplex(0d0, 0d0)

if (mm gt 0d0) then begin

for i=0, 6 do begin
res += coeff(i)*sig^(i+0d0)
endfor

endif else begin

for i=1, 6 do begin
res += coeff(i)*sig^(i-1d0) 
endfor

endelse

return, [real_part(res),imaginary(res)]

end

pro bcrit_theory_mmsn, krange=krange, legend=legend, mmode=mmode, yrange=yrange, xrange=xrange, label=label

common params, hsmall, qsmall, adia, mm, kx, molec 

gmma = 1.4
mu   = 2.0 
smallq = -3./7. 


qsmall = smallq
adia   = gmma 
molec  = mu 


!p.font = 0

ii = dcomplex(0d0, 1d0)

if not keyword_set(mmode) then mmode=0d0
title = 'M='+string(mmode,format='(I2)')
mm = mmode 

if not keyword_set(krange) then krange=[10,100]+0d0 
kmax = krange(1)+0d0
kmin = krange(0)+0d0
if not keyword_set(kswitch) then kswitch = kmax 
nk = 81920
dlogk = alog10(kmax/kmin)/(nk-1d0)
kaxis = krange(0) + (krange(1)-krange(0))*dindgen(nk)/(nk-1d0)

bvals = dblarr(nk)
freqs = dblarr(nk)

bvals1 = dblarr(nk)
freqs1 = dblarr(nk)

xvec  = [0.74826817,-0.00010001318];kx=1d4
;xvec = [2.2055799,-0.0010002441];kx=1000
;xvec=[5.0,-0.1];kx=100 

for j=0, nk-1 do begin
;for j=0, 0 do begin 

kx = kaxis(j)

Result = NEWTON(xvec, 'func', /DOUBLE, TOLF=1d-9, tolx=1d-9, itmax=50000)
bvals(j) = result(0)
freqs(j) = result(1)
xvec = result 


;Result1 = NEWTON(xvec1, 'func', /DOUBLE, TOLF=1d-9, tolx=1d-9, itmax=50000)
;bvals1(j) = result1(0)
;freqs1(j) = result1(1)
;xvec1 = result1
endfor
;bvals is the critical distance for each kx. calculate from this the critical cooling time

zeta = bvals 
bcrit_arr = (4.4d5/mu/(gmma-1d0))*zeta^(-57./14.)
bcrit_arr*= 8.3d-9*zeta^(33./7.) + 1d0/kaxis^2

  set_plot, 'ps'
  device, filename='bcrit_theory.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

  plot, kaxis, bvals,xmargin=[8.2,1.8],ymargin=[3.5,2] $
        ,ytitle=textoidl('\zeta'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title,/xlog ,/ylog

  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, 2 do begin
;        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
        xyouts, x0, legend(j+1),textoidl(label(j)),charsize=1.5
     endfor
  endif
  device,/close


nrad = 81920
rmin=min(bvals)
rmax=max(bvals)
zeta = 10.^(alog10(rmin) + alog10(rmax/rmin)*dindgen(nrad)/(nrad-1d0))

bcrit = zeta 
for i=0, nrad-1 do begin
temp = min(abs(zeta(i)-bvals),kgrid)

bcrit(i) = (4.4d5/mu/(gmma-1d0))*zeta(i)^(-57./14.)
bcrit(i)*= 8.3d-9*zeta(i)^(33./7.) + 1d0/kaxis(kgrid)^2
endfor

kx = [1d0, 1d1, 1d2]

bcool = dblarr(3, nrad)

for k=0, 2 do begin

bcool(k,*)  = 4.4d5/mu/(gmma-1d0)*zeta^(-57./14.)
bcool(k,*) *= 8.3d-9*zeta^(33./7.) + 1d0/kx(k)^2
endfor

  set_plot, 'ps'
  device, filename='bcrit_mmsn.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

  plot, zeta, bcrit,xmargin=[7.8,2.2],ymargin=[3.5,2] $
        ,ytitle=textoidl('t_c\Omega_k'),xtitle=textoidl('\zeta=r/AU') $
        ,charsize=1.5, thick=4, xrange=[1,1d3],xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title,/xlog ,/ylog

kx = [1d0, 1d1, 1d2]

bcool = dblarr(3, nrad)

rmin=1d0
rmax=1d3

zeta = 10.^(alog10(rmin) + alog10(rmax/rmin)*dindgen(nrad)/(nrad-1d0))

for k=0, 2 do begin

bcool(k,*)  = 4.4d5/mu/(gmma-1d0)*zeta^(-57./14.)
bcool(k,*) *= 8.3d-9*zeta^(33./7.) + 1d0/kx(k)^2
endfor

  for k=0, 2 do begin
  oplot,zeta, bcool(k,*), thick=4, linestyle=k+1
  endfor

  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)

     for j=0, n_elements(label)-1 do begin
        ynew = 10.^(alog10(y0) - dy*j)
        oplot, [x0,x1], [ynew,ynew], thick=4, linestyle=j
        xyouts, x1, ynew,textoidl(label(j)),charsize=1.5
     endfor
  endif
  device,/close


end
