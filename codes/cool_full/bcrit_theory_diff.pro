function func0, x
common params, hsmall, qsmall, adia, mm, kx

eps = hsmall
smallq= qsmall
gmma = adia
khat = kx
mmode = mm

b0 = x(0)
om0= x(1)

res1 = gmma*b0^2*(1d0-gmma) + 2d0*eps*smallq*gmma*b0*om0 - om0^2 + b0^2*gmma^2*om0^4
res2 = b0*(1d0-gmma) + eps*smallq*om0 - eps*smallq*gmma^2*b0^2*om0^3

return, [res1,res2]
end

function func, x
common params, hsmall, qsmall, adia, mm, kx 

ii = dcomplex(0d0, 1d0)

eps = hsmall
smallq= qsmall 
gmma = adia 
khat = kx
mmode = mm 

bcool = x(0)/khat^2
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

pro bcrit_theory_diff, eps=eps, smallq=smallq, gmma=gmma, krange=krange, legend=legend, mmode=mmode, yrange=yrange, xrange=xrange, label=label

common params, hsmall, qsmall, adia, mm, kx 

hsmall = eps 
qsmall = smallq
adia    = gmma 

!p.font = 0

ii = dcomplex(0d0, 1d0)

if not keyword_set(mmode) then mmode=0d0
title = 'M='+string(mmode,format='(I2)')
mm = mmode 

if not keyword_set(krange) then krange=[10,100]+0d0 
kmax = krange(1)+0d0
kmin = krange(0)+0d0
if not keyword_set(kswitch) then kswitch = kmax 
nk = 8192
dlogk = alog10(kmax/kmin)/(nk-1d0)
kaxis = krange(0) + (krange(1)-krange(0))*dindgen(nk)/(nk-1d0)

bvals = dblarr(nk)
freqs = dblarr(nk)

bvals1 = dblarr(nk)
freqs1 = dblarr(nk)

xvec = [0.45017755,-0.71502877];sol for kx=1, gmma=1.4, eps=0.05, mmode=0, smallq=-1 (bcrit for isothermal vsi)
xvec1= [1.5690316, -0.73509775];sol for kx=1, gmma=1.4, eps=0.05, mmode=0, smallq=-1 (bcrit for adiabatic vsi) 

;xvec = [0.12557000, -0.010043542];sol for kx=100, gmma=1.4, eps=0.05, mmode=0, smallq=-1 (bcrit for isothermal vsi)
;xvec1= [102.47491, -0.073093245 ];sol for kx=100, gmma=1.4, eps=0.05, mmode=0, smallq=-1 (bcrit for adiabatic vsi)

for j=0, nk-1 do begin
;for j=0, 0 do begin 

kx = kaxis(j)

Result = NEWTON(xvec, 'func', /DOUBLE, TOLF=1d-9, tolx=1d-9, itmax=50000)
bvals(j) = result(0)
freqs(j) = result(1)
xvec = result 


Result1 = NEWTON(xvec1, 'func', /DOUBLE, TOLF=1d-9, tolx=1d-9, itmax=50000)
bvals1(j) = result1(0)
freqs1(j) = result1(1)
xvec1 = result1
endfor

etavals = 1d0/bvals 

etavals1 = 1d0/bvals1



;xvec1 = [10,1]
;Result2 = NEWTON(xvec1, 'func0', /DOUBLE, TOLF=1d-9, tolx=1d-9, itmax=50000)
;b0=result2(0)
;om0=result2(1)
;print, 'b0,om0', b0, om0


  set_plot, 'ps'
  device, filename='bcrit_theory.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

  plot, kaxis, bvals,xmargin=[8.2,1.8],ymargin=[3.5,2] $
        ,ytitle=textoidl('\beta_{diff}'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title,/xlog ,/ylog
  oplot, kaxis, bvals1, thick=4

;  oplot, kaxis, b0*sqrt(kaxis), linestyle=1,thick=4 ;asymptotic limit for fund adia vsi 
  oplot, kaxis, abs(eps*smallq)/(gmma-1d0)*kaxis^2,  linestyle=1,thick=4; asymptotic limit for fund iso vsi

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
end
