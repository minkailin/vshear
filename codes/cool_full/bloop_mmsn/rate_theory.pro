pro rate_theory, rad=rad,  krange=krange,  legend=legend, mmode=mmode, yrange=yrange, xrange=xrange, root=root

!p.font = 0

;MMSN

mu = 2.0
smallq = -0.42857143
gmma   = 1.4
eps    = 0.022*rad^(2./7.)

ii = dcomplex(0d0, 1d0)

if not keyword_set(mmode) then mmode=0d0
title = 'M='+string(mmode,format='(I2)')

if not keyword_set(krange) then krange=[0.1d0,100d0]
kmax = krange(1)
kmin = krange(0)
nk = 8192*10
kaxis = krange(0) + (krange(1)-krange(0))*dindgen(nk)/(nk-1d0)

rates = dblarr(nk)
freqs = dblarr(nk)

for j=0, nk-1 do begin
khat = kaxis(j)
ksq  = khat^2

bcool = (4.4d5/mu/(gmma-1d0))*rad^(-57./14.)
bcool*= 8.3d-9*rad^(33./7.) + 1d0/ksq 

bigA = 1d0 + ii*eps*smallq*khat 

c0 = mmode*(mmode+1d0)*bigA^2

c1 = ii*bcool*( (1d0-gmma) - ksq*(1d0+2d0*mmode)^2*(gmma-1d0) + 4d0*ii*eps*smallq*khat*mmode*(mmode+1d0)*(gmma-1d0) - 2d0*bigA^2*gmma*mmode*(mmode+1d0) ) 

c2 = (ksq+1d0)*bigA + bcool^2*( (1d0-gmma) + gmma*(1d0-gmma)*( ksq*(1d0+2d0*mmode)^2 - 4d0*ii*eps*smallq*khat*mmode*(mmode+1d0) ) - gmma^2*bigA^2*mmode*(mmode+1d0) )

c3 = bcool*( -3d0*ii - 2d0*ii*ksq + eps*khat*smallq + gmma*(ii + smallq*eps*khat*(1d0+2d0*ksq)) )

c4 = -(1d0+ksq)^2 + bcool^2*(1d0+gmma*ksq)*(-2d0 + gmma - ii*eps*smallq*khat*gmma) 

c5 = 2d0*ii*bcool*(1d0+ksq)*(1d0+gmma*ksq) 

c6 = bcool^2*(1d0+gmma*ksq)^2 

coeff = [c0, c1, c2, c3, c4, c5, c6]

result = fz_roots(coeff,eps=1d-9,/double)

if(j eq 0) then begin
sig = result(root(0)-1)
endif else begin
sigold = dcomplex(freqs(j-1),rates(j-1))
temp = min(abs(result-sigold),grid)
sig = result(grid)
endelse

rates(j) = imaginary(sig)
freqs(j) = real_part(sig)

endfor

rates /= eps 
freqs /= eps

loadct,6
  set_plot, 'ps'
  device, filename='rate_theory_grow.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color

  plot, kaxis, 10d0*rates(*),xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('10\times\nu/(\epsilon\Omega_k)'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title,/xlog 

  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, ncases-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
  endif
  device,/close

  set_plot, 'ps'
  device, filename='rate_theory_freq.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color

  plot, kaxis, freqs(*),xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('\omega/(\epsilon\Omega_k)'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title, /xlog

  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, ncases-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
  endif
  device,/close

end
