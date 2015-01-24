pro rate_theory, bmin=bmin, eps=eps, smallq=smallq, gmma=gmma, krange=krange, bmax=bmax, ncases=ncases, legend=legend, mmode=mmode, yrange=yrange, xrange=xrange, root=root, kswitch=kswitch, nkout=nkout

!p.font = 0

ii = dcomplex(0d0, 1d0)

if not keyword_set(mmode) then mmode=0d0
title = 'M='+string(mmode,format='(I2)')

if not keyword_set(bmin) then bmin = 0.001
if not keyword_set(ncases) then ncases = 6
dlogb = alog10(bmax/bmin)/(ncases-1d0)
baxis = 10d0^(alog10(bmin) + dlogb*dindgen(ncases))

label = '\beta='+string(baxis,format='(f6.1)')
color_arr = dindgen(ncases)
for i=0, ncases-1 do color_arr(i) *= 256/ncases

if not keyword_set(krange) then krange=[0.1d0,100d0]
kmax = krange(1)
kmin = krange(0)
if not keyword_set(kswitch) then kswitch = kmax 
nk = 8192*10
kaxis = krange(0) + (krange(1)-krange(0))*dindgen(nk)/(nk-1d0)

rates = dblarr(ncases, nk)
freqs = dblarr(ncases, nk)
alphap= dblarr(ncases, nk)
alpham= dblarr(ncases, nk)

for i=0, ncases-1 do begin
bcool = baxis(i)
for j=0, nk-1 do begin
khat = kaxis(j)
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

result = fz_roots(coeff,eps=1d-9,/double)

if(j eq 0) then begin
sig = result(root(0)-1)
endif else begin
sigold = dcomplex(freqs(i,j-1),rates(i,j-1))
temp = min(abs(result-sigold),grid)
sig = result(grid)
endelse

rates(i,j) = imaginary(sig)
freqs(i,j) = real_part(sig)

chi = (1d0 - ii*sig*bcool)/(1d0 - ii*sig*bcool*gmma)
bigC = (1d0-chi)*(ksq - ii*eps*smallq*khat)

alphap(i,j) =0.5d0*( bigA + sqrt(bigA^2 + 4d0*bigC) )
alpham(i,j) =0.5d0*( bigA - sqrt(bigA^2 + 4d0*bigC) )

endfor
endfor

rates /= eps 
freqs /= eps




if not keyword_set(nkout) then nkout = nk
dlogkx = alog10((kmax+0d0)/(kmin+0d0))/(nkout-1d0)
new_kaxis = 10d0^(alog10(kmin+0d0) + dlogkx*dindgen(nkout))


openw,1,'rate_theory.dat'
for i=0, ncases-1 do begin
   new_freqs    = spline(kaxis, freqs(i,*), new_kaxis)
   new_rates    = spline(kaxis, rates(i,*), new_kaxis)
   for j=0, nkout-1 do begin 
      printf, 1, new_freqs(j), new_rates(j), new_kaxis(j),format='(3(e22.15,x))'
   endfor
endfor
close,1


loadct,6
  set_plot, 'ps'
  device, filename='rate_theory_grow.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color

  plot, kaxis, 10d0*rates(0,*),xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('10\times\nu/(\epsilon\Omega_k)'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title, color=color_arr(0),/xlog 

  for i=1, ncases-1 do begin
     oplot, kaxis, 10d0*rates(i,*), thick=4,color=color_arr(i)
  endfor
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

  plot, kaxis, freqs(0,*),xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('\omega/(\epsilon\Omega_k)'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title, color=color_arr(0),/xlog

  for i=1, ncases-1 do begin
     oplot, kaxis, freqs(i,*), thick=4,color=color_arr(i)
  endfor
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


  sig2_abs = (freqs*eps)^2 + (rates*eps)^2
  set_plot, 'ps'
  device, filename='rate_theory_sig.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color

  plot, kaxis, sig2_abs(0,*),xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('|\sigma|^2/\Omega_k^2'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title, color=color_arr(0),/xlog

  for i=1, ncases-1 do begin
     oplot, kaxis, sig2_abs(i,*), thick=4,color=color_arr(i)
  endfor
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
  device, filename='rate_theory_alphap.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color

  plot, kaxis, alphap(0,*),xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('Re(\alpha)'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=[min(alphap),max(alphap)],ytickinterval=ytickinterval, title=title, color=color_arr(0),/xlog

  for i=1, ncases-1 do begin
     oplot, kaxis, alphap(i,*), thick=4,color=color_arr(i)
  endfor
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
  device, filename='rate_theory_alpham.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color

  plot, kaxis, alpham(0,*),xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('Re(\alpha)'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=[min(alpham),max(alpham)],ytickinterval=ytickinterval, title=title, color=color_arr(0),/xlog

  for i=1, ncases-1 do begin
     oplot, kaxis, alpham(i,*), thick=4,color=color_arr(i)
  endfor
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
