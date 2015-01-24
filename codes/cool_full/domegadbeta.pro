pro domegadbeta, eps=eps, smallq=smallq, gmma=gmma, krange=krange, yrange=yrange, xrange=xrange, legend=legend, label=label

!p.font = 0

ii = dcomplex(0d0, 1d0)

if not keyword_set(krange) then krange=[0.1,100]
nk = 4096
kaxis = krange(0) + (krange(1)-krange(0))*dindgen(nk)/(nk-1d0)

domega_dbeta = dblarr(nk)
dnu_dbeta    = dblarr(nk)

for j=0, nk-1 do begin
khat = kaxis(j)
ksq  = khat^2

temp = -ii*(gmma-1d0)*ksq*(ii*eps*smallq - khat)^2
temp/= 2d0*(1d0+ii*eps*smallq*khat)*(1d0+ksq)^2

domega_dbeta(j) = real_part(temp)
dnu_dbeta(j)    = imaginary(temp)

endfor

  set_plot, 'ps'
  device, filename='domegadbeta.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

  plot, kaxis, 10.*domega_dbeta,xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('10\times\Omega_k^{-2}\partial\sigma/\partialt_c at t_c=0'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title,/xlog 
  oplot, kaxis, 10.*dnu_dbeta, thick=4, linestyle=1
  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, numcases-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
  endif
  device,/close

end
