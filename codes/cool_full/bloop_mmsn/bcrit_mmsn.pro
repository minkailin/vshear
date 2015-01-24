
pro bcrit_mmsn, legend=legend, label=label, xrange=xrange, yrange=yrange 

!p.font = 0


nrad = 1024

rmin=1
rmax=1d3

zeta = 10.^(alog10(rmin) + alog10(rmax/rmin)*dindgen(nrad)/(nrad-1d0))

gmma=1.4
mu = 2.0

kx = [1d0, 1d1, 1d2]

bcool = dblarr(3, nrad)

for k=0, 2 do begin

bcool(k,*)  = 4.4d5/mu/(gmma-1d0)*zeta^(-57./14.)
bcool(k,*) *= 8.3d-9*zeta^(33./7.) + 1d0/kx(k)^2
endfor

bcrit = 9.4d-3/(gmma-1d0)*zeta^(2./7.) 


  set_plot, 'ps'
  device, filename='bcrit_mmsn.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches

  plot, zeta, bcrit,xmargin=[7.8,2.2],ymargin=[3.5,2] $
        ,ytitle=textoidl('t_c\Omega_k'),xtitle=textoidl('\zeta=r/AU') $
        ,charsize=1.5, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title,/xlog ,/ylog
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
