pro compare_eigen, legend=legend, mmode=mmode, yrange=yrange, xrange=xrange, cases=cases 

  !p.font = 0

  params = dblarr(7,1)
  openr,1,'params.dat'
  readf, 1, params, format='(7(e22.15,x))'
  close,1
  gmma  = params(0,0)
  bgmma = params(1,0)
  eps   = params(2,0)
  smallh= params(3,0)
  smalls= params(4,0)

  nz    = file_lines('basic.dat')
  basic = dblarr(8,nz)
  openr,1,'basic.dat'
  readf,1,basic
  close,1
  zaxis = basic(0,*)*smallh/eps
  logrho= basic(1,*)
  omega2 = basic(4,*)
  kappa2= basic(6,*)
  csq   = basic(7,*)

  nk = file_lines('kaxis.dat')
  kaxis = dblarr(nk)
  openr,1,'kaxis.dat'
  readf,1,kaxis
  close,1
  
  nb = file_lines('gaxis.dat')
  baxis = dblarr(4,nb)
  openr,1,'gaxis.dat'
  readf,1,baxis
  close,1 

  nmodes = file_lines('eigenvalues.dat')
  eigenvalues = dblarr(3,nmodes)
  openr,1,'eigenvalues.dat'
  readf,1,eigenvalues
  close,1

  ncases = n_elements(cases)
  bvals  = dblarr(ncases)
  freq   = dblarr(ncases, nk)
  growth = dblarr(ncases, nk) 

  for i=0, ncases-1 do begin
  temp = min(abs(baxis(0,*) - cases(i)),grid)
  bvals(i) = baxis(0,grid)
  
  ibeg = grid*nk 
  freq(i,*)  = eigenvalues(0,ibeg:ibeg+nk-1)
  growth(i,*)= eigenvalues(1,ibeg:ibeg+nk-1) 

  endfor  

  label = '\gamma='+string(bvals,format='(f4.2)')
  color_arr = dindgen(ncases)
  for i=0, ncases-1 do color_arr(i) *= 256/ncases
 
  loadct,6
  set_plot, 'ps'
  device, filename='compare_eigen_imag.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color

  plot, kaxis, 10d0*abs(growth(0,*)) ,xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('10\times|\nu|/(\epsilon\Omega_k)'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5,xminor=10, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title, color=color_arr(0), /xlog
  for i=1, ncases-1 do begin
     oplot, kaxis, 10d0*abs(growth(i,*)), thick=4,color=color_arr(i)
  endfor
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

  set_plot, 'ps'
  device, filename='compare_eigen_real.ps'$
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color

  plot, kaxis, abs(freq(0,*)) ,xmargin=[7,3],ymargin=[3.5,2] $
        ,ytitle=textoidl('|\omega|/(\epsilon\Omega_k)'),xtitle=textoidl('k_xH_{iso}') $
        ,charsize=1.5,xminor=10, thick=4, xrange=xrange,xtickinterval=xtickinterval $
        ,yrange=yrange,ytickinterval=ytickinterval, title=title, color=color_arr(0),/xlog
  for i=1, ncases-1 do begin
     oplot, kaxis, abs(freq(i,*)), thick=4,color=color_arr(i)
  endfor
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
