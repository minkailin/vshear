pro result, mmode=mmode, xrange=xrange, yrange=yrange, scale=scale

  
  !p.font = 0
  
  ii = dcomplex(0d0, 1d0)
  if not keyword_set(scale) then scale=1d0

  params = dblarr(7,1)
  openr,1,'params.dat'
  readf, 1, params, format='(7(e22.15,x))'
  close,1
  gmma  = params(0,0)
  bgmma = params(1,0)
  eps   = params(2,0)
  smallh= params(3,0)
  smalls= params(4,0)
  bcool = params(5,0)
  kx    = params(6,0)
 

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

  set_plot, 'ps'
  device, filename='omega2.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, omega2,xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min(omega2),max(omega2)] $
        , xtitle=('z/\epsilonr'), ytitle=textoidl('\Omega^2/\Omega_k^2') $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close

  set_plot, 'ps'
  device, filename='kappa2.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, kappa2,xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min(kappa2),max(kappa2)] $
        , xtitle='z/\epsilonr', ytitle=textoidl('\kappa^2/\Omega_k^2') $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close
  
 
  logrho = exp(logrho) 
  set_plot, 'ps'
  device, filename='logrho.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, logrho,xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min(logrho),max(logrho)] $
        , xtitle='z/\epsilonr', ytitle=textoidl('\rho/\rho_0') $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  ;; if keyword_set(legend) then begin
  ;;    x0=legend(0)
  ;;    x1=legend(1)
  ;;    y0=legend(2)
  ;;    dy=legend(3)
  ;;    for j=0, n_elements(label)-1 do begin
  ;;       oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
  ;;       xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
  ;;    endfor
  ;; endif
  device,/close

  set_plot, 'ps'
  device, filename='csq.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, csq,xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min(csq),max(csq)] $
        , xtitle='z/\epsilonr', ytitle=textoidl('c_s^2/c_s^2(0)') $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close


  nmodes = file_lines('eigenvalues.dat')
  eigenvalues = dblarr(3,nmodes)
  openr,1,'eigenvalues.dat'
  readf,1,eigenvalues
  close,1
  freq  = eigenvalues(0,*)
  growth= eigenvalues(1,*) 
  
  set_plot, 'ps'
  device, filename='eigenvalues.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, abs(freq), abs(growth)*10d0,xmargin=[6,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=textoidl('|\omega|/(\epsilon\Omega_k)'), ytitle=textoidl('10\times|\nu|/(\epsilon\Omega_k)') $
        ,ytickinterval=ytickinterval, psym=7, symsize=1  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close


  temp = min(abs(freq),minloc)
  temp = max(abs(growth),maxloc)
  print, 'nmodes, min_real, max_imag', nmodes, minloc+1,maxloc+1
  
  if not keyword_set(mmode) then mmode = minloc+1 ;plot eigenfunction closest to pure growth
  print, 'mode, eigen=', minloc+1, eigenvalues(0,minloc), eigenvalues(1,minloc)
  print, 'mode, eigen=', maxloc+1, eigenvalues(0,maxloc), eigenvalues(1,maxloc)

  array=dblarr(14,long(nz)*nmodes)
  openr,1,'eigenvectors.dat'
  readf,1,array
  close,1
  eigenW_real = array(0, (mmode-1)*nz : mmode*nz -1)
  eigenW_imag = array(1, (mmode-1)*nz : mmode*nz -1)
  eigenW = scale*dcomplex(eigenW_real,eigenW_imag) 
  eigenW_real = real_part(eigenW)
  eigenW_imag = imaginary(eigenW)

  eigenQ_real = array(4, (mmode-1)*nz : mmode*nz -1)
  eigenQ_imag = array(5, (mmode-1)*nz : mmode*nz -1)  
  eigenQ = scale*dcomplex(eigenQ_real,eigenQ_imag)
  eigenQ_real = real_part(eigenQ)
  eigenQ_imag = imaginary(eigenQ)

  if not keyword_set(xrange) then xrange = [min(zaxis),max(zaxis)]
  
  eigen_abs = abs(eigenW) 
  
  eigenW_real /= max(eigen_abs)
  eigenW_imag /= max(eigen_abs)

  eigenQ_real /= max(eigen_abs)
  eigenQ_imag /= max(eigen_abs)

  realfreq   = (eigenvalues(0,mmode-1))
  growthrate = (eigenvalues(1,mmode-1))
  title = '\omega='+string(realfreq,format='(f6.2)')+'\epsilon\Omega_k, '
  title+= '\nu='+string(growthrate,format='(f4.2)')+'\epsilon\Omega_k, '
  title+= 't_c\Omega_k='+string(bcool,format='(f5.3)')
  title = textoidl(title)
  set_plot, 'ps'
  device, filename='eigenvectorW.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, eigenW_real,xmargin=[8,2],ymargin=[3.5,2], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min([min(eigenW_real),min(eigenW_imag)]),max([max(eigenW_real),max(eigenW_imag)])] $
        , xtitle=textoidl('z/\epsilonr'), ytitle=textoidl('W/max|W|') $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1, title=title
  oplot, zaxis, eigenW_imag, thick=4, linestyle=1
  device,/close
 
  set_plot, 'ps'
  device, filename='eigenvectorQ.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, eigenQ_real, xmargin=[8,2],ymargin=[3.5,2], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min([min(eigenW_real),min(eigenW_imag)]),max([max(eigenW_real),max(eigenW_imag)])] $
        , xtitle=textoidl('z/\epsilonr'), ytitle=textoidl('Q/max|W|') $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1, title=title
  oplot, zaxis, eigenQ_imag, thick=4, linestyle=1
  device,/close

  eigenvz_real = array(10, (mmode-1)*nz : mmode*nz -1)
  eigenvz_imag = array(11, (mmode-1)*nz : mmode*nz -1)
  eigenvz = scale*dcomplex(eigenvz_real, eigenvz_imag)
  eigenvz_real = real_part(eigenvz)
  eigenvz_imag = imaginary(eigenvz) 
  eigen_abs = abs(eigenvz)

  eigenvz_real /= max(eigen_abs)
  eigenvz_imag /= max(eigen_abs)

  set_plot, 'ps'
  device, filename='eigenvectorvz.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, eigenvz_real, xmargin=[8,2],ymargin=[3.5,2], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min([min(eigenvz_real),min(eigenvz_imag)]),max([max(eigenvz_real),max(eigenvz_imag)])] $
        , xtitle=textoidl('z/\epsilonr'), ytitle=textoidl('\deltav_z/max|\deltav_z|') $
        ,ytickinterval=ytickinterval  $
        ,xtickinterval=xtickinterval, xstyle=1, title=title
  oplot, zaxis, eigenvz_imag, thick=4, linestyle=1
  device,/close


  ;; set_plot,'x'
  ;; window,0
  ;; plot, zaxis, eigen_real, thick=4, linestyle=0, xstyle=1, xrange=xrange
  ;; oplot,zaxis, eigen_imag,thick=4,linestyle=1
  
  ;; window,1
  ;; plot, abs(eigenvalues(0,*)), abs(eigenvalues(1,*)), psym=1, symsize=1
  ;; xyouts, abs(eigenvalues(0,mmode-1)), abs(eigenvalues(1,mmode-1)), 'mode',charthick=2, charsize=1.5
;basic = dblarr(7,nz)
;openr,1,'basic.dat'
;readf,1,basic
;close,1
;stop
;window,1
;plot, basic(0,*), basic(2,*), thick=4, xrange=[-1,-0.95]
end
