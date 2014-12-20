pro result, mmode=mmode, xrange=xrange, yrange=yrange

  !p.font = 0

  nz    = file_lines('basic.dat')
  basic = dblarr(7,nz)
  openr,1,'basic.dat'
  readf,1,basic
  close,1
  zaxis = basic(0,*)
  logrho= basic(1,*)
  omega = sqrt(basic(4,*))
  kappa2= basic(6,*)

  set_plot, 'ps'
  device, filename='omega.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, omega,xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min(omega),max(omega)] $
        , xtitle='z/H', ytitle=textoidl('\Omega/\Omega_k') $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close

  set_plot, 'ps'
  device, filename='kappa2.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, kappa2,xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min(kappa2),max(kappa2)] $
        , xtitle='z/H', ytitle=textoidl('\kappa^2/\Omega_k^2') $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close
  
  set_plot, 'ps'
  device, filename='omega.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, omega,xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min(omega),max(omega)] $
        , xtitle='z/H', ytitle=textoidl('\Omega/\Omega_k') $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close
  

  set_plot, 'ps'
  device, filename='logrho.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, logrho,xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min(logrho),max(logrho)] $
        , xtitle='z/H', ytitle=textoidl('ln(\rho/\rho_0)') $
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


  nmodes = file_lines('eigenvalues.dat')
  eigenvalues = dblarr(2,nmodes)
  openr,1,'eigenvalues.dat'
  readf,1,eigenvalues
  close,1
  freq  = eigenvalues(0,*)
  growth= eigenvalues(1,*) 
  
  set_plot, 'ps'
  device, filename='eigenvalues.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, abs(freq), abs(growth),xmargin=[6,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=textoidl('|\omega|/(\epsilon\Omega_k)'), ytitle=textoidl('|\nu|/(\epsilon\Omega_k)') $
        ,ytickinterval=ytickinterval, psym=7, symsize=1  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close


  temp = min(abs(freq),minloc)
  temp = max(abs(growth),maxloc)
  print, 'nmodes, min_real, max_imag', nmodes, minloc+1,maxloc+1
  
  if not keyword_set(mmode) then mmode = minloc+1 ;plot eigenfunction closest to pure growth
  print, 'mode, eigen=', minloc+1, eigenvalues(0,minloc), eigenvalues(1,minloc)
  print, 'mode, eigen=', maxloc+1, eigenvalues(0,maxloc), eigenvalues(1,maxloc)
  array=dblarr(3,long(nz)*nmodes)
  openr,1,'eigenvectors.dat'
  readf,1,array
  close,1
  eigen_real = array(1, (mmode-1)*nz : mmode*nz -1)
  eigen_imag = array(2, (mmode-1)*nz : mmode*nz -1)
  
  eigen_abs = sqrt(eigen_real^2 + eigen_imag^2)
  realfreq   = abs(eigenvalues(0,mmode-1))
  growthrate = abs(eigenvalues(1,mmode-1))
  title = '|\omega|='+string(realfreq,format='(f4.2)')+'\epsilon\Omega_k, '
  title+= '|\nu|='+string(growthrate,format='(f4.2)')+'\epsilon\Omega_k'
  title = textoidl(title)


   if not keyword_set(xrange) then xrange = [-10,10]
  if not keyword_set(yrange) then yrange = [min([min(eigen_real),min(eigen_imag)]),max([max(eigen_real),max(eigen_imag)])]/max(eigen_abs)


  set_plot, 'ps'
  device, filename='eigenvector.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zaxis, eigen_real/max(eigen_abs),xmargin=[8,2],ymargin=[3.5,2], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=textoidl('z/\epsilonr'), ytitle=textoidl('W/max|W|') $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1, title=title
  oplot, zaxis, eigen_imag/max(eigen_abs), thick=4, linestyle=1
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
