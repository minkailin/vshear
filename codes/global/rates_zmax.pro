pro rates_zmax, xrange=xrange, yrange=yrange

  !p.font = 0

  nz    = file_lines('rates_zmax.dat')
  basic = dblarr(3,nz)
  openr,1,'rates_zmax.dat'
  readf,1,basic
  close,1
  zmax = basic(0,*)
  freq = basic(1,*)
  growth = abs(basic(2,*))

  max_theory = 0.5*zmax 
 
  set_plot, 'ps'
  device, filename='rates_zmax.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, zmax, growth,xmargin=[8,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=[min(growth),max(growth)] $
        , xtitle=textoidl('z_{max}/H_{iso}'), ytitle=textoidl('max|\nu|/(\epsilon\Omega_k)') $
        ,ytickinterval=ytickinterval  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  oplot, zmax, max_theory, linestyle=1, thick=4
  device,/close

end
