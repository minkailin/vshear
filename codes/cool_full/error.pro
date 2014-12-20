pro error, mmode=mmode, xrange=xrange, yrange=yrange

  ii = dcomplex(0d0, 1d0) 

  nz    = file_lines('basic.dat')
  basic = dblarr(8,nz)
  openr,1,'basic.dat'
  readf,1,basic, format='(8(e22.15,x))'
  close,1
  zaxis = basic(0,*)
  logrho= basic(1,*)
  dlogrho= basic(2,*) 
  d2logrho= basic(3,*) 
  omega = sqrt(basic(4,*))
  domega2= basic(5,*)  
  kappa2= basic(6,*)
  csq   = basic(7,*)

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

  nmodes = file_lines('eigenvalues.dat')
  eigenvalues = dblarr(2,nmodes)
  openr,1,'eigenvalues.dat'
  readf,1,eigenvalues,format='(2(e22.15,x))'
  close,1
  freq  = eigenvalues(0,*)
  growth= eigenvalues(1,*) 
  sigma = dcomplex(freq,growth)

  temp = min(abs(freq),minloc)
  temp = max(abs(growth),maxloc)
  print, 'nmodes, min_real, max_imag', nmodes, minloc+1,maxloc+1
  
  if not keyword_set(mmode) then mmode = minloc+1 ;plot eigenfunction closest to pure growth
  print, 'mode, eigen=', mmode, eigenvalues(0,mmode-1), eigenvalues(1,mmode-1)
  sig = sigma(mmode-1) 

  array=dblarr(14,long(nz)*nmodes)
  openr,1,'eigenvectors.dat'
  readf,1,array, format='(14(e22.15,x))'
  close,1
  bigW  = dcomplex( array(0, (mmode-1)*nz : mmode*nz -1), $
                    array(1, (mmode-1)*nz : mmode*nz -1)  )
  dbigW = dcomplex( array(2, (mmode-1)*nz : mmode*nz -1), $
                    array(3, (mmode-1)*nz : mmode*nz -1)  )
  bigQ = dcomplex( array(4, (mmode-1)*nz : mmode*nz -1), $
                    array(5, (mmode-1)*nz : mmode*nz -1)  )

  vx  = dcomplex( array(6, (mmode-1)*nz : mmode*nz -1), $
                    array(7, (mmode-1)*nz : mmode*nz -1)  )  
  vy  = dcomplex( array(8, (mmode-1)*nz : mmode*nz -1), $
                    array(9, (mmode-1)*nz : mmode*nz -1)  )
  vz  = dcomplex( array(10, (mmode-1)*nz : mmode*nz -1), $
                    array(11, (mmode-1)*nz : mmode*nz -1)  )
  dvz  = dcomplex( array(12, (mmode-1)*nz : mmode*nz -1), $
                    array(13, (mmode-1)*nz : mmode*nz -1)  )


;  lagp = -(dbigW + dlogrho*bigW)*csq*dlogrho + csq*dlogrho^2*bigQ - sig*sig*smallh^2*bigW

   lagp = bigW + ii*vz*csq*dlogrho/(sig*smallh)
   lagp = abs(lagp)

;  vz = dbigW + dlogrho*(bigW-bigQ) 
;  vz/= ii*sig*smallh 
;  dvz = d2bigW + d2logrho*(bigW-bigQ) + dlogrho*(dbigW - dbigQ) 
;  dvz/= ii*sig*smallh 

;  vx = (-1d0/bigD)*(sig*eps*kx*bigW/smallh + domega2*vz/smallh) 

  test1 = -ii*sig*smallh*bigQ/csq + ii*kx*vx/eps + dvz + dlogrho*vz 
  test2 = ii*eps*sig*vx + 2d0*omega*vy - ii*kx*bigW/smallh
  test3 = -kappa2*vx + 2d0*ii*eps*sig*omega*vy - domega2*vz/smallh
  test4 = -ii*smallh*sig*vz + dbigW + dlogrho*(bigW-bigQ)
  test5 = ii*sig*bigW - (bigW - bigQ/bgmma)/(eps*bcool) - ii*sig*gmma*bigQ/bgmma + csq*dlogrho*(gmma/bgmma-1d0)*vz/smallh

  print, 'max error', max([max(abs(test1)),max(abs(test2)),max(abs(test3)),max(abs(test4)),max(abs(test5))])
  print, 'bc error (nolagp)', max([lagp(0),lagp(nz-1)])
  print, 'bc error (nozvel)', max([abs(vz(0)),abs(vz(nz-1))])

  set_plot,'x'
  window,0
  plot, zaxis, abs(lagp), xrange=xrange,xstyle=1
;  oplot, zaxis, abs(test2), linestyle=1 
end
