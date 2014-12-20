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
  bigD= kappa2 ;- sig^2*eps^2 

  array=dblarr(10,long(nz)*nmodes)
  openr,1,'eigenvectors.dat'
  readf,1,array, format='(10(e22.15,x))'
  close,1
  bigW  = dcomplex( array(0, (mmode-1)*nz : mmode*nz -1), $
                    array(1, (mmode-1)*nz : mmode*nz -1)  )
  dbigW = dcomplex( array(2, (mmode-1)*nz : mmode*nz -1), $
                    array(3, (mmode-1)*nz : mmode*nz -1)  )
  bigQ = dcomplex( array(4, (mmode-1)*nz : mmode*nz -1), $
                    array(5, (mmode-1)*nz : mmode*nz -1)  )

  vz   = dcomplex( array(6, (mmode-1)*nz : mmode*nz -1), $
                   array(7, (mmode-1)*nz : mmode*nz -1)  )
  dvz   = dcomplex( array(8, (mmode-1)*nz : mmode*nz -1), $
                    array(9, (mmode-1)*nz : mmode*nz -1)  )


lagp = -(dbigW + dlogrho*bigW)*csq*dlogrho + csq*dlogrho^2*bigQ - sig*sig*smallh^2*bigW

;  vz = dbigW + dlogrho*(bigW-bigQ) 
;  vz/= ii*sig*smallh 
;  dvz = dcomplex(deriv(zaxis,real_part(vz)),deriv(zaxis,imaginary(vz)))

;  dvz = d2bigW + d2logrho*(bigW-bigQ) + dlogrho*(dbigW - dbigQ) 
;  dvz/= ii*sig*smallh 

  vx = (-1d0/bigD)*(sig*eps*kx*bigW/smallh + domega2*vz/smallh) 

  test = -ii*sig*smallh*bigQ/csq + ii*kx*vx/eps + dvz + dlogrho*vz 

  set_plot,'x'
  window,0
  plot, zaxis, abs(test), xrange=xrange,xstyle=1 
stop 
end
