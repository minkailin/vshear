pro find_fund

  
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

  nmodes = file_lines('eigenvalues.dat')
  nnodes = dblarr(nmodes)
  array=dblarr(14,long(nz)*nmodes)
  openr,1,'eigenvectors.dat'
  readf,1,array
  close,1

  for k=0, nmodes-1 do begin 

  eigenvz_real = array(10, (k-1)*nz : k*nz -1)
  eigenvz_imag = array(11, (k-1)*nz : k*nz -1)

  nnodes(k) = 0 
  for i=0, nz-2 do begin
  test1 = eigenvz_real(i)*eigenvz_real(i+1)
  test2 = eigenvz_imag(i)*eigenvz_imag(i+1) 
  if((test1 lt 0d0) and (test2 lt 0d0) then nodes+=1
  endfor
  
  temp = where(nnodes eq 0)
  print, 'mode numbers with no nodes in vz are'
  print, temp+1
    
end
