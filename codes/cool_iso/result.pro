pro result, mmode=mmode, xrange=xrange 

nz    = file_lines('basic.dat')
basic = dblarr(7,nz)
openr,1,'basic.dat'
readf,1,basic
close,1
zaxis = basic(0,*)

nmodes = file_lines('eigenvalues.dat')
eigenvalues = dblarr(2,nmodes)
openr,1,'eigenvalues.dat'
readf,1,eigenvalues
close,1
freq  = eigenvalues(0,*)
growth= eigenvalues(1,*) 

temp = min(abs(freq),minloc)
temp = max(abs(growth),maxloc)
print, 'nmodes, min_real, max_imag', nmodes, minloc+1,maxloc+1

if not keyword_set(mmode) then mmode = minloc+1 ;plot eigenfunction closest to pure growth
print, 'mode, eigen=', mmode, eigenvalues(0,mmode-1), eigenvalues(1,mmode-1)
array=dblarr(4,long(nz)*nmodes)
openr,1,'eigenvectors.dat'
readf,1,array
close,1
eigen_real = array(0, (mmode-1)*nz : mmode*nz -1)
eigen_imag = array(1, (mmode-1)*nz : mmode*nz -1)

if not keyword_set(xrange) then xrange = [min(zaxis),max(zaxis)]

set_plot,'x'
window,0
plot, zaxis, eigen_real, thick=4, linestyle=0, xstyle=1, xrange=xrange
oplot,zaxis, eigen_imag,thick=4,linestyle=1

window,1
plot, abs(eigenvalues(0,*)), abs(eigenvalues(1,*)), psym=2, symsize=1.5

;basic = dblarr(7,nz)
;openr,1,'basic.dat'
;readf,1,basic
;close,1
;stop
;window,1
;plot, basic(0,*), basic(2,*), thick=4, xrange=[-1,-0.95]

end
