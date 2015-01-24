pro reduce, nk=nk

  nb = file_lines('baxis.dat')
  baxis = dblarr(4,nb)
  openr,1,'baxis.dat'
  readf,1,baxis
  close,1

 nmodes = file_lines('eigenvalues.dat')
 eigenvalues = dblarr(3,nmodes)
 openr,1,'eigenvalues.dat'
 readf,1,eigenvalues, format='(3(e22.15,x))'
 close,1
 freq  = eigenvalues(0,*)
 growth= eigenvalues(1,*)


 nkold = nmodes/nb

 kaxis = eigenvalues(2,0:nkold-1)

 kmin = kaxis(0)
 kmax = kaxis(nkold-1)


 new_kaxis = dblarr(nk)
 dlogkx = alog10(kmax/kmin)/(nk-1d0)
 for i=0, nk-1 do new_kaxis(i) = 10d0^(alog10(kmin) + dlogkx*i)


 openw,1,'eigenvalues.dat'
 for j=0, nb-1 do begin

 new_freq    = spline(kaxis, freq(j*nkold:(j+1)*nkold-1), new_kaxis)
 new_growth  = spline(kaxis, growth(j*nkold:(j+1)*nkold-1), new_kaxis)

 
 for i=0, nk-1 do begin
 printf,1,new_freq(i),new_growth(i),new_kaxis(i), format='(3(e22.15,x))'
 endfor
 endfor
 close,1

end
