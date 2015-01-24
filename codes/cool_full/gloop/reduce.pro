pro reduce, nk=nk

 nmodes = file_lines('eigenvalues.dat')
 eigenvalues = dblarr(3,nmodes)
 openr,1,'eigenvalues.dat'
 readf,1,eigenvalues, format='(3(e22.15,x))'
 close,1
 freq  = eigenvalues(0,*)
 growth= eigenvalues(1,*)
 kaxis = eigenvalues(2,*)

 kmin = kaxis(0)
 kmax = kaxis(nmodes-1)

; new_kaxis = kmin + (kmax - kmin)*dindgen(nk)/(nk-1d0)

 new_kaxis = dblarr(nk)
 dlogkx = alog10(kmax/kmin)/(nk-1d0)
 for i=0, nk-1 do new_kaxis(i) = 10d0^(alog10(kmin) + dlogkx*i)

 new_freq    = spline(kaxis, freq, new_kaxis)
 new_growth  = spline(kaxis, growth, new_kaxis)

 
 openw,1,'eigenvalues.dat'
 for i=0, nk-1 do begin
 printf,1,new_freq(i),new_growth(i),new_kaxis(i), format='(3(e22.15,x))'
 endfor
 close,1

end
