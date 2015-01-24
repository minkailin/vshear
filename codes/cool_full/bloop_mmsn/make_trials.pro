pro make_trials, krep=krep, brep=brep   
  
; table of cooling times
  nb = file_lines('baxis.dat')
  baxis = dblarr(4,nb)
  openr,1,'baxis.dat'
  readf,1,baxis
  close,1

; eigenvalues calculated by vsi.exe
  nmodes = file_lines('eigenvalues.dat')
  eigenvalues = dblarr(3,nmodes)
  openr,1,'eigenvalues.dat'
  readf,1,eigenvalues
  close,1

; theoretical eigenvalues
  nmodes = file_lines('rate_theory.dat')
  eigenvalues_theory = dblarr(3,nmodes)
  openr,1,'rate_theory.dat'
  readf,1,eigenvalues_theory
  close,1

  nk = nmodes / nb 

; replace selected eigenvlaues
; brep is the cooling time, krep is the radial wavenumber beyond which
; eigenvalues are replaced by theory

  temp = min(abs(baxis(0,*) - brep),bgrid)
  
  kbeg = bgrid*nk
  kend = kbeg + nk - 1


  kaxis = eigenvalues(2,*)
  growth= eigenvalues(1,kbeg:kend)
  growth_theory = eigenvalues_theory(1,kbeg:kend)


  for j=kbeg, kend do begin
     k = eigenvalues(2,j)
     if(k ge krep) then begin
        eigenvalues(0,j) = eigenvalues_theory(0,j)
        eigenvalues(1,j) = eigenvalues_theory(1,j)
     endif
  endfor
  
  openw,1,'eigenvalues.dat'
  for i=0, nmodes-1 do begin
     printf,1,eigenvalues(0:2,i), format='(3(e22.15,x))'
  endfor
  close,1
  
  
end
