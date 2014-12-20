pro theory, mmax=mmax, xrange=xrange, yrange=yrange, eps=eps, smallq=smallq, kx=kx
  !p.font = 0
  
  ii = dcomplex(0d0, 1d0)
  kx/= eps

  eigenvalues = dcomplexarr(mmax)
  for i=0, mmax-1 do begin
     sig2 = (i+1d0)*(1d0 + ii*smallq*eps*kx)/(1d0 + kx^2)
     eigenvalues(i) = sqrt(sig2)
  endfor

  freq  = real_part(eigenvalues(*))/eps
  growth= imaginary(eigenvalues(*))/eps
  
  set_plot, 'ps'
  device, filename='eigenvalues_theory.ps' $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches
  plot, abs(freq), abs(growth),xmargin=[6,2],ymargin=[3.5,0.5], ystyle=0  $
        ,charsize=1.5, thick=4, xrange=xrange, yrange=yrange $
        , xtitle=textoidl('|\omega|/(\epsilon\Omega_k)'), ytitle=textoidl('|\nu|/(\epsilon\Omega_k)') $
        ,ytickinterval=ytickinterval, psym=7, symsize=1  $ 
        ,xtickinterval=xtickinterval, xstyle=1
  device,/close

end
