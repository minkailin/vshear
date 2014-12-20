module global
  character*6 :: vbc
  integer :: nz, lwork, info, nzmid, nzeff
  integer, allocatable :: ipiv(:), eigen_test(:)
  real*8, parameter :: pi = 2d0*acos(0d0)
  real*8 :: smallq, smallp, smallh, kx, zmax, zbar,lmode, T_l, dT_l, d2T_l, T_lp2, dT_lp2, d2T_lp2
  real*8 :: dz
  real*8, allocatable :: zaxis(:), rwork(:), freq(:), growth(:)
  real*8, allocatable :: logrho(:), dlogrho(:), omega2(:), domega2(:), kappa2(:)
  complex*16, parameter :: ii = (0d0,1d0)
  complex*16 :: sig2, sig
  complex*16, allocatable :: bigpi(:),  T(:,:), Tp(:,:), Tpp(:,:), matrix(:,:), w(:), vl(:,:), vr(:,:), rhs(:,:)
  complex*16, allocatable :: work(:), alpha(:), beta(:)
end module global

program vsi
  use global
  implicit none
  character*1, parameter :: JOBVL = 'N', JOBVR = 'V'
  integer :: i,j,k, loc(1), lmax, ip1
  real*8, external :: dlogrho_dz, domega2_dz, omega2_z, shear_z, logrho_z
  real*8 :: m
  namelist /params/ smallq, smallp, smallh, kx, zmax, vbc
  namelist /grid/ nz
  
  !read input parameters.
  open(7, file="input")
  read(7, nml=params)
  read(7, nml=grid)
  close(7)
  
  if(mod(nz,2).eq.0) then
     print*, 'Nz needs to be odd but Nz=', nz
     stop
  else
     Nzmid = (Nz+1)/2
  endif
  
  !allocate grids
  allocate(zaxis(nz))
  allocate(logrho(nz))
  allocate(dlogrho(nz))
  allocate(omega2(nz))
  allocate(domega2(nz))
  allocate(kappa2(nz))
  allocate(bigpi(nz))
  
  allocate(T(nz,nz))
  allocate(Tp(nz,nz))
  allocate(Tpp(nz,nz))
  
  
  if(vbc.ne.'nolagp') then
     nzeff = nz-2
  else
     nzeff = nz
  endif

  allocate(matrix(nzeff,nzeff))
  allocate(rhs(nzeff,nzeff))
  allocate(vr(nzeff,nzeff))
  allocate(vl(nzeff,nzeff))
  allocate(w(nzeff))
  allocate(growth(nzeff))
  allocate(freq(nzeff))
  allocate(eigen_test(nzeff))   

  lwork = 8*nzeff
  allocate(work(lwork))
  allocate(rwork(2*nzeff))
  allocate(ipiv(nzeff))
  
  !setup Z axis. these are extrema of T_lmax(Z/zmax) plus end points
  lmax = nz-1
  zaxis(Nzmid) = 0d0
  do j=Nzmid+1, nz
     zaxis(j) = -zmax*cos(pi*(j-1d0)/lmax)
  enddo
  do j=1, Nzmid-1
     zaxis(j) = -zaxis(Nz - j + 1)
  enddo
  
  !basic state
   !setup basic state 
  do i=1, nz
     
     dlogrho(i)  = dlogrho_dz(zaxis(i))
     domega2(i)  = domega2_dz(zaxis(i))
     omega2(i)   = omega2_z(zaxis(i))
     if(omega2(i).le.0d0) then
     print*, 'omega2<=0'
     stop
     endif
     kappa2(i)   = 2d0*sqrt(omega2(i))*(2d0*sqrt(omega2(i)) + shear_z(zaxis(i)))
     logrho(i)  = logrho_z(zaxis(i))
  enddo

  !output basic state
  open(10,file='basic.dat')
  do i=1, nz
     write(10,fmt='(7(e22.15,x))'), zaxis(i), logrho(i), dlogrho(i), 0.0, omega2(i), domega2(i), kappa2(i)
  enddo
  close(10)

  !set up basis matrices

  if (vbc.eq.'nopert') then
     do j=1, nz !jth physical grid (endpoints not actually used until reconstruction)
        zbar = zaxis(j)/zmax
        do k=1, nzeff !kth basis (only T_2 ... T_lmax used)
           lmode = dble(k) + 1d0
           call chebyshev_poly(lmode, zbar, T_l, dT_l, d2T_l)
           
           if(mod(k+1,2).eq.0) then
              
              T(j,k)  = T_l - 1d0
              Tp(j,k) = dT_l
              Tpp(j,k)= d2T_l
              
           else
              T(j,k)  = T_l - zbar
              Tp(j,k) = dT_l - 1d0 
              Tpp(j,k)= d2T_l
           endif
           Tp(j,k)  = Tp(j,k)/zmax          !convert to deriv wrt physical grid
           Tpp(j,k) = Tpp(j,k)/zmax**2d0    !convert to deriv wrt physical grid
           
        enddo
     enddo
  endif
  
  if (vbc.eq.'nograd') then
     do j=1, nz !jth physical grid  (endpoints not actually used until reconstruction)
        zbar = zaxis(j)/zmax

        k=1
        T(j,k)  =   1d0
        Tp(j,k) =   0d0
        Tpp(j,k)=   0d0

        do k=2, nzeff !kth basis (only T_0 ... T_lmax-2 used)
           m = dble(k-1d0)
           call chebyshev_poly(m,     zbar, T_l, dT_l, d2T_l)
           call chebyshev_poly(m+2d0, zbar, T_lp2, dT_lp2, d2T_lp2)

           if(mod(k-1,2).eq.0) then
              
              
              T(j,k)  =   T_l - ((m/2d0)/(m/2d0+1d0))**2d0*T_lp2
              Tp(j,k) =  dT_l - ((m/2d0)/(m/2d0+1d0))**2d0*dT_lp2
              Tpp(j,k)= d2T_l - ((m/2d0)/(m/2d0+1d0))**2d0*d2T_lp2
              

           else
              T(j,k)  =   T_l - (m/(m+2d0))**2d0*T_lp2
              Tp(j,k) =  dT_l - (m/(m+2d0))**2d0*dT_lp2
              Tpp(j,k)= d2T_l - (m/(m+2d0))**2d0*d2T_lp2
              
            
           endif
           Tp(j,k)  = Tp(j,k)/zmax          !convert to deriv wrt physical grid
           Tpp(j,k) = Tpp(j,k)/zmax**2d0    !convert to deriv wrt physical grid
           
        enddo
     enddo
  endif
  
  if ((vbc.eq.'nolagp').or.(vbc.eq.'noxvel')) then
     do j=1, nz !jth physical grid 
        zbar = zaxis(j)/zmax
        do k=1, nzeff !kth basis
           m = dble(k-1d0)
           call chebyshev_poly(m,     zbar, T_l, dT_l, d2T_l)
           T(j,k)  =   T_l
           Tp(j,k) =  dT_l
           Tpp(j,k)= d2T_l
           Tp(j,k)  = Tp(j,k)/zmax          !convert to deriv wrt physical grid
           Tpp(j,k) = Tpp(j,k)/zmax**2d0    !convert to deriv wrt physical grid
        enddo
     enddo
  endif
  
  !set up matrices
  do i=1, nzeff
     if((vbc.eq.'nolagp').or.(vbc.eq.'noxvel')) then
        ip1 = i+1 !for nograd or nopert vbc (boundary points not used in matrix)
     else
        ip1 = i   !for nolagp vbc (impost boundary condition at end points)
     endif
     matrix(i,:) =  Tpp(ip1,:) + (dlogrho(ip1) - ii*kx*domega2(ip1)/kappa2(ip1)/smallh**2d0)*Tp(ip1,:) 
     matrix(i,:) =  -matrix(i,:)/smallh**2d0
     matrix(i,:) =   matrix(i,:)/(1d0 + kx*kx/kappa2(ip1)/smallh**2d0)
     rhs(i,:)    =   T(ip1,:)
  enddo
 
  !explicit boundary condition for 'nolagp'
  if (vbc.eq.'nolagp') then
     do i=1,nzeff,nzeff-1
        matrix(i,:) = -Tp(i,:)*dlogrho(i)/smallh**2d0
        rhs(i,:)    =   T(i,:)
     enddo
  endif
  
  if (vbc.eq.'noxvel') then
     do i=1,nzeff,nzeff-1
        matrix(i,:) = -Tp(i,:)*domega2(i)/(ii*kx*smallh**2d0)
        rhs(i,:)    =   T(i,:)
     enddo
  endif



  call zgetrf (nzeff, nzeff, rhs, nzeff, IPIV, INFO)
  call ZGETRI(nzeff, rhs, nzeff, IPIV, WORK, LWORK, INFO )
  
  matrix = matmul(rhs,matrix)
  
  call zgeev (JOBVL, JOBVR, nzeff, matrix, nzeff, W, vl, nzeff, VR, nzeff, WORK, LWORK, RWORK, INFO)

  !w is sig^2 
  w = sqrt(w)
  !filter out unwanted eigenvalues (if too large)
  do i=1, nzeff 
!  if(dimag(w(i)) .le. 0d0) w(i) = 0d0
  if((abs(w(i)) .gt. 1d0/smallh).or.(abs(dimag(w(i))).le.smallh)) then 
  w(i) = 0d0 
  eigen_test(i) = 0
  else 
  eigen_test(i) = 1
  endif 
  enddo

  if(maxval(eigen_test).le.0) then
  print*, 'no good eigenvalues found'
  stop
  else
  print*, 'total number of  modes=', sum(eigen_test)
  endif


  freq = dble(w)
  growth=dimag(w)
!  loc = maxloc(abs(growth))
!  write(6,fmt='((e22.15,x),I04)') growth(loc(1)), loc(1)
  
  open(10,file='eigenvectors.dat')
  do k=1, nzeff
     if(eigen_test(k) .gt. 0) then  
     bigpi = matmul(T,vr(:,k))
     do i=1, nz
        write(10,fmt='(3(e22.15,x))'), zaxis(i), dble(bigpi(i)), dimag(bigpi(i))
     enddo
     endif 
  enddo
  close(10)
  
  open(20,file='eigenvalues.dat')
  do i=1, nzeff
     if(eigen_test(i) .gt. 0) then
     write(20,fmt='(2(e22.15,x))'), dble(w(i)), dimag(w(i))
     endif 
  enddo
  close(20)
end program vsi
 
real*8 function logrho_z(zhat)
  use global
  real*8, intent(in) :: zhat
  
  logrho_z = -(1d0/smallh**2d0)*(1d0 - 1d0/sqrt(1d0 + smallh**2d0*zhat**2d0))
end function logrho_z

real*8 function dlogrho_dz(zhat)
  use global
  real*8, intent(in) :: zhat
  
  dlogrho_dz = -zhat*(1d0 + smallh**2d0*zhat**2d0)**(-1.5d0)
end function dlogrho_dz

real*8 function omega2_z(zhat)
  use global
  real*8, intent(in) :: zhat

  omega2_z = 1d0 + (smallp+smallq)*smallh**2d0
  omega2_z = omega2_z + smallq*(1d0 - 1d0/sqrt(1d0+smallh**2d0*zhat**2d0))
end function omega2_z

real*8 function domega2_dz(zhat)
  use global
  real*8, intent(in) :: zhat
  
  domega2_dz = smallq*smallh**2d0*zhat
  domega2_dz = domega2_dz*(1d0 + smallh**2d0*zhat**2d0)**(-1.5d0)
end function domega2_dz

real*8 function shear_z(zhat)
  use global
  real*8, intent(in) :: zhat
  
  shear_z = -1.5d0*(1d0 + (2d0-smallq)*(smallp+smallq)*smallh**2d0/3d0 + &
       smallq*(1d0 - (1d0 + (2d0/3d0)*smallh**2d0*zhat**2d0)*(1d0 + smallh**2d0*zhat**2d0)**(-1.5d0))) 
  shear_z = shear_z/sqrt(1d0 + (smallp + smallq)*smallh**2d0 + &
       smallq*(1d0 - 1d0/sqrt(1d0 + smallh**2d0*zhat**2d0)))
end function shear_z

  
subroutine chebyshev_poly(l, zbar, T_l, dT_l, d2T_l)
  implicit none
  real*8, intent(in) :: l, zbar
  real*8, intent(out):: T_l, dT_l, d2T_l
  real*8 :: t, lt, lsq
  lsq = l*l
  t = acos(zbar)
  lt = l*t
  T_l = cos(lt)
  if(abs(zbar).lt.1d0) then
     dT_l = l*sin(lt)/sin(t)
     d2T_l= -lsq*cos(lt)/sin(t)**2 + l*cos(t)*sin(lt)/sin(t)**3
  else
     dT_l =lsq
     d2t_l =lsq*(lsq-1d0)/3d0
     if(zbar.eq.-1d0) then
        dT_l = (-1d0)**(l+1d0)*dT_l
        d2T_l= (-1d0)**(l+2d0)*d2T_l
     endif
  endif
  return
end subroutine chebyshev_poly


  
  
