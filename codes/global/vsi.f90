module global
  character*6 :: vbc
  integer :: nz, nzeff, bignz
  real*8, parameter :: pi = 2d0*acos(0d0)
  real*8 :: smallq, smallp, smallh, kx, zmax 
  real*8, allocatable :: zaxis(:), logrho(:), dlogrho(:), d2logrho(:), omega2(:), domega2(:), kappa2(:), freq(:), growth(:)
  complex*16, parameter :: ii = (0d0,1d0)
  complex*16, allocatable :: bigW(:) 
end module global

program vsi
  use global
  implicit none
  character*1, parameter :: JOBVL = 'N', JOBVR = 'V'
  integer, allocatable :: ipiv(:), eigen_test(:)
  integer :: i,ip1, j,k, loc(1), lmax, nzmid, lwork, info
  real*8 :: T_l, dT_l, d2T_l, T_lp2, dT_lp2, d2T_lp2, m, zbar,lmode
  real*8, allocatable :: rwork(:)
  complex*16, allocatable :: T(:,:), Tp(:,:), Tpp(:,:)
  complex*16, allocatable :: L1(:,:), L2(:,:), L3(:,:)
  complex*16,allocatable :: work(:), matrix(:,:), w(:), vl(:,:), &
       vr(:,:), rhs(:,:)
  real*8, external :: dlogrho_dz, domega2_dz, shear_z, omega2_z, logrho_z, kappa2_z, d2logrho_dz2
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

  if(vbc.ne.'nolagp') then
  nzeff = nz-2
  else
  nzeff = nz
  endif 

  bignz = 2*nzeff

  !allocate grids
  allocate(zaxis(nz))
  allocate(logrho(nz))
  allocate(dlogrho(nz))
  allocate(d2logrho(nz))
  allocate(omega2(nz))
  allocate(domega2(nz))
  allocate(kappa2(nz))

  allocate(bigW(nz))
  
  allocate(T(nz,nzeff))
  allocate(Tp(nz,nzeff))
  allocate(Tpp(nz,nzeff))

  !sub-matrices for linear operators
  allocate(L1(nzeff,nzeff))
  allocate(L2(nzeff,nzeff))
  allocate(L3(nzeff,nzeff))

  !big matrix and work matrices
  allocate(matrix(bignz,bignz))
  allocate(rhs(bignz,bignz))
  allocate(vr(bignz,bignz))
  allocate(vl(bignz,bignz))
  allocate(w(bignz))
  allocate(growth(bignz))
  allocate(freq(bignz))
  allocate(eigen_test(bignz))

  lwork = 8*bignz
  allocate(work(lwork))
  allocate(rwork(2*bignz))
  allocate(ipiv(bignz))
  
  !setup z axis.
  lmax = nz-1
  zaxis(Nzmid) = 0d0
  do j=Nzmid+1, nz
     zaxis(j) = -zmax*cos(pi*(j-1d0)/lmax)
  enddo
  do j=1, Nzmid-1
     zaxis(j) = -zaxis(Nz - j + 1)
  enddo
  
  !setup basic state 
  do i=1, nz
     logrho(i)  =  logrho_z(zaxis(i))
     dlogrho(i) = dlogrho_dz(zaxis(i))
     d2logrho(i)= d2logrho_dz2(zaxis(i))
     omega2(i)  = omega2_z(zaxis(i))
     if(omega2(i).le.0d0) then
     print*, 'omega2<=0'
     stop
     endif
     domega2(i) = domega2_dz(zaxis(i))
     kappa2(i)  = kappa2_z(zaxis(i))
  enddo
  
  !output basic state
  open(10,file='basic.dat')
  do i=1, nz
     write(10,fmt='(7(e22.15,x))'), zaxis(i), logrho(i), dlogrho(i), d2logrho(i),omega2(i), domega2(i), kappa2(i)
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
  if (vbc.eq.'nolagp') then
     do j=1, nz !jth physical grid 
        zbar = zaxis(j)/zmax
        do k=1, nzeff !kth basis (only T_0 ... T_lmax-2 used)
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

  !set up sub-matrices for linear operators
  !only do interior points
  !remember basis matrices size (nz, nzeff) but operators size (nzeff, nzeff) with nzeff = nz-2
  do i=1, nzeff
     
     if(vbc.ne.'nolagp') then
     ip1 = i+1 !for nograd or nopert vbc (boundary points not used in matrix)
     else
     ip1 = i   !for nolagp vbc (impost boundary condition at end points)
     endif

     L1(i,:) = kappa2(ip1)*Tpp(ip1,:) + (kappa2(ip1)*dlogrho(ip1) - ii*kx*domega2(ip1)/smallh**2d0)*Tp(ip1,:)
     L2(i,:) = smallh**2d0*Tpp(ip1,:) + smallh**2d0*dlogrho(ip1)*Tp(ip1,:) - (kx*kx+smallh**2d0*kappa2(ip1))*T(ip1,:) 
     L3(i,:) = smallh**4d0*T(ip1,:)
  enddo
  
  !explicit boundary condition for 'nolagp'
  if (vbc.eq.'nolagp') then
  do i=1,nzeff,nzeff-1
  L1(i,:) = -Tp(i,:)*dlogrho(i)
  L2(i,:) =  smallh**2d0*T(i,:)
  L3(i,:) = 0d0 
  enddo
  endif 

  !set up big matrices
  matrix(:,:) = 0d0
  rhs(:,:)    = 0d0 

  matrix(1:nzeff,1:nzeff)                       = L1
  if (vbc.ne.'nolagp') then  
  matrix(nzeff+1:bignz,nzeff+1:bignz)           = T(2:nz-1,:)
  else
  matrix(nzeff+1:bignz,nzeff+1:bignz)           = T(:,:)
  endif

  rhs(1:nzeff,1:nzeff)                          = L2
  rhs(1:nzeff,nzeff+1:bignz)                    = L3
  if (vbc.ne.'nolagp') then 
  rhs(nzeff+1:bignz,1:nzeff)                    = T(2:nz-1,:)
  else
  rhs(nzeff+1:bignz,1:nzeff)                    = T(:,:)  
  endif 

  !invert rhs and get final matrix
  call zgetrf(bignz, bignz, rhs, bignz, IPIV, INFO)
!  print*, 'inversion success?', info
  call ZGETRI(bignz, rhs, bignz, IPIV, WORK, LWORK, INFO )
!  print*, 'inversion success?', info
  matrix = matmul(rhs,matrix)

  !eigenvalue problem
  call zgeev (JOBVL, JOBVR, bignz, matrix, bignz, W, vl, bignz, VR, bignz, WORK, LWORK, RWORK, INFO) 
  print*, 'eigen success?',info

  w = sqrt(w)
  freq = dble(w)
  growth=dimag(w)
  
  open(20,file='eigenvalues.dat')
  open(30,file='eigenvectors.dat')
  do i=1, bignz
     !discard modes that grow too slowly
     if((abs(w(i)).lt.1d0/smallh).and.(abs(growth(i)).ge.smallh)) then
        write(20,fmt='(2(e22.15,x))'), freq(i), growth(i)
        bigW = matmul(T,vr(1:nzeff,i))
        do j=1,nz
           write(30,fmt='(3(e22.15,x))') zaxis(j), dble(bigW(j)), dimag(bigW(j))
        enddo
     endif
  enddo
  close(20)
  
  
end program vsi

real*8 function logrho_z(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  
  logrho_z = -(1d0/smallh**2d0)*(1d0 - 1d0/sqrt(1d0 + smallh**2d0*zhat**2d0))
end function logrho_z

real*8 function dlogrho_dz(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  
  dlogrho_dz = -zhat*(1d0 + smallh**2d0*zhat**2d0)**(-1.5d0)
end function dlogrho_dz

real*8 function d2logrho_dz2(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  
  d2logrho_dz2 =(2d0*smallh**2d0*zhat**2d0-1d0)*(1d0+smallh**2d0*zhat**2d0)**(-2.5d0)
end function d2logrho_dz2

real*8 function omega2_z(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  
  omega2_z = 1d0 + (smallp+smallq)*smallh**2d0
  omega2_z = omega2_z + smallq*(1d0 - 1d0/sqrt(1d0+smallh**2d0*zhat**2d0))
end function omega2_z

real*8 function domega2_dz(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  
  domega2_dz = smallq*smallh**2d0*zhat
  domega2_dz = domega2_dz*(1d0 + smallh**2d0*zhat**2d0)**(-1.5d0)
end function domega2_dz

real*8 function shear_z(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  
  shear_z = -1.5d0*(1d0 + (2d0-smallq)*(smallp+smallq)*smallh**2d0/3d0 + &
       smallq*(1d0 - (1d0 + 2d0*smallh**2d0*zhat**2d0/3d0)*(1d0 + smallh**2d0*zhat**2d0)**(-1.5d0)))
  shear_z = shear_z/sqrt(1d0 + (smallp + smallq)*smallh**2d0 + &
       smallq*(1d0 - 1d0/sqrt(1d0 + smallh**2d0*zhat**2d0)))
end function shear_z

real*8 function kappa2_z(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  real*8 :: omega
  real*8, external :: omega2_z, shear_z
  
  omega = sqrt(omega2_z(zhat))
  
  kappa2_z = 2d0*omega*(2d0*omega + shear_z(zhat))
end function kappa2_z



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


  
  
