module global
  character*6 :: vbc
  integer :: nz, nzeff, bignz
  real*8, parameter :: pi = 2d0*acos(0d0)
  real*8 :: smallq, smallp, eps, smallh, kx, zmax, gmma, bgmma, bcool, smalls  
  real*8, allocatable :: zaxis(:), logrho(:), dlogrho(:), d2logrho(:), omega2(:), domega2(:), kappa2(:), csq(:), freq(:), growth(:)
  complex*16, parameter :: ii = (0d0,1d0)
  complex*16, allocatable :: bigW(:), dbigW(:), bigQ(:), dbigQ(:), d2bigW(:), vz(:), dvz(:)
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
  complex*16, allocatable :: L1(:,:), L2(:,:), L3(:,:), L4(:,:), &
       L1bar(:,:), L2bar(:,:), L3bar(:,:), L4bar(:,:), L5bar(:,:), &
       L6bar(:,:)
  complex*16,allocatable :: work(:), matrix(:,:), w(:), vl(:,:), &
       vr(:,:), rhs(:,:)
  real*8, external :: dlogrho_dz, domega2_dz, omega2_z, logrho_z, kappa2_z, d2logrho_dz2, csq_z 
  namelist /params/ smallq, smallp, eps, gmma, bgmma, bcool, kx, zmax, vbc
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

  smallh = 1d0/(1d0 - eps**2d0/(bgmma-1d0))**2d0 - 1d0
  smallh = sqrt(smallh)  

  smalls = smallq + smallp*(1d0 - bgmma) 

  !convert input kx to code units
  kx = kx*2d0*pi*eps*smallh 

  print*, 'smallh, smalls=', smallh, smalls 

  open(10,file='params.dat') 
  write(10,fmt='(7(e22.15,x))'), gmma, bgmma, eps, smallh, smalls, bcool, kx 
  close(10)

!!$  if((vbc.ne.'nolagp').and.(vbc.ne.'nozvel')) then
!!$  print*, 'only valid bc is "nolagp" or "nozvel"'
!!$  stop
!!$  endif 

  nzeff = nz

  bignz = 4*nzeff

  !allocate grids
  allocate(zaxis(nz))
  allocate(logrho(nz))
  allocate(dlogrho(nz))
  allocate(d2logrho(nz))
  allocate(omega2(nz))
  allocate(domega2(nz))
  allocate(kappa2(nz))
  allocate(csq(nz))

  allocate(bigQ(nz))
  allocate(bigW(nz))
  allocate(dbigW(nz))  
  allocate(dbigQ(nz))
  allocate(d2bigW(nz))
  allocate(vz(nz))
  allocate(dvz(nz))

  allocate(T(nz,nzeff))
  allocate(Tp(nz,nzeff))
  allocate(Tpp(nz,nzeff))

  !sub-matrices for linear operators
  allocate(L1(nzeff,nzeff))
  allocate(L2(nzeff,nzeff))
  allocate(L3(nzeff,nzeff))
  allocate(L4(nzeff,nzeff))

  allocate(L1bar(nzeff,nzeff))
  allocate(L2bar(nzeff,nzeff))
  allocate(L3bar(nzeff,nzeff))
  allocate(L4bar(nzeff,nzeff))
  allocate(L5bar(nzeff,nzeff))
  allocate(L6bar(nzeff,nzeff))

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
     csq(i)     = csq_z(zaxis(i))
  enddo
  
  !output basic state
  open(10,file='basic.dat')
  do i=1, nz
     write(10,fmt='(8(e22.15,x))'), zaxis(i), logrho(i), dlogrho(i), d2logrho(i),omega2(i), domega2(i), kappa2(i), csq(i)
  enddo
  close(10)
  
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

  !set up sub-matrices for linear operators
  
  do i=1, nzeff
     L1(i,:) = Tpp(i,:) + (2d0*dlogrho(i) - ii*kx*domega2(i)/smallh/eps/kappa2(i))*Tp(i,:) &
          +((dlogrho(i) - ii*kx*domega2(i)/smallh/eps/kappa2(i))*dlogrho(i) + d2logrho(i))*T(i,:) 
     L1(i,:) = -L1(i,:)

     L2(i,:) = (d2logrho(i) + (dlogrho(i) - ii*kx*domega2(i)/smallh/eps/kappa2(i))*dlogrho(i))*T(i,:) &
          +dlogrho(i)*Tp(i,:)
     L3(i,:) = (kx*kx/kappa2(i))*T(i,:)
     L4(i,:) = smallh**2d0*T(i,:)/csq(i) 

!     L1bar(i,:)= csq(i)*(eps/smallh)**2d0*(gmma/bgmma-1d0)*dlogrho(i)*(Tp(i,:) + dlogrho(i)*T(i,:))
!     L2bar(i,:)=-csq(i)*(eps/smallh)**2d0*(gmma/bgmma-1d0)*dlogrho(i)**2d0*T(i,:)
!     L3bar(i,:)= eps**2d0*T(i,:)
!     L4bar(i,:)=-eps**2d0*gmma*T(i,:)/bgmma
!     L5bar(i,:)= eps*ii*T(i,:)/bcool
!     L6bar(i,:)=-eps*ii*T(i,:)/bcool/bgmma

     L1bar(i,:)= csq(i)*eps*(1d0/smallh)**2d0*(gmma/bgmma-1d0)*dlogrho(i)*(Tp(i,:) + dlogrho(i)*T(i,:))
     L2bar(i,:)=-csq(i)*eps*(1d0/smallh)**2d0*(gmma/bgmma-1d0)*dlogrho(i)**2d0*T(i,:)
     L3bar(i,:)= eps*T(i,:)
     L4bar(i,:)=-eps*gmma*T(i,:)/bgmma
     L5bar(i,:)= ii*T(i,:)/bcool
     L6bar(i,:)=-ii*T(i,:)/bcool/bgmma

  enddo

  !set up big matrices
  matrix(:,:) = 0d0
  rhs(:,:)    = 0d0 

  matrix(1:nzeff,1:nzeff)                       = L1
  matrix(1:nzeff,nzeff+1:2*nzeff)               = L2
  matrix(nzeff+1:2*nzeff,1:nzeff)               = L1bar
  matrix(nzeff+1:2*nzeff,nzeff+1:2*nzeff)       = L2bar
  matrix(2*nzeff+1:3*nzeff,2*nzeff+1:3*nzeff)   = T
  matrix(3*nzeff+1:bignz,3*nzeff+1:bignz)       = T

  rhs(1:nzeff,2*nzeff+1:3*nzeff)                = L3
  rhs(1:nzeff,3*nzeff+1:bignz)                  = L4
  rhs(nzeff+1:2*nzeff,1:nzeff)                  = L5bar
  rhs(nzeff+1:2*nzeff,nzeff+1:2*nzeff)          = L6bar
  rhs(nzeff+1:2*nzeff,2*nzeff+1:3*nzeff)        = L3bar
  rhs(nzeff+1:2*nzeff,3*nzeff+1:bignz)          = L4bar
  rhs(2*nzeff+1:3*nzeff,1:nzeff)                = T
  rhs(3*nzeff+1:bignz,nzeff+1:2*nzeff)          = T

  !explicit boundary condition for 'nolagp'
  if (vbc.eq.'nolagp') then
     do i=1,nzeff,nzeff-1
        matrix(i,:) = 0d0 
        matrix(i,1:nzeff) = -dlogrho(i)*csq(i)*(Tp(i,:) + dlogrho(i)*T(i,:))
        matrix(i,nzeff+1:2*nzeff) = csq(i)*dlogrho(i)**2d0*T(i,:)
        
        rhs(i,:) = 0d0 
        rhs(i,2*nzeff+1:3*nzeff) = smallh**2d0*T(i,:)   
     enddo
  endif

  !explicit bc for vz=0 (set vz=0 in energy eq)
  if (vbc.eq.'nozvel') then
     do i=1,nzeff,nzeff-1
        matrix(i,:) = 0d0
        rhs(i,:) = 0d0

! energy eq with vz=0 
!        matrix(i,1:nzeff) = T(i,:)
!        matrix(i,nzeff+1:2*nzeff) = -T(i,:)/bgmma
        
!        rhs(i,1:nzeff) = eps*ii*T(i,:)*bcool
!        rhs(i,nzeff+1:2*nzeff) = -eps*ii*gmma*T(i,:)*bcool/bgmma 


! energy eq with vz=0, then mult. by i*sigma  (all terms on rhs) 
         rhs(i,1:nzeff) = ii*T(i,:) 
         rhs(i,nzeff+1:2*nzeff) = -ii*T(i,:)/bgmma
         rhs(i,2*nzeff+1:3*nzeff)  = eps*T(i,:)*bcool 
         rhs(i,3*nzeff+1:bignz)    = -eps*T(i,:)*bcool*gmma/bgmma  
     enddo
  endif
 
 
  !invert rhs and get final matrix
  call zgetrf(bignz, bignz, rhs, bignz, IPIV, INFO)
  print*, 'inversion success?', info
  call ZGETRI(bignz, rhs, bignz, IPIV, WORK, LWORK, INFO )
  print*, 'inversion success?', info
  matrix = matmul(rhs,matrix)

  !eigenvalue problem
  call zgeev (JOBVL, JOBVR, bignz, matrix, bignz, W, vl, bignz, VR, bignz, WORK, LWORK, RWORK, INFO) 
  print*, 'eigen success?',info

  freq = dble(w)
  growth=dimag(w)
  
  eigen_test(:) = 0 

  open(20,file='eigenvalues.dat')
  open(30,file='eigenvectors.dat')
  do i=1, bignz
     !discard modes with eigenvalues too large (approx eqns assume |sigma^2|<<kappa2
     !discard modes that grow too slowly
     if((abs(w(i)).le.1d0/eps).and.(growth(i).ge.eps**2d0)) then
        eigen_test(i) = 1
        write(20,fmt='(2(e22.15,x))'), freq(i), growth(i)
        bigW = matmul(T,vr(1:nzeff,i))
        dbigW= matmul(Tp,vr(1:nzeff,i))
        d2bigW=matmul(Tpp,vr(1:nzeff,i)) 

        bigQ = matmul(T,vr(nzeff+1:2*nzeff,i))
       dbigQ = matmul(Tp,vr(nzeff+1:2*nzeff,i))

        vz = (dbigW + dlogrho*(bigW-bigQ))/(ii*smallh*w(i))
        dvz= (d2bigW + d2logrho*(bigW-bigQ) + dlogrho*(dbigW - dbigQ))/(ii*smallh*w(i))

        do j=1,nz
           write(30,fmt='(10(e22.15,x))') dble(bigW(j)), dimag(bigW(j)), dble(dbigW(j)), dimag(dbigW(j)), dble(bigQ(j)), dimag(bigQ(j)), &
                                         dble(vz(j)), dimag(vz(j)), dble(dvz(j)), dimag(dvz(j)) 
        enddo
     endif
  enddo
  if (maxval(eigen_test).eq.0) then
  print*, 'no good eigenvalues found'
  stop
  endif 
  close(20)
  

end program vsi

real*8 function logrho_z(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  real*8 :: zsq  

  zsq = (smallh*zhat)**2d0 

  logrho_z = 1d0 + ( (bgmma-1d0)/eps**2d0 )*(1d0/sqrt(1d0+zsq) - 1d0)
  logrho_z = log(logrho_z)/(bgmma - 1d0)
end function logrho_z

real*8 function csq_z(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  real*8 :: zsq
  
  zsq = (smallh*zhat)**2d0

  csq_z = 1d0 + ( (bgmma-1d0)/eps**2d0 )*(1d0/sqrt(1d0+zsq) - 1d0)
end function csq_z

real*8 function dlogrho_dz(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  real*8 :: temp

  temp = 1d0 + smallh**2d0*zhat**2d0  

  dlogrho_dz = -smallh**2d0*zhat*temp**(-1.5d0)
  dlogrho_dz = dlogrho_dz/( eps**2d0 + (bgmma-1d0)*(1d0/sqrt(temp) - 1d0) )
end function dlogrho_dz

real*8 function d2logrho_dz2(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  real*8 :: temp, fac, zsq

  zsq = (smallh*zhat)**2d0
  temp = 1d0 + zsq
  fac = eps**2d0 +(bgmma-1d0)*(1d0/sqrt(temp) - 1d0)  

  d2logrho_dz2 = (1d0 - 2d0*zsq)*fac &
                + (bgmma-1d0)*zsq/sqrt(temp)
  d2logrho_dz2 = -d2logrho_dz2*smallh**2d0*temp**(-2.5d0)
  d2logrho_dz2 = d2logrho_dz2/fac**2d0
end function d2logrho_dz2

real*8 function omega2_z(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  real*8 :: temp 
  real*8, external :: csq_z   
  
  temp = 1d0 + (smallh*zhat)**2d0

  omega2_z = eps**2d0*(smallp + smalls/bgmma) + (smalls/bgmma)*(1d0 - 1d0/sqrt(temp)) &
            + 1d0 
end function omega2_z

real*8 function domega2_dz(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  real*8 :: zsq
  
  zsq = (smallh*zhat)**2d0  

  domega2_dz = (smalls/bgmma)*smallh**2d0*zhat*(1d0 + zsq)**(-1.5d0) 
end function domega2_dz

real*8 function kappa2_z(zhat)
  use global
  implicit none
  real*8, intent(in) :: zhat
  real*8 :: temp, zsq 
  real*8 :: rdomega2_dr, rdlogrho_dr
  real*8, external :: omega2_z, csq_z 
 
  zsq   = (smallh*zhat)**2d0
  temp = 1d0 + zsq 

  rdomega2_dr = eps**2d0*(smallq - 2d0)*(smallp+smalls/bgmma) &
               -3d0*(smalls/bgmma)*(1d0 - 1d0/sqrt(temp)) &
               -(smalls/bgmma)*zsq*temp**(-1.5d0) &
               -3d0   
  kappa2_z = 4d0*omega2_z(zhat) + rdomega2_dr 
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


  
  
