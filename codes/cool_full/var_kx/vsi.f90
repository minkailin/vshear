module global
  character*6 :: vbc
  logical :: gterms
  integer :: nz, nzeff, bignz, nkx
  real*8, parameter :: pi = 2d0*acos(0d0)
  real*8 :: smallq, smallp, eps, smallh, kx, zmax, gmma, bgmma, bcool, smalls
  real*8 :: kxmin, kxmax, dlogkx
  real*8, allocatable :: zaxis(:), logrho(:), dlogrho(:), d2logrho(:), omega2(:), domega2(:), kappa2(:), csq(:), freq(:), growth(:)
  real*8, allocatable :: kaxis(:)
  complex*16, parameter :: ii = (0d0,1d0)
  complex*16, allocatable :: bigW(:), dbigW(:), bigQ(:), vx(:), vy(:), vz(:), dvz(:) 
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
       L5(:,:), L6(:,:), L7(:,:), L8(:,:), &
       L9(:,:), L10(:,:), L11(:,:), &
       L1bar(:,:), L2bar(:,:), L3bar(:,:), L4bar(:,:), L5bar(:,:), &
       L6bar(:,:)
  complex*16,allocatable :: work(:), matrix(:,:), w(:), vl(:,:), &
       vr(:,:), rhs(:,:)
  real*8, external :: dlogrho_dz, domega2_dz, omega2_z, logrho_z, kappa2_z, d2logrho_dz2, csq_z 
  namelist /params/ smallq, smallp, eps, gmma, bgmma, bcool, zmax, vbc, gterms
  namelist /loop/ kxmin, kxmax, nkx 
  namelist /grid/ nz
  
  !read input parameters.
  open(7, file="input")
  read(7, nml=params)
  read(7, nml=loop)
  read(7, nml=grid)
  close(7)

  if(mod(nz,2).eq.0) then
     print*, 'Nz needs to be odd but Nz=', nz
     stop
  else
     Nzmid = (Nz+1)/2
  endif

  if(gterms.eqv..true.) print*, 'adding extra terms due to rad entrop grad'

  smallh = (1d0 - eps**2d0/(bgmma-1d0))**(-2d0) - 1d0
  smallh = sqrt(smallh)  

  smalls = smallq + smallp*(1d0 - bgmma) 

  print*, 'smallh, smalls=', smallh, smalls 

  open(10,file='params.dat') 
  write(10,fmt='(7(e22.15,x))'), gmma, bgmma, eps, smallh, smalls, bcool, kx 
  close(10)

  if((vbc.ne.'nolagp').and.(vbc.ne.'nozvel')) then
  print*, 'only valid bc is "nolagp" or "nozvel"'
  stop
  endif 

  nzeff = nz
  bignz = 5*nzeff

  !allocate grids
  allocate(kaxis(nkx))

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
  allocate(vx(nz))
  allocate(vy(nz))
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
  allocate(L5(nzeff,nzeff))
  allocate(L6(nzeff,nzeff))
  allocate(L7(nzeff,nzeff))
  allocate(L8(nzeff,nzeff))
  allocate(L9(nzeff,nzeff))
  allocate(L10(nzeff,nzeff))
  allocate(L11(nzeff,nzeff))

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
  
  !setup kaxis. input is kxHiso. then output
  dlogkx = log10(kxmax/kxmin)/(nkx-1d0)
  do i=1, nkx
  kaxis(i) =10**(log10(kxmin) + dlogkx*(i-1d0))
  enddo
  open(10,file='kaxis.dat')
  do i=1, nkx
     write(10,fmt='(e22.15,x)'), kaxis(i)
  enddo
  close(10)

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

  if ((vbc.eq.'nolagp').or.(vbc.eq.'nozvel')) then !vanishing lagragian pressure pert, vanishing vertical vel 
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

  open(20,file='eigenvalues.dat')
  open(30,file='eigenvectors.dat')
  
  do k=1, nkx !loop over kx values
     kx = kaxis(k)*smallh !convert input into code units 
     print*, 'iteration, kxHiso=', k, kaxis(k)
     
     !set up sub-matrices for linear operators
     do i=1, nzeff
        L1(i,:) = -T(i,:)/bcool
        L2(i,:) =  T(i,:)/bcool/bgmma 
        L3(i,:) =  csq(i)*(eps/smallh)*dlogrho(i)*(gmma/bgmma-1d0)*T(i,:) 
        
        L4(i,:) = ii*kx*T(i,:)/eps 
        L5(i,:) = Tp(i,:) + dlogrho(i)*T(i,:) 
        
        L6(i,:) = ii*kx*T(i,:)/smallh 
        L7(i,:) = -2d0*sqrt(omega2(i))*T(i,:)  
        
        L8(i,:) = kappa2(i)*T(i,:)/(2d0*sqrt(omega2(i))) 
        L9(i,:) = domega2(i)*T(i,:)/(2d0*smallh*sqrt(omega2(i)))
        
        L10(i,:) = Tp(i,:) + dlogrho(i)*T(i,:)
        L11(i,:) = -dlogrho(i)*T(i,:) 
        
        
        L1bar(i,:) = -ii*eps*T(i,:)
        L2bar(i,:) =  ii*eps*(gmma/bgmma)*T(i,:) 
        L3bar(i,:) =  ii*smallh*T(i,:)/csq(i)
        L4bar(i,:) =  ii*eps*T(i,:) 
        L5bar(i,:) =  ii*eps*T(i,:) 
        L6bar(i,:) =  ii*smallh*T(i,:)
     enddo
     
     !set up big matrices
     matrix(:,:) = 0d0
     rhs(:,:)    = 0d0 
     
     matrix(1:nzeff,1:nzeff)                       = L1
     matrix(1:nzeff,nzeff+1:2*nzeff)               = L2
     if(gterms.eqv..true.) then
        do i=1,nzeff !extra term for global radial entropy gradient 
           matrix(i,2*nzeff+1:3*nzeff) =  -eps*csq(i)*smalls*T(i,:)/bgmma  
        enddo
     endif
     matrix(1:nzeff,4*nzeff+1:bignz)               = L3 
     
     matrix(nzeff+1:2*nzeff, 2*nzeff+1:3*nzeff)    = L4 
     matrix(nzeff+1:2*nzeff, 4*nzeff+1:bignz)      = L5
     
     matrix(2*nzeff+1:3*nzeff,1:nzeff)             = L6
     if(gterms.eqv..true.) matrix(2*nzeff+1:3*nzeff,nzeff+1:2 *nzeff)     =smalls*T(:,:)*eps/bgmma !extra term for global radial entropy gradient 
     matrix(2*nzeff+1:3*nzeff,3*nzeff+1:4*nzeff)   = L7
     
     
     matrix(3*nzeff+1:4*nzeff,2*nzeff+1:3*nzeff)   = L8
     matrix(3*nzeff+1:4*nzeff,4*nzeff+1:bignz)     = L9
     
     matrix(4*nzeff+1:bignz,1:nzeff)               = L10
     matrix(4*nzeff+1:bignz,nzeff+1:2*nzeff)       = L11
     
     rhs(1:nzeff,1:nzeff)                          = L1bar
     rhs(1:nzeff,nzeff+1:2*nzeff)                  = L2bar 
     
     rhs(nzeff+1:2*nzeff,nzeff+1:2*nzeff)          = L3bar 
     
     rhs(2*nzeff+1:3*nzeff,2*nzeff+1:3*nzeff)      = L4bar
     
     rhs(3*nzeff+1:4*nzeff,3*nzeff+1:4*nzeff)      = L5bar 
     
     rhs(4*nzeff+1:bignz,4*nzeff+1:bignz)           = L6bar
     
     !explicit boundary condition for 'nolagp'
     if (vbc.eq.'nolagp') then
        do i=1,nzeff,nzeff-1
           matrix(i,:) = 0d0 
           matrix(i,4*nzeff+1:bignz) = csq(i)*dlogrho(i)*T(i,:) 
           
           rhs(i,:) = 0d0 
           rhs(i,1:nzeff) = ii*smallh*T(i,:) 
        enddo
     endif
     
     !explicit bc for vz=0 
     if (vbc.eq.'nozvel') then
        do i=1,nzeff,nzeff-1
           matrix(i,:) = 0d0
           rhs(i,:) = 0d0
           
           !set vz =0 in energy eq 
           matrix(i,1:nzeff) = T(i,:)
           matrix(i,nzeff+1:2*nzeff) = - T(i,:)/bgmma 
           rhs(i,1:nzeff) = ii*eps*T(i,:)*bcool
           rhs(i,nzeff+1:2*nzeff) = - ii*eps*T(i,:)*gmma*bcool/bgmma
           
        enddo
     endif
     
     
     !invert rhs and get final matrix
     call zgetrf(bignz, bignz, rhs, bignz, IPIV, INFO)
     if(info.ne.0) print*, 'inversion fail'
     call ZGETRI(bignz, rhs, bignz, IPIV, WORK, LWORK, INFO )
     if(info.ne.0) print*, 'inversion fail'
     matrix = matmul(rhs,matrix)
     
     !eigenvalue problem
     call zgeev (JOBVL, JOBVR, bignz, matrix, bignz, W, vl, bignz, VR, bignz, WORK, LWORK, RWORK, INFO) 
     print*, 'eigen success?',info
     
     freq = dble(w)
     growth=dimag(w)

     !discard modes with eigenvalues too large or too small
     !discard modes that decay
     !discard modes with too small imaginary part
     do i=1, bignz
        if((abs(w(i)).gt.1d0/eps).or.(growth(i).lt.0d0).or.(abs(w(i)).lt.eps**2d0)) then
           w(i) = (1d6,1d6)
        endif
     enddo

     freq = dble(w)
     growth=dimag(w)
     
     !now pick the mode with smallest |real freq| (fundamental mode)
     loc = minloc(abs(freq))
     i = loc(1)
     write(20,fmt='(2(e22.15,x))'), freq(i), growth(i)
     bigW = matmul(T,vr(1:nzeff,i))
     dbigW= matmul(Tp,vr(1:nzeff,i))
     
     bigQ = matmul(T,vr(nzeff+1:2*nzeff,i))
     
     vx   = matmul(T,vr(2*nzeff+1:3*nzeff,i))
     vy   = matmul(T,vr(3*nzeff+1:4*nzeff,i))
     vz   = matmul(T,vr(4*nzeff+1:bignz,i))
     dvz  = matmul(Tp,vr(4*nzeff+1:bignz,i)) 
     do j=1,nz
        write(30,fmt='(14(e22.15,x))') dble(bigW(j)), dimag(bigW(j)), dble(dbigW(j)), dimag(dbigW(j)), dble(bigQ(j)), dimag(bigQ(j)), & 
             dble(vx(j)), dimag(vx(j)), dble(vy(j)), dimag(vy(j)), dble(vz(j)), dimag(vz(j)), dble(dvz(j)), dimag(dvz(j)) 
     enddo
  enddo
  
  
  close(20)
  close(30)

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


  
  
