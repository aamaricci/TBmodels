! DESCRIPTION
!   Solve the non-interacting 2x2 part of the BHZ model and generate the hamiltoniana matrix H(k), 
!   Model parameters are: t_dd,v0/tpd,ep0,ed0.
!   Using the 2d square lattice dispersion. 
!   The matrix is written in the Wannier90 form, as expected by w2CTQMC & ggED code.
!   Output is on file *hkfile.in default.
!
! MODEL Hamiltonian is:
!|m-(Cos{kx}+Cos{ky})         & \lambda*(Sin{kx}-i*Sin{ky}) &               0            &           0            |
!|\lambda*(Sin{kx}+i*Sin{ky}) & -m+(Cos{kx}+Cos{ky})        &               0            &           0            |
!|             0              &          0                  & m-(Cos{kx}+Cos{ky}) & \lambda*(-Sin{kx}-i*Sin{ky})  |
!|             0              &          0                  & \lambda *(-Sin{kx}+i*Sin{ky}) & -m+(Cos{kx}+Cos{ky})|
!=
!|    h^{2x2}(k)      &          0      |
!|        0           &  [h^{2x2}]*(-k) |
!
! OPTIONS
!   nkx=[100]  - Number of points along the x component of the k-grid
!   nky=[100]  - Number of points along the y component of the k-grid
!   m=[1.0]    - 
!   lambda=[0.0]- 
!   beta=[100] - Inverse temperature 
!   file=[hkfile_bhz.in]- output unit file,default is no write
!   dcflag=[T] - output unit file,default is no write

program tdd_pam
  USE SCIFOR
  implicit none

  integer,parameter    :: L=1000,Norb=2
  integer              :: i,j,k,ik,iorb,jorb
  real(8)              :: mh,rh,lambda,delta
  real(8)              :: xmu,beta,eps
  integer              :: dcshift,count
  real(8)              :: epsik,fmesh,n(Norb),meig(Norb)
  integer              :: Nkx,Nky,Nk
  real(8)              :: ix,iy
  real(8)              :: kx,ky
  real(8),dimension(L) :: wm,wr
  complex(8)           :: w,Hk(Norb,Norb)
  complex(8)           :: fg(L,Norb,Norb),fgr(L,Norb,Norb)
  character(len=20)    :: file,nkstring
  logical              :: iexist,ibool,dcflag

  namelist/hkvars/nkx,nky,mh,lambda,delta,Rh,xmu,beta,eps,file,dcflag

  nkx=100
  nky=100
  mh = 3.d0
  lambda=0.3d0
  delta=0.d0
  Rh=0.d0
  xmu=0.d0
  eps=4.d-2
  beta=100.d0
  file="hkfile_bhz.in"
  dcflag=.true.

  inquire(file="inputBHZ.in",exist=iexist)
  if(iexist)then
     open(10,file="inputBHZ.in")
     read(10,nml=hkvars)
     close(10)
  endif

  call parse_cmd_variable(nkx,"NKX")
  call parse_cmd_variable(nky,"NKY")
  call parse_cmd_variable(mh,"MH")
  call parse_cmd_variable(rh,"RH")
  call parse_cmd_variable(lambda,"LAMBDA")
  call parse_cmd_variable(delta,"DELTA")
  call parse_cmd_variable(xmu,"XMU")
  call parse_cmd_variable(eps,"EPS")
  call parse_cmd_variable(beta,"BETA")
  call parse_cmd_variable(file,"FILE")
  call parse_cmd_variable(dcflag,"DCFLAG")

  dcshift=0 ; if(dcflag)dcshift=1

  write(*,nml=hkvars)

  wm = pi/beta*real(2*arange(1,L)-1,8)

  wr = linspace(-10.d0,10.d0,L,mesh=fmesh)

  Nk=Nkx*Nky
  open(50,file=trim(file))
  write(nkstring,*)Nk
  write(50,*)trim(adjustl(trim(Nkstring)))," 1 2 1 1"

  fgr=zero ; fg =zero ;count=0
  call start_progress
  do ix=1,Nkx
     kx = -pi + 2.d0*pi*real(ix-1,8)/real(Nkx,8)
     do iy=1,Nky
        ky = -pi + 2.d0*pi*real(iy-1,8)/real(Nky,8)
        Hk(:,:) = hk_bhz2x2(kx,ky)
        write(50,"(3(F10.7,1x))")kx,ky,pi
        do iorb=1,Norb
           write(50,"(10(2F10.7,1x))")(Hk(iorb,jorb),jorb=1,Norb)
        enddo

        do i=1,L
           w = cmplx(wr(i),eps,8)+xmu
           fgr(i,:,:)=fgr(i,:,:) + inverse_g0k2x2(w,Hk)
           w = xi*wm(i)+xmu
           fg(i,:,:) =fg(i,:,:)  + inverse_g0k2x2(w,Hk)
        enddo
        count=count+1
        call progress_bar(count,Nk)
     enddo
  enddo
  call stop_progress
  fgr= fgr/real(Nk,8)
  fg = fg/real(Nk,8)


  ik = 0
  open(10,file="Eigenbands.dat")
  !From \Gamma=(0,0) to X=(pi,0): 100 steps
  do ix=1,100
     ik=ik+1
     kx = 0.d0 + pi*real(ix-1,8)/100.d0
     ky = 0.d0
     Hk(:,:) = hk_bhz2x2(kx,ky)
     meig = Eigk2x2(Hk)
     write(10,"(I0,4F25.12)")ik,meig(1),meig(2)
  enddo
  !From X=(pi,0) to M=(pi,pi): 100 steps
  do iy=1,100
     ik=ik+1
     kx = pi
     ky = 0.d0 + pi*real(iy-1,8)/100.d0
     Hk(:,:) = hk_bhz2x2(kx,ky)
     meig = Eigk2x2(Hk)
     write(10,"(I0,4F25.12)")ik,meig(1),meig(2)
  enddo
  !From M=(pi,pi) to \Gamma=(0,0): 100 steps
  do ix=1,100
     ik=ik+1
     iy=ix
     kx = pi - pi*real(ix-1,8)/100.d0
     ky = pi - pi*real(iy-1,8)/100.d0
     Hk(:,:) = hk_bhz2x2(kx,ky)
     meig = Eigk2x2(Hk)
     write(10,"(I0,4F25.12)")ik,meig(1),meig(2)
  enddo
  close(10)





contains


  function hk_bhz2x2(kx,ky) result(hk)
    real(8)                   :: kx,ky,epsik
    complex(8),dimension(2,2) :: hk
    epsik   = cos(kx)+cos(ky)
    hk(1,1) = mh - epsik
    hk(2,2) =-mh + epsik
    hk(1,2) = lambda*(sin(kx)-xi*sin(ky))
    hk(2,1) = lambda*(sin(kx)+xi*sin(ky))
  end function hk_bhz2x2


  function inverse_g0k2x2(iw,hk) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: delta,ppi,vmix
    g0k=zero
    delta = iw - hk(1,1)
    ppi   = iw - hk(2,2)
    vmix  = -hk(1,2)
    g0k(1,1) = one/(delta - abs(vmix)**2/ppi)
    g0k(2,2) = one/(ppi - abs(vmix)**2/delta)
    g0k(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    g0k(2,1) = conjg(g0k(1,2))
  end function inverse_g0k2x2


  function Eigk2x2(hk) result(eig)
    complex(8),dimension(2,2) :: hk
    real(8),dimension(2)      :: eig
    ! call eigh(hk,eig)
    eig(1)=hk(1,1)+hk(2,2) + sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eig(2)=hk(1,1)+hk(2,2) - sqrt(abs(hk(1,1)-hk(2,2))**2 + 4.d0*hk(1,2)*hk(2,1) )
    eig = eig/2.d0
  end function Eigk2x2


  function matrix_uk(hk) result(uk)
    complex(8),dimension(2,2) :: hk
    complex(8),dimension(2,2) :: uk
    complex(8) :: delta0,ek,u,v
    delta0 = -(hk(1,1)-hk(2,2))
    ek     = sqrt( delta0**2 + 4.d0*hk(1,2)*hk(2,1) )
    u = sqrt(0.5d0*(1.d0 - delta0/ek))
    v = sqrt(0.5d0*(1.d0 + delta0/ek))
    uk(1,1) = u
    uk(2,2) = u
    uk(1,2) = v
    uk(2,1) =-v
  end function matrix_uk


end program tdd_pam


