program BaCoNiS2
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                 :: Norb
  integer                                 :: Nspin,Nso,Nlso
  integer                                 :: Nk,Nkpath,Nkx,Nky,Nkz,Nkpts
  !
  integer                                 :: i,j,k,ik,ikx,iky,iorb,jorb,ikpoint,Lreal,ispin,info
  real(8),dimension(2)                    :: bk1,bk2,ei1,ei2
  real(8),dimension(:,:),allocatable      :: kpath,ktrims
  !
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk,Eval
  real(8)                                 :: Ez2,Exy,Ez,Ea
  real(8)                                 :: tz22z,tz22a,txy2z,txy2a
  real(8)                                 :: t1,t2,t3,t4,Jz2,Jxy,Jpz2,Jpxy,lam
  real(8)                                 :: dEz2,dExy
  real(8)                                 :: Efermi
  real(8)                                 :: xmu,beta,eps,wmax,wmin
  namelist/params/ Ez2,Exy,Ez,Ea,tz22z,tz22a,txy2z,txy2a
  !
  character(len=64)                       :: finput
  complex(8),dimension(:,:,:),allocatable :: rhoK
  real(8),dimension(:,:),allocatable      :: eK,Kgrid
  real(8)                                 :: filling,Ef,k0,Emin,delta,E0,Eshift

  Norb = 2
  Nspin = 1 
  Nso   = Nspin*Norb
  Nlso  = 2*Nso
  !
  call parse_cmd_variable(finput,"finput",default="inputBaNiS2.conf")
  call parse_input_variable(nkx,"NKX",finput,default=70)
  call parse_input_variable(nky,"NKY",finput,default=70)
  call parse_input_variable(nkpath,"nkpath",finput,default=200)
  !
  call parse_input_variable(delta,"delta",finput,default=0d0)
  call parse_input_variable(eshift,"eshift",finput,default=0d0)
  call parse_input_variable(filling,"filling",finput,default=1d0)
  call parse_input_variable(xmu,"xmu",finput,default=0d0)
  call parse_input_variable(beta,"beta",finput,default=1000d0)
  call parse_input_variable(eps,"eps",finput,default=0.1d0)
  call parse_input_variable(wmax,"wmax",finput,default=10d0)
  call parse_input_variable(wmin,"wmin",finput,default=-10d0)
  !
  call save_input_file(reg(finput))
  !
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")
  !
  !
  Nk = Nkx*Nky

  open(unit=10, file='BaNiS2_params.in')
  read(10,NML=params)
  close(10)

  filling = 4d0*filling


  Exy = Exy - Ez2
  Ez  = Ez-Ez2-Delta
  Ea  = Ea-Ez2-Delta
  Ez2 = 0d0!! Delta
  ! Exy = Exy + Delta
  !
  ! Ez  = Ez-Delta
  ! Ea  = Ea-Delta
  write(201,*)Ez2,Exy,Ez,Ea

  call TB_set_bk([pi2,0d0],[0d0,pi2])
  call TB_print_bk

  !
  Nkpts = 3
  allocate(Kpath(Nkpts,2))
  Kpath(1,:) = [-pi,0d0] !Gamma
  Kpath(2,:) = [0d0, 0d0] !X
  Kpath(3,:) = [pi, 0d0] 
  ! Kpath(4,:) = [0d0, 0d0]
  ! Kpath(5,:) = [-pi, -pi]
  ! !

  Jz2 =  4*tz22z**2/Ez
  Jxy =  4*txy2z**2/Ez
  Jpz2 = Jz2 - 4*tz22a**2/Ea
  Jpxy = Jxy - 4*txy2a**2/Ea
  lam  = -4d0*sqrt(Jz2*Jxy)-tz22a*txy2a/Ea
  dEz2 = Ez2 - 2*Jz2 + Jpz2
  dExy = Exy - 2*Jxy + Jpxy

  write(200,*)delta,Jz2,Jxy,Jpz2,Jpxy,lam,dEz2,dExy


  Nk = (Nkpts-1)*Nkpath

  E0=0d0
  allocate(Hk(5,5,Nk),Eval(5))
  call TB_build_model(Hk,hk_BaNiS2_5orb,5,kpath,Nkpath)
  E0=huge(1d0)
  do ik=1,Nk
     call eigh(Hk(:,:,ik),Eval)
     if(Eval(5)<=E0)E0=Eval(5)
  enddo
  print*,E0
  deallocate(Hk,Eval)  


  Emin = 0d0
  allocate(Hk(Nso,Nso,Nk),Eval(Nso))
  call TB_build_model(Hk,hk_BaNiS2_2bands,Nso,kpath,Nkpath)
  Emin=huge(1d0)
  do ik=1,Nk
     call eigh(Hk(:,:,ik),Eval)
     if(Eval(2)<=Emin)Emin=Eval(2)
  enddo
  print*,Emin
  deallocate(Hk,Eval)  
  !
  E0 = E0-Eshift
  Emin = Emin-Eshift

  !
  call TB_Solve_model(hk_BaNiS2_5orb,5,Kpath,Nkpath,&
       colors_name=[red1,blue1,gray0,gray1,gray2],&
       points_name=[character(len=20) :: '-M', '{/Symbol G}','M'],&
       file="Eigenbands_5b")

  call TB_Solve_model(hk_BaNiS2_2bands,Nso,Kpath,Nkpath,&
       colors_name=[red1,blue1],&
       points_name=[character(len=20) :: '-M', '{/Symbol G}','M'],&
       file="Eigenbands_Dirac")

  Nk = Nkx*Nky
  allocate(Hk(5,5,1),Eval(5))
  allocate(Kgrid(Nk,2))
  call TB_build_kgrid([Nkx,Nky],Kgrid,origin=[-0.5d0,-0.5d0])
  do ikx=1,Nkx
     do iky=1,Nky
        ik = ikx + (iky-1)*Nkx
        Hk(:,:,1) = hk_BaNiS2_5orb(kgrid(ik,:),5)
        call eigh(Hk(:,:,1),Eval)
        write(100,*)kgrid(ik,1),kgrid(ik,2),Eval(4),Eval(5)
     enddo
     write(100,*)""
  enddo

contains


  function hk_BaNiS2_2bands(kpoint,N) result(Mat)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    real(8)                       :: ekz2,ekxy,vkz2xy
    complex(8),dimension(N,N)     :: Mat
    !
    !
    kx=kpoint(1)
    ky=kpoint(2)
    !
    !
    ekz2   = -Jz2*(Cos(kx) + Cos(ky)) - Jpz2*Cos(kx)*Cos(ky)
    ekxy   =  Jxy*(Cos(kx) + Cos(ky)) - Jpxy*Cos(kx)*Cos(ky)
    vkz2xy =  lam*Sin(kx)*Sin(ky)
    !
    Mat(1,1) = dEz2 + ekz2 - Emin
    Mat(2,2) = dExy + ekxy - Emin
    Mat(1,2) = vkz2xy
    Mat(2,1) = vkz2xy
    !
  end function hk_BaNiS2_2bands



  function hk_BaNiS2_5orb(kpoint,N) result(Mat)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: Mat1,Mat2
    complex(8),dimension(N,N) :: Mat
    complex(8)                :: vz22z,vz22x,vz22y,vxy2z,vxy2x,vxy2y
    !
    Mat = zero
    !
    kx=kpoint(1)
    ky=kpoint(2)
    !
    vz22z = -4d0*tz22z*cos(kx/2)*cos(ky/2)
    vz22x =  2d0*tz22a*(sin(kx/2)*cos(ky/2)+sin(ky/2)*cos(kx/2))
    vz22y = -2d0*tz22a*(sin(kx/2)*cos(ky/2)-sin(ky/2)*cos(kx/2))
    !
    vxy2z = -4d0*txy2z*sin(kx/2)*sin(ky/2)
    vxy2x = -2d0*txy2a*(sin(kx/2)*cos(ky/2)+sin(ky/2)*cos(kx/2))
    vxy2y = -2d0*txy2a*(sin(kx/2)*cos(ky/2)-sin(ky/2)*cos(kx/2))
    !
    Mat = diag([Ez2,Exy,Ea,Ea,Ez])-E0*zeye(N)
    Mat(1,3) = vz22x
    Mat(1,4) = vz22y
    Mat(1,5) = vz22z
    Mat(2,3) = vxy2x 
    Mat(2,4) = vxy2y
    Mat(2,5) = vxy2z
    Mat(3,1) = Mat(1,3)
    Mat(3,2) = Mat(2,3)
    Mat(4,1) = Mat(1,4)
    Mat(4,2) = Mat(2,4)
    Mat(5,1) = Mat(1,5)
    Mat(5,2) = Mat(2,5)
    !
    ! call Ensure_Hermiticity(Mat)
    !
  end function hk_BaNiS2_5orb



  subroutine Ensure_Hermiticity(Mat)
    complex(8),dimension(:,:) :: Mat
    integer :: N,i,j
    N = size(Mat,1)
    if(N/=size(Mat,2))stop "Ensure_Hermiticity error: N1!=N2"
    do i=1,N
       do j=i+1,N
          Mat(j,i) = conjg(Mat(i,j))
       enddo
    enddo
  end subroutine Ensure_Hermiticity








  subroutine TB_get_Fermi(Ef)
    real(8)                                 :: Ef
    integer                                 :: ik,Nk,Ns
    real(8)                                 :: mu0,Dmin,Dmax
    !
    Ns = size(Hk,1)
    Nk = size(Hk,3)
    call assert_shape(Hk,[Ns,Ns,Nk],"TB_get_Fermi","Hk")
    allocate(rhoK(Ns,Ns,Nk), eK(Ns,Nk))
    do ik = 1,Nk 
       rhoK(:,:,ik) = Hk(:,:,ik)
       call eigh(rhoK(:,:,ik),eK(:,ik))
    enddo
    Dmin = minval(eK)
    Dmax = maxval(eK)
    mu0 = Dmin
    call fzero(get_dens,mu0,Dmax,info,rguess=Dmin+0.5d0*(Dmax-Dmin))
    if(info/=1)then
       write(*,*)"ERROR TB_get_Fermi: fzero returned info>1 ",info
       stop
    endif
    Ef = -mu0
    write(*,"(A6,12G18.9)")"Ef   =",Ef
    deallocate(rhoK,eK)
  end subroutine TB_get_Fermi

  function get_dens(mu) result(dens)
    real(8),intent(in)                          :: mu
    real(8)                                     :: dens
    real(8),dimension(size(Hk,1))               :: ndens
    real(8),dimension(size(Hk,1),size(Hk,1))    :: Rho
    complex(8),dimension(size(Hk,1),size(Hk,1)) :: Uk_f,diagRho
    ndens = 0d0
    do ik = 1,Nk 
       diagRho = diag(fermi(eK(:,ik)-mu, 1000d0))
       Uk_f    = rhoK(:,:,ik)
       Rho     = (Uk_f .x. diagRho) .x. (conjg(transpose(Uk_f)))
       ndens   = ndens + diagonal(Rho)*Wtk(ik)
    enddo
    write(*,"(A9,3G18.9)")"Ef,N0   =",mu,sum(ndens),filling
    dens = sum(ndens)-filling
  end function get_dens

end program BaCoNiS2
