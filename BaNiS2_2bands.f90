program BaNiS2
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                 :: Norb
  integer                                 :: Nspin,Nso,Nlso
  integer                                 :: Nk,Nkpath,Nkx,Nky,Nkz,Nkpts
  !
  integer                                 :: i,j,k,ik,iorb,jorb,ikpoint,Lreal,ispin
  real(8),dimension(2)                    :: bk1,bk2,ei1,ei2
  real(8),dimension(:,:),allocatable      :: kpath,ktrims
  !
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk,wreal,kgrid
  !
  real(8)                                 :: Ez2,Exy,Ez,Ea
  real(8)                                 :: tz22z,tz22a,txy2z,txy2a
  real(8)                                 :: t1,t2,t3,t4
  real(8)                                 :: bEz2,bExy
  real(8)                                 :: Efermi
  real(8)                                 :: xmu,beta,eps,wmax,wmin
  logical                                 :: logical_folding
  complex(8),dimension(:,:,:,:,:,:),allocatable :: Gkreal
  real(8),dimension(:,:),allocatable            :: Akreal
  !
  character(len=64)                       :: finput
  !
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
  !
  Efermi = -8.82;
  Ez2 = 7.48 + Efermi; 
  Exy = 7.54 + Efermi;
  Ez  = 6.85 + Efermi;
  Ea  = 6.33 + Efermi;
  tz22z = 0.58;
  tz22a = 0.11;
  txy2z = 0.37;
  txy2a = 0.88
  !
  t1 = 4*tz22z**2/Ez
  t2 = 4*txy2z**2/Ez
  t3 = 4*tz22a**2/Ea
  t4 = 4*txy2a**2/Ea
  !
  bEz2 = Ez2 - t1 - t3;
  bExy = Exy - t2 - t4;
  !
  print*,t1,t2
  print*,(t1-t3),(t2-t4)
  print*,-Sqrt(t1*t2) + Sqrt(t3*t4)
  print*,bEz2,bExy
  stop
  Lreal = 300
  !
  call TB_set_bk([pi2,0d0],[0d0,pi2])
  call TB_print_bk
  !
  allocate(Hk(Nso,Nso,Nk))
  ! call TB_set_dos_range([wmin,wmax])
  ! call TB_set_dos_lreal(5000)
  ! call TB_set_dos_eps(eps)
  call TB_build_model(Hk,hk_BaNiS2_2bands,Nso,[Nkx,Nky],wdos=.false.)
  !
  !
  Nkpts = 2
  allocate(Kpath(Nkpts,2))
  Kpath(1,:) = [0d0,0d0] !Gamma
  Kpath(2,:) = [pi, 0d0] !X
  !
  !
  call TB_Solve_model(hk_BaNiS2_2bands,Nso,Kpath,Nkpath,&
       colors_name=[red1,blue1,green1],&
       points_name=[character(len=20) :: '{\Symbol G}', 'M'],&
       file="Eigenbands_Dirac")
  deallocate(Hk)


  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(0.16d0,"wfin")
  call add_ctrl_var(eps,"eps")
  call print_ctrl_list


  Nk = (Nkpts-1)*Nkpath
  allocate(Hk(Nso,Nso,Nk))
  call TB_build_model(Hk,hk_BaNiS2_2bands,Nso,kpath,Nkpath)

  allocate(Gkreal(Nk,Nspin,Nspin,Norb,Norb,Lreal))
  call start_timer
  do ik=1,Nk
     call dmft_gk_realaxis(Hk(:,:,ik),1d0,Gkreal(ik,:,:,:,:,:),zeros(Nspin,Nspin,Norb,Norb,Lreal))
     call eta(ik,Nk)
  enddo
  call stop_timer
  !
  allocate(Akreal(Nk,Lreal))
  Akreal = zero
  do ispin=1,Nspin
     do iorb=1,Norb
        Akreal = Akreal -dimag(Gkreal(:,ispin,ispin,iorb,iorb,:))/pi/Nspin/Norb
     enddo
  enddo
  !
  allocate(wreal(Lreal))
  allocate(kgrid(Nk))
  wreal = linspace(wmin,0d0,Lreal)
  kgrid = linspace(0d0,pi,Nk,iend=.false.)
  call splot3d("Akw_realw.dat",kgrid,wreal,Akreal(:,:))



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
    ekz2   = -t1*(Cos(kx) + Cos(ky)) - (t1 - t3)*Cos(kx)*Cos(ky)
    ekxy   =  t2*(Cos(kx) + Cos(ky)) - (t2 - t4)*Cos(kx)*Cos(ky)
    vkz2xy = -(Sqrt(t1*t2) - Sqrt(t3*t4))*Sin(kx)*Sin(ky)
    !
    Mat(1,1) = bEz2 + ekz2
    Mat(2,2) = bExy + ekxy
    Mat(1,2) = vkz2xy
    Mat(2,1) = vkz2xy
    !
  end function hk_BaNiS2_2bands


end program BaNiS2
