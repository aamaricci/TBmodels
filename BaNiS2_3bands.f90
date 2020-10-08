program BaNiS2
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                 :: Norb
  integer                                 :: Nspin,Nso,Nlso
  integer                                 :: Nk,Nkpath,Nkx,Nky,Nkz,Nkpts
  !
  integer                                 :: i,j,k,ik,iorb,jorb,ikpoint
  real(8),dimension(2)                    :: bk1,bk2,ei1,ei2
  real(8),dimension(:,:),allocatable      :: kgrid,kpath,ktrims,Rgrid,Kcrys,bb
  !
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk
  !
  real(8)                                 :: t,g
  real(8)                                 :: alat,Efermi
  real(8),dimension(:),allocatable        :: Hdiag
  real(8)                                 :: xmu,beta,eps,wmax,wmin
  logical                                 :: logical_folding
  !
  character(len=64)                       :: finput

  Norb = 3

  call parse_cmd_variable(finput,"finput",default="inputBaNiS2.conf")
  call parse_input_variable(nkx,"NKX",finput,default=30)
  call parse_input_variable(nky,"NKY",finput,default=30)
  call parse_input_variable(nkpath,"nkpath",finput,default=100)
  !
  call parse_input_variable(logical_folding,"logical_folding",finput,default=.false.)
  !
  allocate(hdiag(Norb))
  call parse_input_variable(hdiag,"hdiag",finput,default=[3d0,0d0,-4d0])
  !
  call parse_input_variable(t,"t",finput,default=1d0)
  call parse_input_variable(g,"g",finput,default=1d0)
  !
  call parse_input_variable(xmu,"xmu",finput,default=0d0)
  call parse_input_variable(beta,"beta",finput,default=1000d0)
  call parse_input_variable(eps,"eps",finput,default=0.1d0)
  call parse_input_variable(wmax,"wmax",finput,default=10d0)
  call parse_input_variable(wmin,"wmin",finput,default=-10d0)
  !
  call save_input_file(reg(finput))

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")


  Nspin = 1 
  Nso   = Nspin*Norb
  Nlso  = 2*Nso

  Nk = Nkx*Nky

  call TB_set_bk([pi2,0d0],[0d0,pi2])
  call TB_print_bk

  allocate(Hk(Nso,Nso,Nk))

  call TB_set_dos_range([wmin,wmax])
  call TB_set_dos_lreal(5000)
  call TB_set_dos_eps(eps)
  call TB_build_model(Hk,hk_BaNiS2_3orb,Nso,[Nkx,Nky])


  Nkpts = 2
  allocate(Kpath(Nkpts,2))
  Kpath(1,:) = [0d0,0d0] !Gamma
  Kpath(2,:) = [pi, 0d0] !X


  call TB_Solve_model(hk_BaNiS2_3orb,Nso,Kpath,Nkpath,&
       colors_name=[red1,blue1,green1],&
       points_name=[character(len=20) :: 'G', 'X'],&
       file="Eigenbands_3b")




contains




  function hk_BaNiS2_3orb(kpoint,N) result(Mat)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N) :: Mat1,Mat2
    complex(8),dimension(N,N)     :: Mat
    !
    !
    kx=kpoint(1)
    ky=kpoint(2)
    !
    Mat = zero
    Mat(1,1) = Hdiag(1) - 2*t*(cos(kx)+cos(ky))
    Mat(2,2) = Hdiag(2)
    Mat(3,3) = Hdiag(3) 
    Mat(1,3)= -4*g*sin(kx/2)*sin(ky/2)
    Mat(2,3)=  4*g*cos(kx/2)*cos(ky/2)
    !
    call Ensure_Hermiticity(Mat)
    !
  end function hk_BaNiS2_3orb



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

end program BaNiS2
