program BaNiS2
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer                                 :: Norb
  integer                                 :: Nspin,Nso,Nlso
  integer                                 :: Nk,Nkpath,Nkx,Nky,Nkz,Nkpts
  !
  integer                                 :: i,j,k,ik,iorb,jorb,ikpoint
  real(8),dimension(3)                    :: bk1,bk2,ei1,ei2
  real(8),dimension(:,:),allocatable      :: kgrid,kpath,ktrims,Rgrid,Kcrys,bb
  !
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk
  !
  real(8)                                 :: tdz2,tdx2,tdx2x,tdx2y,tdz2x,tdz2y,tdx2x2,tdz2z2
  real(8)                                 :: alat,Efermi
  real(8),dimension(:),allocatable        :: Hdiag
  real(8)                                 :: xmu,beta,eps,wmax
  logical                                 :: logical_folding
  !
  character(len=64)                       :: finput

  call parse_cmd_variable(finput,"finput",default="inputBaNiS2.conf")
  call parse_input_variable(Norb,"norb",finput,default=5)
  call parse_input_variable(nkx,"NKX",finput,default=10)
  call parse_input_variable(nky,"NKY",finput,default=10)
  call parse_input_variable(nkz,"NKZ",finput,default=10)
  call parse_input_variable(nkpath,"nkpath",finput,default=90)
  !
  call parse_input_variable(logical_folding,"logical_folding",finput,default=.false.)
  call parse_input_variable(Efermi,"Efermi",finput,default=8.5d0)
  call parse_input_variable(alat,"alat",finput,default=3.13983813786760132240d0)
  !
  allocate(hdiag(Norb))
  
  call parse_input_variable(hdiag,"hdiag",finput,default=[4.8d0,7.3d0,4.5d0,4.3d0,4.3d0])
  print*,"Hdiag (before):",Hdiag
  Hdiag = Hdiag-Efermi
  print*,"Efermi:",Efermi
  print*,"Hdiag (after):",Hdiag  
  !
  call parse_input_variable(tdx2x2,"tdx2x2",finput,default=0d0)
  call parse_input_variable(tdz2z2,"tdz2z2",finput,default=0d0)
  !
  call parse_input_variable(tdx2,"tdx2",finput,default=0d0)
  call parse_input_variable(tdz2,"tdz2",finput,default=0d0)
  !
  call parse_input_variable(tdx2x,"tdx2x",finput,default=0d0)
  call parse_input_variable(tdx2y,"tdx2y",finput,default=0d0)
  call parse_input_variable(tdz2x,"tdz2x",finput,default=0d0)
  call parse_input_variable(tdz2y,"tdz2y",finput,default=0d0)
  !
  !
  call parse_input_variable(xmu,"xmu",finput,default=0d0)
  call parse_input_variable(beta,"beta",finput,default=1000d0)
  call parse_input_variable(eps,"eps",finput,default=0.1d0)
  call parse_input_variable(wmax,"wmax",finput,default=10d0)
  !
  call save_input_file(reg(finput))

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")


  print*,"WARNING: the bands are not the same as in Michele Casula code, this is due to a difference in the generated path. "
  call sleep(1)

  Nspin = 1 
  Nso   = Nspin*Norb
  Nlso  = 2*Nso

  Nk = Nkx*Nky*Nkz

  ei1 = [ 1d0, 1d0,0d0]
  ei2 = [-1d0, 1d0,0d0]
  call TB_set_ei(ei1,ei2)

  call TB_get_bk(bk1,bk2)

  call TB_set_bk(bk1,bk2)

  call TB_print_bk

  Nkpts = 4
  allocate(Kpath(Nkpts,3))
  allocate(Kcrys(Nkpts,3))
  ! Kpath(1,:) = [0d0,0d0,0d0]
  ! Kpath(2,:) = bk1/2d0
  ! Kpath(3,:) = bk1/2d0 + bk2/2d0
  ! Kpath(4,:) = [0d0,0d0,0d0]
  ! Kpath(1,:) = [0d0*bk1(1)+0d0*bk2(1), 0d0*bk1(2)+0d0*bk2(2)]![0d0,0d0,0d0]
  ! Kpath(2,:) = [0d0*bk1(1)+pi*bk2(1),  0d0*bk1(2)+pi*bk2(2)]  ![0d0, pi,0d0]
  ! Kpath(3,:) = [pi*bk1(1) +pi*bk2(1),  pi*bk1(2)+pi*bk2(2)]    ![pi , pi,0d0]
  ! Kpath(4,:) = [0d0*bk1(1)+0d0*bk2(1), 0d0*bk1(2)+0d0*bk2(2)]![0d0,0d0,0d0]

  ! Kpath(1,:) = [0d0,0d0,0d0] !Gamma
  ! Kpath(2,:) = [0d0, pi,0d0] !X
  ! Kpath(3,:) = [pi , pi,0d0] !M
  ! Kpath(4,:) = [0d0,0d0,0d0] !Gamma


  ! Allocate(bb(3,3))
  ! bb=0d0
  ! bb(:,1) = bk1!*alat
  ! bb(:,2) = bk2!*alat
  ! bb(:,3) = [0d0,0d0,pi2]

  ! do ikpoint=1,4
  !    do j=1,3
  !       Kcrys(ikpoint,j)=0d0
  !       do k=1,3
  !          Kcrys(ikpoint,j)=Kcrys(ikpoint,j)+Kpath(ikpoint,k)*bb(j,k)
  !       enddo
  !    enddo
  !    print*,Kcrys(ikpoint,:)
  ! enddo


  ! Kpath(1,:) = [0d0  ,0d0  ,0d0]
  ! Kpath(2,:) = [-0.25d0,0.25d0,0d0]
  ! Kpath(3,:) = [0d0,0.5d0,0d0]
  ! Kpath(4,:) = [0d0  ,0d0  ,0d0]
  Kpath(1,:) = [0d0,0d0  ,0d0]
  Kpath(2,:) = [1d0,0d0,0d0]
  Kpath(3,:) = [1d0,1d0,0d0]
  Kpath(4,:) = [0d0,0d0  ,0d0]
  Kpath = Kpath*pi


  select case(Norb)
  case (3)
     if(logical_folding)then
        call TB_Solve_model(hk_BaNiS2_3orb_folded,2*Nso,Kpath,Nkpath,&
             colors_name=[red1,blue1,green1,red3,blue3,green3],&
             points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
             file="Eigenbands_3b_folded")
     else
        call TB_Solve_model(hk_BaNiS2_3orb,Nso,Kpath,Nkpath,&
             colors_name=[red1,blue1,green1],&
             points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
             file="Eigenbands_3b")
     endif
     !
  case (5)
     if(logical_folding)then
        call TB_Solve_model(hk_BaNiS2_5orb_Folded,2*Nso,Kpath,Nkpath,&
             colors_name=[red1,blue1,green1,gray0,gray0,blue3,red3,gray3,gray3,green3],&
             points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
             file="Eigenbands_5b_folded")
     else
        call TB_Solve_model(hk_BaNiS2_5orb,Nso,Kpath,Nkpath,&
             colors_name=[red1,blue1,orange1,gray0,gray0],&
             points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
             file="Eigenbands_5b")
     endif
  case default
     stop 
  end select



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
    Mat(1,1) = Hdiag(1) - 2*tdx2x2*(cos(kx)+cos(ky))
    Mat(2,2) = Hdiag(2) - 2*tdz2z2*(cos(kx)+cos(ky))
    Mat(3,3) = Hdiag(3) 
    Mat(1,3)= -4*tdx2*sin(kx/2)*sin(ky/2)
    Mat(2,3)=  4*tdz2*cos(kx/2)*cos(ky/2)
    !
    call Ensure_Hermiticity(Mat)
    !
  end function hk_BaNiS2_3orb


  function hk_BaNiS2_3orb_Folded(kpoint,N) result(Mat)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(Nso,Nso) :: Mat1,Mat2
    complex(8),dimension(N,N)     :: Mat
    !
    if(N/=2*Nso)stop "hk_BaNiS2 error: N!=2*Nso"
    !
    Mat = zero

    kx=kpoint(1)
    ky=kpoint(2)
    !
    Mat1 = zero
    Mat1(1,1) = Hdiag(1) - 2*tdx2x2*(cos(kx)+cos(ky))
    Mat1(2,2) = Hdiag(2) - 2*tdz2z2*(cos(kx)+cos(ky))
    Mat1(3,3) = Hdiag(3) 
    Mat1(1,3)=-4*tdx2*sin(kx/2)*sin(ky/2)
    Mat1(2,3)= 4*tdz2*cos(kx/2)*cos(ky/2)
    Mat(1:3,1:3)   = Mat1
    !
    kpoint = kpoint + bk1
    !
    kx=kpoint(1) 
    ky=kpoint(2)
    Mat1 = zero
    Mat1(1,1) = Hdiag(1) - 2*tdx2x2*(cos(kx)+cos(ky))
    Mat1(2,2) = Hdiag(2) - 2*tdz2z2*(cos(kx)+cos(ky))
    Mat1(3,3) = Hdiag(3) 
    Mat1(1,3)=-4*tdx2*sin(kx/2)*sin(ky/2)
    Mat1(2,3)= 4*tdz2*cos(kx/2)*cos(ky/2)
    Mat(4:6,4:6) = Mat1
    !
    call Ensure_Hermiticity(Mat)
    !
  end function hk_BaNiS2_3orb_Folded




  function hk_BaNiS2_5orb(kpoint,N) result(Mat)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: Mat1,Mat2
    complex(8),dimension(N,N)     :: Mat
    !
    Mat = zero
    !
    kx=kpoint(1)
    ky=kpoint(2)
    !
    Mat = diag(hdiag)
    Mat(1,3)=-4*tdx2*sin(kx/2)*sin(ky/2)
    Mat(1,4)=-2*tdx2x*(sin(kx/2)*cos(ky/2)+sin(ky/2)*cos(kx/2))
    Mat(1,5)= 2*tdx2y*(sin(kx/2)*cos(ky/2)-sin(ky/2)*cos(kx/2))
    Mat(2,3)= 4*tdz2*cos(kx/2)*cos(ky/2)
    Mat(2,4)=-2*tdz2x*(sin(kx/2)*cos(ky/2)+sin(ky/2)*cos(kx/2))
    Mat(2,5)= 2*tdz2y*(sin(kx/2)*cos(ky/2)-sin(ky/2)*cos(kx/2))
    !
    call Ensure_Hermiticity(Mat)
    !
  end function hk_BaNiS2_5orb



  function hk_BaNiS2_5orb_folded(kpoint,N) result(Mat)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(Nso,Nso) :: Mat1,Mat2
    complex(8),dimension(N,N)     :: Mat
    !
    if(N/=2*Nso)stop "hk_BaNiS2 error: N!=2*Nso"
    !
    Mat = zero
    !
    kx=kpoint(1)
    ky=kpoint(2)
    Mat1 = zero
    Mat1 = diag(hdiag)
    ! Mat1(1,1) = Mat1(1,1) - 2*tdx2x2*(cos(kx)+cos(ky))
    ! Mat1(2,2) = Mat1(2,2) - 2*tdz2z2*(cos(kx)+cos(ky))
    Mat1(1,3)=-4*tdx2*sin(kx/2)*sin(ky/2)
    Mat1(2,3)= 4*tdz2*cos(kx/2)*cos(ky/2)
    Mat1(1,4)=-2*tdx2x*(sin(kx/2)*cos(ky/2)+sin(ky/2)*cos(kx/2))
    Mat1(1,5)= 2*tdx2y*(sin(kx/2)*cos(ky/2)-sin(ky/2)*cos(kx/2))
    Mat1(2,4)=-2*tdz2x*(sin(kx/2)*cos(ky/2)+sin(ky/2)*cos(kx/2))
    Mat1(2,5)= 2*tdz2y*(sin(kx/2)*cos(ky/2)-sin(ky/2)*cos(kx/2))
    Mat(:Nso,:Nso)   = Mat1
    !
    kx=kpoint(1) + bk1(1)
    ky=kpoint(2) + bk1(2)
    Mat1 = zero
    Mat1 = diag(hdiag)
    ! Mat1(1,1) = Mat1(1,1) - 2*tdx2x2*(cos(kx)+cos(ky))
    ! Mat1(2,2) = Mat1(2,2) - 2*tdz2z2*(cos(kx)+cos(ky))
    Mat1(1,3)=-4*tdx2*sin(kx/2)*sin(ky/2)
    Mat1(2,3)= 4*tdz2*cos(kx/2)*cos(ky/2)
    Mat1(1,4)=-2*tdx2x*(sin(kx/2)*cos(ky/2)+sin(ky/2)*cos(kx/2))
    Mat1(1,5)= 2*tdx2y*(sin(kx/2)*cos(ky/2)-sin(ky/2)*cos(kx/2))
    Mat1(2,4)=-2*tdz2x*(sin(kx/2)*cos(ky/2)+sin(ky/2)*cos(kx/2))
    Mat1(2,5)= 2*tdz2y*(sin(kx/2)*cos(ky/2)-sin(ky/2)*cos(kx/2))
    Mat(Nso+1:,Nso+1:)   = Mat1
    !
    !Ensure Hermiticity
    call Ensure_Hermiticity(Mat)
    !
  end function hk_BaNiS2_5orb_folded



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
