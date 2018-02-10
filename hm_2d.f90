program hm_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=1024,Norb=1,Nspin=1,Nso=Nspin*Norb
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb,io,ispin
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz,bklen
  real(8),dimension(3) :: e1,e2,e3,bk1,bk2,bk3,d1,d2,d3
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk

  real(8),dimension(:,:,:),allocatable    :: nkgrid
  real(8)                                 :: ts,xmu,beta,eps,wmax
  real(8)                                 :: n(Nso)
  complex(8)                              :: w,Hloc(Nso,Nso)
  complex(8)                              :: Gmats(Nspin,Nspin,Norb,Norb,L)
  complex(8)                              :: Greal(Nspin,Nspin,Norb,Norb,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,iener


  call parse_input_variable(nkx,"NKX","inputHM.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputHM.conf",default=500)
  call parse_input_variable(ts,"TS","inputHM.conf",default=0.5d0)
  call parse_input_variable(xmu,"XMU","inputHM.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputHM.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM.conf",default=1000.d0)
  call parse_input_variable(wmax,"WMAX","inputHM.conf",default=10d0)
  call parse_input_variable(file,"FILE","inputHM.conf",default="hkfile_bhz.in")
  call save_input_file("inputHM.conf")

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")



  !>Direct lattice basis vector
  call TB_set_bk([1d0,0d0],[0d0,1d0])




  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])
  Wtk = 1d0/Nktot
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"Hloc.dat")


  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb + (ispin-1)*Norb
        n(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
     enddo
  enddo

  !solve along the standard path in the 2D BZ.
  Npts=4
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_Gamma
  kpath(2,:)=kpoint_M1
  kpath(3,:)=kpoint_X1
  kpath(4,:)=kpoint_Gamma
  call TB_Solve_model(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[red1],&
       points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
       file="Eigenband.nint")


  !plot observables
  open(10,file="observables.nint")
  write(10,"(10F20.12)")(n(iorb),iorb=1,Nso),sum(n)
  close(10)
  write(*,"(A,10F14.9)")"Occupations =",(n(iorb),iorb=1,Nso),sum(n)

contains




  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = -2d0*ts*(cos(kx)+cos(ky))
  end function hk_model



  subroutine reciprocal_basis(e1,e2,e3, bk1,bk2,bk3)
    !   This routine generates the reciprocal lattice vectors b1,b2,b3
    !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
    real(8),dimension(:),intent(in)                 :: e1
    real(8),dimension(size(e1)),intent(in),optional :: e2
    real(8),dimension(size(e1)),intent(in),optional :: e3
    real(8),dimension(size(e1))                     :: bk1
    real(8),dimension(size(e1)),optional            :: bk2
    real(8),dimension(size(e1)),optional            :: bk3
    !
    real(8),dimension(3)                            :: a1
    real(8),dimension(3)                            :: a2
    real(8),dimension(3)                            :: a3
    real(8),dimension(3)                            :: b1
    real(8),dimension(3)                            :: b2
    real(8),dimension(3)                            :: b3
    real(8)                                         :: den, s
    integer                                         :: iperm, i, j, k, l, ipol
    integer                                         :: N
    real(8),dimension(3,3)                          :: Mat
    real(8),parameter                               :: pi2 = 6.28318530717958623200
    N = size(e1)

    a1=[1d0,0d0,0d0]
    a2=[0d0,1d0,0d0]
    a3=[0d0,0d0,1d0]

    a1(:N)=e1
    if(present(e2))a2(:N)=e2
    if(present(e3))a3(:N)=e3

    den = 0
    i = 1
    j = 2
    k = 3
    s = 1.d0
100 do iperm = 1, 3
       den = den + s * a1 (i) * a2 (j) * a3 (k)
       l = i
       i = j
       j = k
       k = l
    enddo
    i = 2
    j = 1
    k = 3
    s = - s
    if (s.lt.0.d0) goto 100

    Mat(:,1) = a1
    Mat(:,2) = a2
    Mat(:,3) = a3
    print*,"Det(M),den:",det(Mat),den

    !
    !    here we compute the reciprocal vectors
    !
    i = 1
    j = 2
    k = 3
    do ipol = 1, 3
       b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) )/den*pi2
       b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) )/den*pi2
       b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) )/den*pi2
       l = i
       i = j
       j = k
       k = l
    enddo

    bk1(:N)=b1
    if(present(bk2))bk2(:N)=b2
    if(present(bk3))bk3(:N)=b3

    return
  end subroutine reciprocal_basis




end program hm_2d


