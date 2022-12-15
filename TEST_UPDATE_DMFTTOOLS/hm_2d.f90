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
  complex(8)                              :: Gmats(Nspin,Nspin,Norb,Norb,L),Giw(Nso,Nso,L)
  complex(8)                              :: Greal(Nspin,Nspin,Norb,Norb,L),Gwr(Nso,Nso,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,iener


  call parse_input_variable(nkx,"NKX","inputHM.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputHM.conf",default=500)
  call parse_input_variable(ts,"TS","inputHM.conf",default=0.5d0)
  call parse_input_variable(xmu,"XMU","inputHM.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputHM.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM.conf",default=1000.d0)
  call parse_input_variable(wmax,"WMAX","inputHM.conf",default=10d0)
  call save_input_file("inputHM.conf")

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")



  !>Direct lattice basis vector
  call TB_set_bk([pi2,0d0],[0d0,pi2])




  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"Hloc.dat")

  call dmft_gloc_realaxis(Hk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_get_gloc(Hk,Gwr,zeros(Nso,Nso,L),axis="real")
  call dmft_write_gf(Greal,"Greal",axis="real")
  call dmft_write_gf(Gwr,"Gwr",axis="real")


  call dmft_gloc_matsubara(Hk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_get_gloc(Hk,Giw,zeros(Nso,Nso,L),"mats")
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb + (ispin-1)*Norb
        n(io) = fft_get_density(Giw(io,io,:),beta);print*,n(io)
        n(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta);print*,n(io)
     enddo
  enddo
  call dmft_write_gf(Gmats,"Gmats",axis="mats")
  call dmft_write_gf(Giw,"Giw",axis="mats")



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




end program hm_2d


