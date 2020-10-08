program PCO
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=1024,Norb=3,Nspin=1,Nso=Nspin*Norb,Nlso=Nso
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb,io,ispin
  integer                                 :: ix,iy,iz,qr(15)
  real(8)                                 :: kx,ky,kz,bklen
  real(8),dimension(3)                    :: e1,e2,e3,bk1,bk2,bk3,d1,d2,d3
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk

  real(8),dimension(:,:,:),allocatable    :: nkgrid
  real(8)                                 :: ts,xmu,beta,eps,wmax,wmin
  real(8)                                 :: n(Nso)
  complex(8)                              :: w,Hloc(Nso,Nso)
  complex(8)                              :: Gmats(Nspin,Nspin,Norb,Norb,L)
  complex(8)                              :: Greal(Nspin,Nspin,Norb,Norb,L)
  character(len=20)                       :: nkstring
  logical                                 :: iexist,ibool,iener


  call parse_input_variable(Nkx,"NKX","inputHM.conf",default=10)
  call parse_input_variable(nkpath,"NKPATH","inputHM.conf",default=500)
  call parse_input_variable(xmu,"XMU","inputHM.conf",default=11.037956d0)
  call parse_input_variable(eps,"EPS","inputHM.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM.conf",default=1000.d0)
  call parse_input_variable(wmin,"WMIN","inputHM.conf",default=-10d0)
  call parse_input_variable(wmax,"WMAX","inputHM.conf",default=10d0)
  call save_input_file("inputHM.conf")

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")

  
  !>Direct lattice basis vector
  call TB_set_ei([1d0,0d0,0d0],[0d0,1d0,0d0],[0d0,0d0,1d0])
  call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])



  !METHOD 1 (setup W90 --> use internal W90 model)
  call TB_w90_setup("W90_hr_bulk.dat",nlat=1,nspin=Nspin,norb=Norb)

  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  call TB_build_model(Hk,TB_w90_model,Nso,[Nkx,Nkx,Nkx])
  Wtk = 1d0/Nktot
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"w90Hloc.dat")


  call get_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call get_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  call print_gf_matsubara(Gmats,"w90Gloc",iprint=1)
  call print_gf_realaxis(Greal,"w90Gloc",iprint=1)
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb + (ispin-1)*Norb
        n(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
     enddo
  enddo

  !solve along the standard path in the 3D BZ.

  Npts = 8
  Nk=(Npts-1)*Nkpath
  allocate(kpath(Npts,3))
  kpath(1,:)=[0,0,0]*pi
  kpath(2,:)=[1,0,0]*pi
  kpath(3,:)=[1,1,0]*pi
  kpath(4,:)=[0,0,0]*pi
  kpath(5,:)=[1,1,1]*pi
  kpath(6,:)=[1,0,1]*pi
  kpath(7,:)=[0,0,1]*pi
  kpath(8,:)=[0,0,0]*pi
  call TB_Solve_model(TB_w90_model,Nso,kpath,Nkpath,&
       colors_name=[red1,green1,blue1],&
       points_name=[character(len=20) ::'G', 'X', 'M', 'G', 'R', 'A', 'Z','G'],&
       file="w90Eigenband.nint")


  !plot observables
  open(10,file="w90observables.nint")
  write(10,"(10F20.12)")(n(iorb),iorb=1,Nso),sum(n)
  close(10)
  write(*,"(A,10F14.9)")"w90 Occupations =",(n(iorb),iorb=1,Nso),sum(n)


  call TB_w90_delete()


end program PCO


