! DESCRIPTION
!   Solve the non-interacting BHZ model and generate the hamiltoniana matrix H(k), 
!   Using the 2d square lattice dispersion. 
!
! MODEL Hamiltonian is:
!
! | m-2*t*(Cos{kx}+Cos{ky})         & \lambda*(Cos{kx}-Cos{ky}) |
! | \lambda*(Cos{kx}-Cos{ky}) & -m-2*t*(Cos{kx}+Cos{ky})        |
!
program hm_2b_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                           :: Norb=2,Nspin=1,Nso=Nspin*Norb
  !
  integer                                     :: Nk,Nktot,Nkpath,Nkx,Npts,L
  integer                                     :: Nky,Nlat,Nx,Ny
  integer                                     :: i,j,k,ik,iorb,jorb,ispin,io
  integer :: ilat,jlat
  integer                                     :: ix,iy,iz
  real(8)                                     :: kx,ky,kz
  real(8),dimension(:,:),allocatable          :: kgrid,kpath,ktrims,Rgrid
  integer,dimension(:,:),allocatable          :: Links
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:,:,:),allocatable   :: Hlat
  real(8),dimension(:),allocatable            :: Wtk

  real(8)                                     :: mh,lambda,ts,tsgn
  real(8)                                     :: xmu,beta,eps,Eout(2)
  real(8)                                     :: dens(Nso)
  complex(8)                                  :: Hloc(Nso,Nso),arg
  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats,Greal
  character(len=20)                           :: file
  logical                                     :: iexist
  complex(8),dimension(Nso,Nso)               :: Gamma1,Gamma2,Gamma5
  complex(8),dimension(:,:,:),allocatable     :: ftHk
  complex(8),dimension(:,:,:,:),allocatable   :: ftHlat
  real(8),dimension(2)                        :: vecK,vecRi,vecRj

  call parse_input_variable(nkx,"NKX","inputHM2b.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputHM2b.conf",default=500)
  call parse_input_variable(L,"L","inputHM2b.conf",default=2048)
  call parse_input_variable(ts,"TS","inputHM2b.conf",default=0.5d0)
  call parse_input_variable(tsgn,"TSGN","inputHM2b.conf",default=1d0)
  call parse_input_variable(mh,"MH","inputHM2b.conf",default=1d0)
  call parse_input_variable(lambda,"LAMBDA","inputHM2b.conf",default=0.3d0)
  call parse_input_variable(xmu,"XMU","inputHM2b.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputHM2b.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM2b.conf",default=1000.d0)
  call parse_input_variable(file,"FILE","inputHM2b.conf",default="hkfile_hm.in")

  call save_input_file("inputHM2b.conf")

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")


  Nky  = Nkx
  Nktot= Nkx*Nky
  !
  Nx   = Nkx
  Ny   = Nkx
  Nlat = Nx*Ny

  !Define unit basis in real space
  call TB_set_ei([1d0,0d0],[0d0,1d0])

  !Define unit basis in reciprocal space pi2 is in Scifor
  call TB_set_bk([pi2,0d0],[0d0,pi2])


  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:  
  write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])
  Wtk = 1d0/Nktot

  call TB_write_hk(Hk,trim(file),Nso,&
       Nd=Norb,Np=1,Nineq=1,&
       Nkvec=[Nkx,Nkx])


  !GET LOCAL PART OF THE HAMILTONIAN
  Hloc=sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc)


  !SOLVE ALONG A PATH IN THE BZ.
  Npts=5
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_X1
  kpath(2,:)=kpoint_Gamma
  kpath(3,:)=kpoint_M1
  kpath(4,:)=kpoint_X1
  kpath(5,:)=kpoint_Gamma
  call TB_Solve_model(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[red1,blue1],&
       points_name=[character(len=20) :: 'X', 'G', 'M', 'X', 'G'],&
       file="Eigenband.nint")




  !Build the local GF:
  allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
  allocate(Greal(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_matsubara(Gmats,"Gloc",1)
  call dmft_print_gf_realaxis(Greal,"Gloc",1)


  !Occupation:
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb+(ispin-1)*Norb
        dens(io) = 2*fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)   !  *2 spin degeneracy
     enddo
  enddo
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(dens(iorb),iorb=1,Nso),sum(dens)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(dens(iorb),iorb=1,Nso),sum(dens)

  !Kinetic Energy
  call dmft_kinetic_energy(Hk,Wtk,zeros(Nspin,Nspin,Norb,Norb,L))



contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8) :: ek
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    if(N/=2)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -2d0*ts*(cos(kx)+cos(ky))
    Hk(1,1) = Mh + ek
    Hk(1,2) = lambda*(cos(kx)-cos(ky))
    Hk(2,1) = lambda*(cos(kx)-cos(ky))
    Hk(2,2) = -Mh + tsgn/abs(tsgn)*ek
  end function hk_model

end program


