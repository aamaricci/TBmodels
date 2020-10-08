! DESCRIPTION
!   Solve the non-interacting BHZ model and generate the hamiltoniana matrix H(k), 
!   Using the 2d square lattice dispersion. 
!
! MODEL Hamiltonian is:
!
! |     h^{2x2}(k)              &         0                   |
! |      0                      &        [h^{2x2}]*(-k)       |
!
!
! h^{2x2}(k):=
!
! | m-(Cos{kx}+Cos{ky})         & \lambda*(Sin{kx}-i*Sin{ky}) |
! | \lambda*(Sin{kx}+i*Sin{ky}) & -m+(Cos{kx}+Cos{ky})        |
!
program bhz_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                           :: Norb=2,Nspin=2,Nso=Nspin*Norb
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

  real(8)                                     :: chern,z2
  real(8)                                     :: mh,lambda
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

  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(L,"L","inputBHZ.conf",default=2048)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=1d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz.in")
  
  call save_input_file("inputBHZ.conf")
  
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


  !SETUP THE GAMMA MATRICES: 4x4
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)

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
       colors_name=[red1,blue1,red1,blue1],&
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
        dens(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
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
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
  end function hk_model

end program


