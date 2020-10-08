program kanemele
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                           :: Norb=2,Nspin=2,Nlat=1,Nlso=Nlat*Nspin*Norb
  integer                                     :: Nkx,Nky,Nk,Nkpath,Npts,L

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)     :: Hk     !the Hamiltonian H(k)
  complex(8),allocatable,dimension(:,:)       :: pxpyHloc !loc Hamiltonian [Nlso][Nlso]
  complex(8),allocatable,dimension(:,:,:,:)   :: Hloc   !loc Hamiltonian [Nspin][Nspin][Norb][Norb]
  real(8),allocatable,dimension(:)            :: Wtk    !weight of the k-points

  !parameters for the model:
  real(8)                                     :: Vsigma,Vpi,LamISB,LamSOC
  real(8)                                     :: xmu,beta,eps
  !
  character(len=32)                           :: finput !input file name
  !
  !Aux. variables
  integer                                     :: io,ispin,iorb
  real(8),dimension(:,:),allocatable          :: KPath !high symmetry path in the BZ [Npts,Dim]
  complex(8),dimension(:,:,:,:,:),allocatable :: Greal,Gmats !local GFs [Nspin,Nspin,Norb,Norb,L]
  real(8),dimension(Nlso)                     :: dens        !occupations
  complex(8),dimension(4,4)                   :: Gamma0,GammaX,GammaY,GammaZ,GammaS

  !Parse variables
  call parse_cmd_variable(finput,"FINPUT",default='inputPXPY.conf')
  !
  call parse_input_variable(Nkx,"NKX",finput,default=10,comment="# of k points along X")
  call parse_input_variable(Nky,"NKY",finput,default=10,comment="# of k points along Y")
  call parse_input_variable(nkpath,"NKPATH",finput,default=100,comment="# of points along BZ path segment")
  call parse_input_variable(L,"L",finput,default=2048,comment="# of frequencies")
  !
  call parse_input_variable(Vsigma,"Vsigma",finput,default=1d0)
  call parse_input_variable(Vpi,"Vpi",finput,default=-1d0)
  call parse_input_variable(LamISB,"LAMISB",finput,default=0.1d0)
  call parse_input_variable(LamSOC,"LAMSOC",finput,default=0d0)
  !aux. variables
  call parse_input_variable(xmu,"XMU",finput,default=0d0)
  call parse_input_variable(eps,"EPS",finput,default=4d-2)
  call parse_input_variable(beta,"BETA",finput,default=1000d0)
  !
  call save_input_file(finput)  !print out and write used.input.file
  call print_input()

  !PUSH viariables to DMFT_TOOLS memory pool:
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps") !w + i*eps

  Gamma0 = kron_pauli( pauli_tau_0, pauli_sigma_0)
  GammaX = kron_pauli( pauli_tau_x, pauli_sigma_0)
  GammaY = kron_pauli( pauli_tau_y, pauli_sigma_0)
  GammaZ = kron_pauli( pauli_tau_z, pauli_sigma_0)
  GammaS = kron_pauli( pauli_tau_y, pauli_sigma_z)

  !Now the main work starts:

  !STEP 1: set up the direct- and reciprocal-lattice vectors basis (a=1)
  !> setup real space basis
  call TB_set_ei([1d0,0d0],[-0.5d0,sqrt(3d0)/2d0])
  !> get reciprocal lattice basis
  call TB_build_bk(verbose=.true.)

  !STEP 2: build H(k)
  Nk = Nkx*Nky
  allocate(Hk(Nlso,Nlso,Nk))
  allocate(Wtk(Nk))
  !> fill up H(k) matrix by iteratively calling user provided function:
  call TB_build_model(Hk,hk_triang_pxpy,Nlso,[Nkx,Nky],wdos=.false.)
  Wtk = 1d0/Nk

  Npts = 5
  allocate(Kpath(Npts,2))
  KPath(1,:)=[0d0,0d0]
  KPath(2,:)=[0.5d0,0d0]
  Kpath(3,:)=[1d0/3d0,1d0/3d0]
  KPath(4,:)=[-1d0/3d0,2d0/3d0]
  Kpath(5,:)=[0d0,0d0]
  Kpath=Kpath*pi2
  call TB_Solve_model(hk_triang_pxpy,Nlso,KPath,Nkpath,&
       colors_name=[red,blue,orange,green],&
       points_name=[character(len=10) :: "{/Symbol G}","M","K","K`","{/Symbol G}"],&
       file="PxPybands.nint",iproject=.false.)


  !Build the local GF:
  allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
  allocate(Greal(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)


  !Occupation:
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb + (ispin-1)*Norb
        dens(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
     enddo
  enddo
  write(*,"(A,20F14.9)")"Occupations =",(dens(io),io=1,Nlso),sum(dens)


contains

  !Fortran FUNCTION with two arguments mandatory: Kvec (dble, dim=2), N (integer) 
  function hk_triang_pxpy(kvec,N) result(Hk)
    !input
    real(8),dimension(:)      :: kvec
    integer                   :: N
    !output
    complex(8),dimension(N,N) :: Hk !the matrix H(k) for a given k=kvec
    !aux.
    real(8) :: kx,ky
    real(8) :: cx,cy,cxy
    real(8) :: sx,sy,sxy
    !
    if(N/=4) stop "hk_triang_pxpy error: N != 4"
    !
    kx  = kvec(1)
    ky  = kvec(2)
    cx  = cos(kx)
    cy  = cos(ky)
    cxy = cos(kx+ky)
    sx  = sin(kx)
    sy  = sin(ky)
    sxy = sin(kx+ky)
    !
    Hk = (Vsigma+Vpi)*(cx + cy + cxy)*Gamma0 + &
         sqrt(3d0)/2d0*(Vsigma-Vpi)*(cxy-cy)*GammaX + &
         LamISB*(sx + sy - sxy)*GammaY + &
         0.5d0*(Vsigma-Vpi)*(2*cx - cy - cxy)*GammaZ + &
         LamSOC*GammaS
    !
  end function hk_triang_pxpy

end program kanemele
