! DESCRIPTION
!   Solve the non-interacting weyl model and generate the hamiltoniana matrix H(k), 
!   Model parameters are: t_dd,v0/tpd,ep0,ed0.
!   Using the 3d square lattice dispersion. 
!   The matrix is written in the Wannier90 form, as expected by w2CTQMC & ggED code.
!   Output is on file *hkfile.in default.
!
program weyl_3d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                           :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                     :: Nk,Nktot,Nkpath,Nkx,Npts,L,z2(4)
  integer                                     :: i,j,k,ik,iorb,jorb,io,ispin
  integer                                     :: ix,iy,iz
  real(8)                                     :: kx,ky,kz
  real(8),dimension(:,:),allocatable          :: kgrid,kpath,ktrims
  complex(8),dimension(:,:,:),allocatable     :: Hk
  real(8),dimension(:),allocatable            :: Wtk

  real(8)                                     :: chern
  real(8)                                     :: mh,rh,lambda,delta,bx,by,bz,BIA
  real(8)                                     :: xmu,beta,eps,e0
  real(8)                                     :: dens(Nso)
  complex(8)                                  :: Hloc(Nso,Nso)
  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats,Greal,Sfoo !(Nspin,Nspin,Norb,Norb,L)
  complex(8),dimension(:,:,:),allocatable     :: Skin !(Nso,Nso,L)
  character(len=20)                           :: file
  logical                                     :: iexist, gfflag
  complex(8),dimension(Nso,Nso)               :: Gamma1,Gamma2,Gamma3,Gamma5

  call parse_input_variable(gfflag,"GFFLAG","inputweyl.conf",default=.true.)
  call parse_input_variable(nkx,"NKX","inputweyl.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputweyl.conf",default=500)
  call parse_input_variable(L,"L","inputweyl.conf",default=2048)
  call parse_input_variable(mh,"MH","inputweyl.conf",default=1d0)
  call parse_input_variable(e0,"E0","inputweyl.conf",default=1d0)
  call parse_input_variable(rh,"RH","inputweyl.conf",default=0d0)
  call parse_input_variable(lambda,"LAMBDA","inputweyl.conf",default=0.5d0)
  call parse_input_variable(delta,"DELTA","inputweyl.conf",default=0d0)
  call parse_input_variable(xmu,"XMU","inputweyl.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputweyl.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputweyl.conf",default=1000.d0)
  call parse_input_variable(file,"FILE","inputweyl.conf",default="hkfile_weyl.in")
  call parse_input_variable(bx,"BX","inputweyl.conf",default=0.3.d0)
  call parse_input_variable(by,"BY","inputweyl.conf",default=0.d0)
  call parse_input_variable(bz,"BZ","inputweyl.conf",default=0.d0)
  call parse_input_variable(BIA,"BIA","inputweyl.conf",default=0.d0)
  call save_input_file("inputweyl.conf")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")



  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma3=kron_pauli( pauli_tau_x,-pauli_sigma_x)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)


  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx*Nkx
  write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2]) 
  call TB_build_model(Hk,hk_weyl,Nso,[Nkx,Nkx,Nkx])
  Wtk = 1d0/Nktot

  call TB_write_hk(Hk,trim(file),Nso,&
       Nd=Norb,Np=1,Nineq=1,&
       Nkvec=[Nkx,Nkx,Nkx])


  !GET LOCAL PART OF THE HAMILTONIAN
  Hloc=sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc)

  !solve along the standard path in the 3D BZ.
  Npts=8
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_m1
  kpath(2,:)=kpoint_x2
  kpath(3,:)=kpoint_gamma
  kpath(4,:)=kpoint_x1
  kpath(5,:)=kpoint_m2
  kpath(6,:)=kpoint_r
  kpath(7,:)=kpoint_x3
  kpath(8,:)=kpoint_gamma
  call TB_solve_model(hk_weyl,Nso,kpath,Nkpath,&
       colors_name=[red1,blue1,red1,blue1],&
       points_name=[character(len=10) :: "M","X","G","X1","A","R","Z","G"],&
       file="Eigenband.nint")

  if (gfflag) then
     !Build the local GF:
     allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
     allocate(Greal(Nspin,Nspin,Norb,Norb,L))
     allocate(Sfoo(Nspin,Nspin,Norb,Norb,L))
     Gmats=zero
     Greal=zero
     Sfoo =zero
     call dmft_gloc_matsubara(Hk,Wtk,Gmats,Sfoo,iprint=3)
     Sfoo =zero
     call dmft_gloc_realaxis(Hk,Wtk,Greal,Sfoo,iprint=3)


     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb
           dens(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
        enddo
     enddo
     !plot observables
     open(10,file="observables.nint")
     write(10,"(20F20.12)")(dens(io),io=1,Nso),sum(dens)
     close(10)
     write(*,"(A,20F14.9)")"Occupations =",(dens(iorb),iorb=1,Nso),sum(dens)

     Sfoo = zero
     call dmft_kinetic_energy(Hk,Wtk,Sfoo)
  endif

  !Evaluate the Z2 index:
  !STRONG TI
  z2(1) = z2_number(reshape( [ [0,0,0] , [1,0,0] , [1,1,0] , [0,1,0] , [0,1,1] , [0,0,1] , [1,0,1] , [1,1,1] ] , [3,8] )*pi)
  !WEAK TI
  !K=1: n_1=1, n_2,3=0,1
  z2(2) = z2_number(reshape( [ [1,0,0] , [1,1,0] , [1,1,1] , [1,0,1] ] , [3,4])*pi)
  !K=2: n_2=1, n_1,2=0,1
  z2(3) = z2_number(reshape( [ [0,1,0] , [0,1,1] , [1,1,1] , [1,1,0] ] , [3,4])*pi)
  !k=3: n_3=1, n_1,2=0,1
  z2(4) = z2_number(reshape( [ [0,0,1] , [0,1,1] , [1,1,1] , [1,0,1] ] , [3,4])*pi)
  open(100,file="z2_invariant.nint")
  write(100,*)z2
  close(100)


contains




  function hk_weyl(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: kx,ky,kz
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_weyl: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    kz=kpoint(3)
    !
    Hk          = zero
    Hk(1:2,1:2) = &
         (Mh - e0*(cos(kx) + cos(ky) + cos(kz)) )*pauli_tau_z +&
         lambda*sin(kx)*pauli_tau_x + lambda*sin(ky)*pauli_tau_y +&
         by*pauli_tau_y + bz*pauli_tau_z
    Hk(3:4,3:4) = conjg( &
         (Mh - e0*(cos(-kx) + cos(-ky) + cos(-kz)) )*pauli_tau_z +&
         lambda*sin(-kx)*pauli_tau_x + lambda*sin(-ky)*pauli_tau_y +&
         by*pauli_tau_y - bz*pauli_tau_z)
    Hk(1:2,3:4) = lambda*sin(kz)*pauli_tau_x - BIA*pauli_tau_x + bx*pauli_tau_z
    Hk(3:4,1:2) = lambda*sin(kz)*pauli_tau_x + BIA*pauli_tau_x + bx*pauli_tau_z
    !
    Hk = Hk
    !
  end function hk_weyl

  function z2_number(ktrims) result(z2)
    real(8),dimension(:,:),intent(in)       :: ktrims
    complex(8),dimension(:,:,:),allocatable :: Htrims
    complex(8),dimension(:),allocatable     :: Delta
    integer                                 :: z2
    integer                                 :: i,j,Ntrim,itrim,Nocc
    !
    Ntrim=size(Ktrims,2)
    allocate(Htrims(Nso,Nso,Ntrim))
    allocate(Delta(Ntrim))
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_weyl(Ktrims(:,itrim),Nso)
       Delta(itrim)=-sign(1d0,dreal(Htrims(1,1,itrim)))
    enddo
    !
    z2=product(Delta(:))
    if(z2>0)then
       z2=0
    else
       z2=1
    end if
    !
  end function z2_number


end program weyl_3d
