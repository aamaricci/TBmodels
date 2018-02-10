program bhz_fcc_stripe
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                             :: Norb=2,Nspin=2,Ncell=2,Nso=Nspin*Norb,Ncso=Ncell*Nso
  integer                                       :: Lfreq
  integer                                       :: Nkpath,Nk,Ly,Npts,Nslat
  integer                                       :: Nx,Ny
  integer                                       :: i,j,k,ik,iorb,ilat,ispin
  integer                                       :: ix,iy,iz
  real(8)                                       :: kx,ky,kz
  real(8),dimension(:,:),allocatable            :: kpath
  real(8),dimension(2)                          :: bk1,bk2,kvec
  complex(8),dimension(:,:,:),allocatable       :: Hkr
  real(8),dimension(:),allocatable              :: Wtk
  complex(8),dimension(Nso,Nso)                 :: Gamma1,Gamma2,Gamma5
  type(rgb_color),dimension(:,:),allocatable    :: colors
  real(8)                                       :: mh,lambda,e0
  real(8)                                       :: xmu,beta,eps,wmax
  real(8),dimension(:,:,:),allocatable          :: dens
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sfoo
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gkmats
  character(len=20)                             :: file
  logical                                       :: nogf

  call parse_input_variable(nk,"NK","inputBHZ.conf",default=25)
  call parse_input_variable(Ly,"Ly","inputBHZ.conf",default=11)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(e0,"e0","inputBHZ.conf",default=1d0)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(Lfreq,"Lmats","inputBHZ.conf",default=1024)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(wmax,"wmax","inputBHZ.conf",default=-10d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=100.d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz.in")
  call parse_input_variable(nogf,"NOGF","inputBHZ.conf",default=.false.)
  call save_input_file("inputBHZ.conf")

  !if(mod(Ly,2)/=0)stop "Error. Please choose use Ly%2=0"

  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)

  bk1 = pi*[1,-1]
  bk2 = 2*pi*[0,1]


  print*,"Build H(k,R) model"
  allocate(Hkr(Ncell*Ly*Nso,Ncell*Ly*Nso,Nk))
  allocate(Wtk(Nk))
  call TB_set_bk(bk1,bk2)
  call TB_build_model(Hkr,bhz_afm2_edge_model,Ly,Ncell*Nso,Nkvec=[Nk,1,1],pbc=.false.)
  Wtk = 1d0/Nk  




  print*,"Solve H(k,R) along -pi:pi"
  Npts=3
  allocate(Kpath(Npts,1))
  kpath(1,:)=[-1]*pi
  kpath(2,:)=[ 0]*pi
  kpath(3,:)=[ 1]*pi
  allocate(colors(Ly,Ncell*Nso))
  colors = gray88
  colors(1,:) = [red1,gray88,blue1,gray88,red1,gray88,blue1,gray88]
  colors(Ly,:) =[blue1,gray88,red1,gray88,blue1,gray88,red1,gray88]
  call TB_solve_model(bhz_afm2_edge_model,Ly,Ncell*Nso,kpath,Nkpath,&
       colors_name=colors,& 
       points_name=[character(len=10) :: "-pi","0","pi"],&
       file="Eigenbands_afm.nint",pbc=.false.)



  if(nogf) stop "Called with NOGF=TRUE! Got bands now exiting."


  !Build the local GF:
  allocate(Gmats(Ncell*Ly,Nspin,Nspin,Norb,Norb,Lfreq))
  allocate(Greal(Ncell*Ly,Nspin,Nspin,Norb,Norb,Lfreq))
  Gmats=zero
  Greal=zero
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")
  call dmft_gloc_matsubara(Hkr,Wtk,Gmats,zeros(Ncell*Ly,Nspin,Nspin,Norb,Norb,Lfreq))
  call dmft_print_gf_matsubara(pi/beta*(2*arange(1,Lfreq)-1),Gmats,"Gloc",4)
  call dmft_gloc_realaxis(Hkr,Wtk,Greal,zeros(Ncell*Ly,Nspin,Nspin,Norb,Norb,Lfreq))
  call dmft_print_gf_realaxis(linspace(-wmax,wmax,Lfreq),Greal,"Gloc",4)


  allocate(Gkmats(Nk,Ncell*Ly,Nspin,Nspin,Norb,Norb,Lfreq))
  do ik=1,Nk
     call dmft_gk_matsubara(Hkr(:,:,ik),1d0,Gkmats(ik,:,:,:,:,:,:),zeros(Ncell*Ly,Nspin,Nspin,Norb,Norb,Lfreq))
  enddo


  !Build local observables:
  allocate(dens(Ncell*Ly,Nspin,Norb))
  open(10,file="density.nint")
  do iy=1,Ncell*Ly
     do ispin=1,Nspin
        do iorb=1,Norb
           dens(iy,ispin,iorb) = fft_get_density(Gmats(iy,ispin,ispin,iorb,iorb,:),beta)
        enddo
     enddo
     write(10,"(I4,1000F20.12)")iy,((dens(iy,ispin,iorb),iorb=1,Norb),ispin=1,Nspin)
  enddo
  close(10)






contains


  !BHZ on a stripe geometry;
  function bhz_afm2_edge_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx
    integer                             :: Nlat,N
    complex(8),dimension(N,N)           :: Hmat,Tmat,TmatH
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    complex(8),dimension(Nso,Nso)       :: M
    complex(8),dimension(Nso,Nso)       :: tx,ty,thx,thy
    !
    if(N/=Ncell*Nso)stop "hk_model error: N != Ncell*Nso" 
    !
    kx=kpoint(1)
    !
    M  = Mh*Gamma5
    tx = -0.5d0*e0*Gamma5 - xi*0.5d0*lambda*Gamma1
    thx= -0.5d0*e0*Gamma5 + xi*0.5d0*lambda*Gamma1
    !
    ty = -0.5d0*e0*Gamma5 - xi*0.5d0*lambda*Gamma2
    thy= -0.5d0*e0*Gamma5 + xi*0.5d0*lambda*Gamma2
    !
    Hmat(1:Nso,1:Nso)             = M
    Hmat(1:Nso,Nso+1:2*Nso)       = tx  + thx*exp(xi*2*kx)
    Hmat(Nso+1:2*Nso,1:Nso)       = thx + tx*exp(-xi*2*kx)
    Hmat(Nso+1:2*Nso,Nso+1:2*Nso) = M
    !
    Tmat(1:Nso,1:Nso)             = zero
    Tmat(1:Nso,Nso+1:2*Nso)       = thy*exp(xi*kx)
    Tmat(Nso+1:2*Nso,1:Nso)       = thy*exp(-xi*kx)
    Tmat(Nso+1:2*Nso,Nso+1:2*Nso) = zero
    !
    TmatH=conjg(transpose(Tmat))
    !
    Hrk=zero
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat
    enddo
    do i=1,Nlat-1
       Idmin=1 + (i-1)*N
       Idmax=        i*N
       Itmin=1 +     i*N
       Itmax=    (i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=Tmat
       Hrk(Itmin:Itmax,Idmin:Idmax)=TmatH
    enddo
    if(pbc)then
       Itmin=1+(Nlat-1)*N
       Itmax=0+Nlat*N
       Hrk(1:N,Itmin:Itmax)=TmatH
       Hrk(Itmin:Itmax,1:N)=Tmat
    endif
  end function bhz_afm2_edge_model



end program bhz_fcc_stripe
