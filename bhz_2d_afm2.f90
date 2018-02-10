! DESCRIPTION
!   Solve the non-interacting BHZ model with AFM 2x2 basis 
!   generate the hamiltoniana matrix H(k), 

program bhz_afm
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                         :: L=2048,Norb=2,Nspin=2,Nineq=2,Nso=Nspin*Norb,Nlso=Nineq*Nso
  integer                                   :: Nk,Nktot,Nkpath,Nkx,Nky,Npts,Nx,Ny,Nlat
  integer                                   :: i,j,k,ik,iorb,jorb,ilat,jlat
  real(8),dimension(:,:),allocatable        :: kpath,Kgrid,Rgrid
  real(8),dimension(2)                      :: bk1,bk2,ei1,ei2
  real(8),dimension(2)                      :: vecK,vecRi,vecRj
  complex(8),dimension(:,:,:),allocatable   :: Hk
  complex(8),dimension(:,:,:,:),allocatable :: Hlat
  complex(8),dimension(:,:,:),allocatable   :: ftHk
  complex(8),dimension(:,:,:,:),allocatable :: ftHlat
  real(8),dimension(:),allocatable          :: Wtk
  integer,dimension(:,:),allocatable        :: Links

  real(8)                                   :: mh,rh,lambda,delta,wmax
  real(8)                                   :: xmu,beta,eps,Ekin,Eloc
  real(8),dimension(L)                      :: wm,wr
  real(8)                                   :: n(Nlso)
  complex(8)                                :: w,arg
  complex(8)                                :: Hloc(Nlso,Nlso)
  complex(8)                                :: Gmats(Nineq,Nspin,Nspin,Norb,Norb,L)
  complex(8)                                :: Greal(Nineq,Nspin,Nspin,Norb,Norb,L)
  complex(8),allocatable                    :: Gkmats(:,:,:,:,:,:)
  character(len=20)                         :: file,nkstring
  logical                                   :: iexist,ibool
  integer                                   :: isite,jsite,ispin,jspin,row,col
  !Dirac matrices:
  complex(8),dimension(4,4)                 :: Gamma1,Gamma2,Gamma3,Gamma4,Gamma5
  complex(8),dimension(2,2)                 :: Theta11,Theta12,Theta21,Theta22

  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nky,"NKY","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3.d0)
  call parse_input_variable(rh,"RH","inputBHZ.conf",default=0.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(delta,"DELTA","inputBHZ.conf",default=0.d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(wmax,"WMAX","inputBHZ.conf",default=10d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz_afm.in")
  call save_input_file("inputBHZ.conf")

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")


  Nktot=Nkx*Nky
  !
  Nx   = Nkx
  Ny   = Nkx
  Nlat = Nx*Ny

  !>Gamma matrices:
  Gamma1 = kron_pauli(pauli_z,pauli_x)
  Gamma2 =-kron_pauli(pauli_0,pauli_y)
  Gamma3 = kron_pauli(pauli_x,pauli_x)
  Gamma4 = kron_pauli(pauli_y,pauli_x)
  Gamma5 = kron_pauli(pauli_0,pauli_z)
  !
  !>Theta matrices:
  Theta11 = transpose( reshape([one,zero,zero,zero],shape(Theta11)))
  Theta12 = transpose( reshape([zero,one,zero,zero],shape(Theta12)))
  Theta21 = transpose( reshape([zero,zero,one,zero],shape(Theta21)))
  Theta22 = transpose( reshape([zero,zero,zero,one],shape(Theta22)))

  !>Direct lattice basis vector
  ei1 = [1d0,1d0]
  ei2 = [2d0,0d0]
  call TB_set_ei(ei1,ei2)

  !>Reciprocal lattice basis vector  
  bk1 = pi*[1,-1]
  bk2 = 2*pi*[0,1]
  call TB_set_bk(bk1,bk2)

  !>Get TB Hamiltonian matrix
  allocate(Hk(Nlso,Nlso,Nktot))
  allocate(Wtk(Nktot))
  call TB_build_model(Hk,hk_model,Nlso,[Nkx,Nky])
  Wtk = 1d0/Nktot
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"Hloc.dat")




  !solve along the standard path in the 2D BZ.
  Npts=4
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_Gamma
  kpath(2,:)=kpoint_X1
  kpath(3,:)=kpoint_M1
  kpath(4,:)=kpoint_Gamma
  call TB_solve_model(hk_model,Nlso,kpath,Nkpath,&
       colors_name=[red1,blue1,red1,blue1,red1,blue1,red1,blue1],&
       points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
       file="Eigenbands_afm.nint")


  !Build the local GF:  
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nineq,Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nineq,Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_matsubara(pi/beta*(2*arange(1,L)-1),Gmats,"Gloc",4)
  call dmft_print_gf_realaxis(linspace(-wmax,wmax,L),Greal,"Gloc",4)


  do i=1,Nineq
     do ispin=1,Nspin
        do iorb=1,Norb
           n(iorb+(ispin-1)*Nspin+(i-1)*Nspin*Norb) = fft_get_density(Gmats(i,ispin,ispin,iorb,iorb,:),beta)
        enddo
     enddo
  enddo

  !plot observables
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(n(iorb),iorb=1,Nlso),sum(n)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(n(iorb),iorb=1,Nlso),sum(n)


  allocate(Gkmats(Nktot,Nspin,Nspin,Norb,Norb,L))
  do ik=1,Nktot
     call dmft_gk_matsubara(Hk(:,:,ik),1d0,Gkmats(ik,:,:,:,:,:),zeros(Nspin,Nspin,Norb,Norb,L))
  enddo








contains




  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    complex(8),dimension(Nso,Nso) :: M
    complex(8),dimension(Nso,Nso) :: tx,ty,thx,thy
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = kpoint(1)
    ky = kpoint(2)
    !
    !
    M  = Mh*Gamma5
    tx = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    thx= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    !
    ty = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma2
    thy= -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma2
    !
    ! H2 =  | m1                       tx + tx^+.e^i.2.kx + ty^+.e^i.(kx+ky) + ty^+.e^i.(kx-ky) |
    !       | tx^+ + tx.e^-i.2.kx + ty.e^-i.(kx+ky)+ ty^+.e^-i.(kx-ky)          m2              |
    !
    hk(1:4,1:4)    = M
    hk(1:4,5:8)    = tx  + thx*exp(xi*2*kx) + thy*exp(xi*(kx+ky)) + ty*exp(xi*(kx-ky))
    !
    hk(5:8,1:4)    = thx + tx*exp(-xi*2*kx) + ty*exp(-xi*(kx+ky)) + thy*exp(-xi*(kx-ky))
    hk(5:8,5:8)    = M
    !
  end function hk_model


end program bhz_afm







