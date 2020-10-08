program kanemele
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                             :: Norb=1,Nspin=2,Nlat=2,Nlso=Nlat*Nspin*Norb
  integer                                       :: Nk,Lk,Nkpath,Npts,L

  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk
  complex(8),allocatable,dimension(:,:)         :: kmHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),allocatable,dimension(:)              :: Wtk

  real(8),dimension(2)                          :: e1,e2   !real-space lattice basis
  real(8),dimension(2)                          :: bk1,bk2 !reciprocal space lattice basis
  real(8),dimension(2)                          :: d1,d2,d3
  real(8),dimension(2)                          :: a1,a2,a3
  real(8),dimension(2)                          :: pointK,pointKp,bklen

  !variables for the model:
  real(8)                                       :: t1,t2,phi,Delta,xmu,beta,eps
  character(len=32)                             :: finput
  integer                                       :: io,ispin,ilat
  real(8),dimension(:,:),allocatable            :: KPath
  complex(8),dimension(:,:,:,:,:,:),allocatable :: Greal,Gmats
  real(8)                                       :: dens(Nlso)


  !Parse additional variables && read Input && read H(k)^2x2
  call parse_cmd_variable(finput,"FINPUT",default='inputKM.conf')
  call parse_input_variable(nk,"NK",finput,default=25)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(L,"L",finput,default=2048)
  call parse_input_variable(t1,"T1",finput,default=1d0)
  call parse_input_variable(t2,"T2",finput,default=0d0)
  call parse_input_variable(phi,"PHI",finput,default=pi/2d0)
  call parse_input_variable(Delta,"Delta",finput,default=0d0)
  call parse_input_variable(xmu,"XMU",finput,default=0.d0)
  call parse_input_variable(eps,"EPS",finput,default=4.d-2)
  call parse_input_variable(beta,"BETA",finput,default=1000.d0)
  call save_input_file(finput)

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")


  !Lattice basis (a=1; a0=sqrt3*a) is:
  !e_1 = a0 [ sqrt3/2 , 1/2 ] = 3/2a[1, 1/sqrt3]
  !e_2 = a0 [ sqrt3/2 ,-1/2 ] = 3/2a[1,-1/sqrt3]
  e1 = 3d0/2d0*[1d0, 1d0/sqrt(3d0)]
  e2 = 3d0/2d0*[1d0,-1d0/sqrt(3d0)]

  !LATTICE BASIS: nearest neighbor: A-->B, B-->A
  d1= [  1d0/2d0 , sqrt(3d0)/2d0 ]
  d2= [  1d0/2d0 ,-sqrt(3d0)/2d0 ]
  d3= [ -1d0     , 0d0           ]

  !next nearest-neighbor displacements: A-->A, B-->B, cell basis
  a1 = d1-d3                    !3/2*a[1,1/sqrt3]
  a2 = d2-d3                    !3/2*a[1,-1/sqrt3]
  a3 = d1-d2

  pointK = [2*pi/3, 2*pi/3/sqrt(3d0)]
  pointKp= [2*pi/3,-2*pi/3/sqrt(3d0)]



  !RECIPROCAL LATTICE VECTORS:
  bklen=2d0*pi/3d0
  bk1=bklen*[ 1d0, sqrt(3d0)] 
  bk2=bklen*[ 1d0,-sqrt(3d0)]
  call TB_set_bk(bkx=bk1,bky=bk2)


  !Build the Hamiltonian on a grid or on path
  Lk= Nk*Nk
  write(*,*)"Build H(k) Kane-Mele:",Lk
  write(*,*)"# of SO-bands     :",Nlso
  !
  allocate(Hk(Nlso,Nlso,Lk));Hk=zero
  allocate(wtk(Lk));Wtk=0d0
  !
  !
  call TB_build_model(Hk,hk_kanemele_model,Nlso,[Nk,Nk],wdos=.false.)
  Wtk = 1d0/Lk
  !
  !
  allocate(kmHloc(Nlso,Nlso))
  kmHloc = sum(Hk(:,:,:),dim=3)/Lk
  where(abs(dreal(kmHloc))<1.d-6)kmHloc=0d0
  call TB_write_Hloc(kmHloc)
  call TB_write_Hloc(kmHloc,'Hloc.txt')
  !
  !
  Npts = 4
  allocate(Kpath(Npts,2))
  KPath(1,:)=[0d0,0d0]
  KPath(2,:)=pointK
  Kpath(3,:)=pointKp
  KPath(4,:)=[0d0,0d0]
  call TB_Solve_model(hk_kanemele_model,Nlso,KPath,Nkpath,&
       colors_name=[red1,blue1,red1,blue1],&
       points_name=[character(len=10) :: "G","K","K`","G"],&
       file="KMbands.nint",iproject=.false.)
  !
  !Build the local GF:
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,L))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4)
  !


  !Occupation:
  do ilat=1,Nlat
     do ispin=1,Nspin
        io = ispin + (ilat-1)*Nspin
        dens(io) = fft_get_density(Gmats(ilat,ispin,ispin,1,1,:),beta)
     enddo
  enddo
  write(*,"(A,20F14.9)")"Occupations =",(dens(io),io=1,Nlso),sum(dens)



contains



  !--------------------------------------------------------------------!
  !Kane-Mele HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_kanemele_model(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(2,2)       :: hk11,hk22
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz
    real(8)                         :: kdotd(3),kdota(3)
    !(k.d_j)
    kdotd(1) = dot_product(kpoint,d1)
    kdotd(2) = dot_product(kpoint,d2)
    kdotd(3) = dot_product(kpoint,d3)
    !(k.a_j)
    kdota(1) = dot_product(kpoint,a1)
    kdota(2) = dot_product(kpoint,a2)
    kdota(3) = dot_product(kpoint,a3)
    !
    h0 = 2*t2*cos(phi)*sum( cos(kdota(:)) )
    hx =-t1*sum( cos(kdotd(:)) )
    hy =-t1*sum( sin(kdotd(:)) )
    hz = 2*t2*sin(phi)*sum( sin(kdota(:)) )
    !
    hk11 = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z + Delta*pauli_z
    !
    hk22 = h0*pauli_0 + hx*pauli_x - hy*pauli_y - hz*pauli_z + Delta*pauli_z
    !
    hk          = zero
    hk(1:2,1:2) = hk11
    hk(3:4,3:4) = hk22
    !
  end function hk_kanemele_model














end program kanemele



