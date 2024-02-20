! DESCRIPTION
!   Solve the non-interacting Hubbard model with AFM 2 atoms in the basis 

program hm_afm_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none
  integer,parameter                       :: L=2048,Norb=1,Nspin=1,Nlat=2,Nso=Nspin*Norb,Nlso=Nlat*Nso
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb
  integer                                 :: ix,iy,iz
  real(8)                                 :: kx,ky,kz
  real(8)                                 :: a,a0
  real(8),dimension(:),allocatable        :: kxgrid
  real(8),dimension(:,:),allocatable      :: kpath
  real(8),dimension(2)                    :: d1,d2,d3
  real(8),dimension(2)                    :: a1,a2,a3
  real(8),dimension(2)                    :: bk1,bk2,kvec
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk
  real(8)                                 :: ts
  real(8)                                 :: xmu,beta,eps,Ekin,Eloc,wmax
  real(8)                                 :: n(Nlso)
  complex(8)                              :: Hloc(Nlso,Nlso)
  complex(8)                              :: Hlocsite(Nlat,Nso,Nso)
  complex(8)                              :: Gmats(Nlso,Nlso,L),Greal(Nlso,Nlso,L),Sfoo(Nlso,Nlso,L)
  character(len=32)                       :: file,nkstring
  logical                                 :: iexist,ibool,dcflag
  integer                                 :: icount,jcount,isite,jsite,ispin,jspin,row,col


  call parse_input_variable(nkx,"NKX","inputHM.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputHM.conf",default=500)
  call parse_input_variable(ts,"TS","inputHM.conf",default=0.5d0)
  call parse_input_variable(xmu,"XMU","inputHM.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputHM.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM.conf",default=1000.d0)
  call parse_input_variable(wmax,"WMAX","inputHM.conf",default=10d0)
  call parse_input_variable(file,"FILE","inputHM.conf",default="hkfile_bhz.in")
  call save_input_file("inputHM.conf")


  Nktot=Nkx*Nkx
  write(*,*) "Using Nk_total="//txtfy(Nktot)



  !RECIPROCAL LATTICE VECTORS:
  bk1=  pi*[ 1d0, -1d0 ]
  bk2=2*pi*[ 0d0,  1d0 ]
  call TB_set_bk(bkx=bk1,bky=bk2)

  allocate(Hk(Nlso,Nlso,Nktot))
  allocate(Wtk(Nktot))

  ik=0
  do iy=1,Nkx
     ky = dble(iy-1)/Nkx
     do ix=1,Nkx
        ik=ik+1
        kx = dble(ix-1)/Nkx
        kvec = kx*bk1 + ky*bk2
        Hk(:,:,ik) = hk_model(kvec,Nlso)
     enddo
  enddo
  Wtk = 1d0/Nktot


  Hloc=sum(Hk(:,:,:),dim=3)/Nktot


  !Build the local GF:
  Gmats=zero
  Greal=zero
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-wmax,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")
  Sfoo =zero
  call get_gloc(Hk,Gmats,Sfoo,axis='matsubara')
  Sfoo =zero
  call get_gloc(Hk,Greal,Sfoo,axis='realaxis')

  do iorb=1,Nlso
     n(iorb) = fft_get_density(Gmats(iorb,iorb,:),beta)
  enddo
  !plot observables
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(n(iorb),iorb=1,Nlso),sum(n)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(n(iorb),iorb=1,Nlso),sum(n)



  !solve along the standard path in the 2D BZ.
  Npts=4
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_Gamma
  kpath(2,:)=kpoint_M1
  kpath(3,:)=kpoint_X1
  kpath(4,:)=kpoint_Gamma
  call TB_solve_model(Hk_model,Nlso,kpath,Nkpath,&
       colors_name=[red,blue],&
       points_name=[character(len=20) :: '$\Gamma$', 'M', 'X', '$\Gamma$'],&
       file="Eigenbands_afm.nint")


contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = kpoint(1)
    ky = kpoint(2)
    !
    ! Hk =  -t * | 0                  1 + e^ikx(e^ikx + e^iky) |
    !            | 1 + e^-ikx(e^-ikx + e^-iky)   0             |
    !
    hk=zero
    hk(1,2) = -ts*(one+exp(xi*2*kx)+exp(xi*(kx+ky))+exp(xi*(kx-ky)))
    hk(2,1) = -ts*(one+exp(-xi*2*kx)+exp(-xi*(kx+ky))+exp(-xi*(kx-ky)))
  end function hk_model




end program hm_afm_2d







