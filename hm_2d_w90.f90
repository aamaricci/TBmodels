program hm_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=1024,Norb=1,Nspin=1,Nso=Nspin*Norb
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,iorb,jorb,io,ispin
  integer                                 :: ix,iy,iz,qr(15)
  real(8)                                 :: kx,ky,kz,bklen
  real(8),dimension(3) :: e1,e2,e3,bk1,bk2,bk3,d1,d2,d3
  real(8),dimension(:,:),allocatable      :: kpath
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk

  real(8),dimension(:,:,:),allocatable    :: nkgrid
  real(8)                                 :: ts,xmu,beta,eps,wmax
  real(8)                                 :: n(Nso)
  complex(8)                              :: w,Hloc(Nso,Nso)
  complex(8)                              :: Gmats(Nspin,Nspin,Norb,Norb,L)
  complex(8)                              :: Greal(Nspin,Nspin,Norb,Norb,L)
  character(len=20)                       :: file,nkstring
  logical                                 :: iexist,ibool,iener


  call parse_input_variable(nkx,"NKX","inputHM.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputHM.conf",default=500)
  call parse_input_variable(xmu,"XMU","inputHM.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputHM.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputHM.conf",default=1000.d0)
  call parse_input_variable(wmax,"WMAX","inputHM.conf",default=10d0)
  call parse_input_variable(file,"FILE","inputHM.conf",default="hkfile_bhz.in")
  call save_input_file("inputHM.conf")

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-5d0,"wini")
  call add_ctrl_var(5d0,"wfin")
  call add_ctrl_var(eps,"eps")



  !>Direct lattice basis vector
  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])



  !METHOD 1 (setup W90 --> use internal W90 model)
  call TB_w90_setup("square_2d.rHam")

  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  call TB_build_model(Hk,TB_w90_model,Nso,[Nkx,Nkx])
  Wtk = 1d0/Nktot
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"w90Hloc.dat")


  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_matsubara(Gmats,"w90Gloc",iprint=1)
  call dmft_print_gf_realaxis(Greal,"w90Gloc",iprint=1)
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb + (ispin-1)*Norb
        n(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
     enddo
  enddo

  !solve along the standard path in the 2D BZ.
  Npts=4
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_Gamma
  kpath(2,:)=kpoint_M1
  kpath(3,:)=kpoint_X1
  kpath(4,:)=kpoint_Gamma
  call TB_Solve_model(TB_w90_model,Nso,kpath,Nkpath,&
       colors_name=[red1],&
       points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
       file="w90Eigenband.nint")


  !plot observables
  open(10,file="w90observables.nint")
  write(10,"(10F20.12)")(n(iorb),iorb=1,Nso),sum(n)
  close(10)
  write(*,"(A,10F14.9)")"w90 Occupations =",(n(iorb),iorb=1,Nso),sum(n)



  call TB_w90_delete()
  Hk=zero
  Wtk=0d0
  Gmats=zero
  Greal=zero
  ts=1d0




  !METHOD 2 (setup W90 --> use TB_build TB_solve)
  call TB_w90_setup("square_2d.rHam")

  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nktot)
  call TB_build_model(Hk,Nso,[Nkx,Nkx])
  Wtk = 1d0/Nktot
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"w90bisHloc.dat")


  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_matsubara(Gmats,"w90bisGloc",iprint=1)
  call dmft_print_gf_realaxis(Greal,"w90bisGloc",iprint=1)
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb + (ispin-1)*Norb
        n(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
     enddo
  enddo

  !solve along the standard path in the 2D BZ.
  Npts=4
  kpath(1,:)=kpoint_Gamma
  kpath(2,:)=kpoint_M1
  kpath(3,:)=kpoint_X1
  kpath(4,:)=kpoint_Gamma
  call TB_Solve_model(Nso,kpath,Nkpath,&
       colors_name=[red1],&
       points_name=[character(len=20) :: 'G', 'M', 'X', 'G'],&
       file="w90bisEigenband.nint")


  !plot observables
  open(10,file="w90Bisobservables.nint")
  write(10,"(10F20.12)")(n(iorb),iorb=1,Nso),sum(n)
  close(10)
  write(*,"(A,10F14.9)")"w90bis Occupations =",(n(iorb),iorb=1,Nso),sum(n)




  Hk=zero
  Wtk=0d0
  Gmats=zero
  Greal=zero
  ts=1d0


  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nktot)
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])
  Wtk = 1d0/Nktot
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"Hloc.dat")


  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb + (ispin-1)*Norb
        n(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
     enddo
  enddo

  !solve along the standard path in the 2D BZ.
  Npts=4
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


