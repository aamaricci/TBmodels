program LaOFeAs
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                       :: L=1024,Norb=5,Nspin=1,Nlat=2,Nlso=Nlat*Nspin*Norb
  integer                                 :: Nk,Nktot,Nkpath,Nkx,Npts
  integer                                 :: i,j,k,ik,ilat,iorb,jorb,io,ispin
  real(8),dimension(3)                    :: e1,e2,e3
  real(8),dimension(:,:),allocatable      :: kpath,kgrid
  complex(8),dimension(:,:,:),allocatable :: Hk
  real(8),dimension(:),allocatable        :: Wtk
  real(8)                                 :: xmu,beta,eps,wmax,wmin
  real(8)                                 :: n(Nlso),Ef,deltak,cutoff
  integer                                 :: FSiter
  complex(8)                              :: Hloc(Nlso,Nlso)
  complex(8)                              :: Gmats(Nlat,Nspin,Nspin,Norb,Norb,L)
  complex(8)                              :: Greal(Nlat,Nspin,Nspin,Norb,Norb,L)
  character(len=20)                       :: w90file

  call parse_input_variable(w90file,"w90file","inputLOFS.conf",default="pwscf_hr.w90")
  call parse_input_variable(Nkx,"NKX","inputLOFS.conf",default=10)
  call parse_input_variable(nkpath,"NKPATH","inputLOFS.conf",default=500)
  call parse_input_variable(xmu,"XMU","inputLOFS.conf",default=12.0906213346969d0)
  call parse_input_variable(eps,"EPS","inputLOFS.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputLOFS.conf",default=1000.d0)
  call parse_input_variable(wmin,"WMIN","inputLOFS.conf",default=-10d0)
  call parse_input_variable(wmax,"WMAX","inputLOFS.conf",default=10d0)

   call parse_input_variable(deltak,"deltak","inputLOFS.conf",default=0.13d0)
    call parse_input_variable(cutoff,"cutoff","inputLOFS.conf",default=0.2d0)
    call parse_input_variable(FSiter,"FSiter","inputLOFS.conf",default=3)
  call save_input_file("inputLOFS.conf")
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wmin,"wini")
  call add_ctrl_var(wmax,"wfin")
  call add_ctrl_var(eps,"eps")


  !>Direct lattice basis vector
  e1 = [7.625044946193571d0,0d0,0d0]
  e2 = [0d0,7.625044946193571d0,0d0]
  e3 = [0d0,0d0,16.51809612755341d0]  
  call TB_set_ei(e1,e2,e3)
  call TB_build_bk(verbose=.true.)

  !METHOD 1 (setup W90 --> use internal W90 model)
  call TB_w90_setup(reg(w90file),nlat=Nlat,nspin=Nspin,norb=Norb)
  call TB_w90_FermiLevel([Nkx,Nkx,Nkx],dble(Nlso)+1,Ef)

  !solve along a path in the 3D BZ.
  Npts = 9
  Nk=(Npts-1)*Nkpath
  allocate(kpath(Npts,3))
  kpath(1,:)=[0.5d0,0.5d0,0d0]
  kpath(2,:)=[0.5d0,0.5d0,0.5d0]
  kpath(3,:)=[0d0,0d0,0d0]
  kpath(4,:)=[0.5d0,0d0,0d0]
  kpath(5,:)=[0.5d0,0.5d0,0d0]
  kpath(6,:)=[0d0,0d0,0d0]
  kpath(7,:)=[0d0,0d0,0.5d0]
  kpath(8,:)=[0.5d0,0d0,0.5d0]
  kpath(9,:)=[0.5d0,0.5d0,0.5d0]
  call TB_Solve_model(TB_w90_model,Nlso,kpath,Nkpath,&
       colors_name=[red,green,blue,magenta,black,red,green,blue,magenta,black],&
       points_name=[character(len=40) ::'M', 'R', 'G', 'X', 'M', 'G', 'Z','A', 'R'],&
       file="w90Eigenband.nint",iproject=.true.)

  
  call TB_fsurface(TB_w90_model,Nlso,xmu,[Nkx,Nkx],&
       colors_name=[red,green,blue,magenta,black,red,green,blue,magenta,black],&
       file="LaOFeAs_FSurf",&
       cutoff=cutoff,deltak=deltak,max_order=FSiter,iwrite=.true.)


  stop
  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:
  Nktot=Nkx*Nkx*Nkx ;   write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nlso,Nlso,Nktot))
  allocate(Wtk(Nktot))
  call TB_build_model(Hk,TB_w90_model,Nlso,[Nkx,Nkx,Nkx])
  Wtk = 1d0/Nktot
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc,"w90Hloc.dat")


  call get_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
  call get_gloc_realaxis(Hk,Wtk,Greal,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
  call print_gf_matsubara(Gmats,"w90Gloc",iprint=4)
  call print_gf_realaxis(Greal,"w90Gloc",iprint=4)
  do ilat=1,Nlat
     do ispin=1,Nspin
        do iorb=1,Norb
           io = iorb + (ispin-1)*Norb + (ilat-1)*Nspin*Norb
           n(io) = fft_get_density(Gmats(Nlat,ispin,ispin,iorb,iorb,:),beta)
        enddo
     enddo
  enddo



  !plot observables
  open(10,file="w90observables.nint")
  write(10,"(20F20.12)")(n(iorb),iorb=1,Nlso),sum(n)
  close(10)
  write(*,"(A,20F14.9)")"w90 Occupations =",(n(iorb),iorb=1,Nlso),sum(n)


  call TB_w90_delete()


end program LaOFeAs


