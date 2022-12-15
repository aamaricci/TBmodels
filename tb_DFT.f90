program tb_DFT
  USE SCIFOR
  USE DMFT_TOOLS
#ifdef _MPI
  USE MPI
#endif
  implicit none

  integer                                       :: Nlat,Norb,Nspin,Lmats,Lreal
  integer                                       :: Nktot,Nkpath,Nkvec(3),Npts,Nlso
  integer                                       :: i,j,k,ik,ilat,iorb,jorb,io,ispin
  real(8),dimension(3)                          :: e1,e2,e3
  real(8),dimension(:,:),allocatable            :: kpath
  complex(8),dimension(:,:,:),allocatable       :: Hk
  complex(8),dimension(:,:),allocatable         :: Hloc
  complex(8),dimension(:,:,:,:,:,:),allocatable :: Gmats,Greal
  real(8),allocatable                           :: Dens(:)
  character(len=60)                             :: w90file,InputFile,latfile,kpathfile,ineqfile,hkfile,OrderFile
  character(len=40),allocatable                 :: points_name(:)
  real(8)                                       :: ef,beta,xmu,eps,wini,wfin,filling
  logical                                       :: FSflag,EFflag,Spinor
  logical                                       :: master=.true.,bool
  logical                                       :: bool_hk
  logical                                       :: bool_lat
  logical                                       :: bool_kpath
  logical                                       :: bool_ineq
  logical                                       :: bool_order
  integer                                       :: unit
  integer,allocatable,dimension(:)              :: ineq_sites
  integer,dimension(3)                          :: Nin_w90
  character(len=5),dimension(3)                 :: OrderIn_w90


#ifdef _MPI
  call init_MPI
  master = get_master_MPI()
#endif


  call parse_cmd_variable(InputFile,"INPUTFILE",default="inputDFT.conf")
  call parse_input_variable(Nlat,"NLAT",InputFile,default=1,comment="Number of inequivalent sites.")
  call parse_input_variable(Norb,"NORB",InputFile,default=1,comment="Number of impurity orbitals.")
  call parse_input_variable(Nspin,"NSPIN",InputFile,default=1,comment="Number of spin degeneracy")
  call parse_input_variable(Filling,"FILLING",InputFile,default=0d0,comment="Total filling of the local problem: 0:Norb")
  call parse_input_variable(beta,"BETA",InputFile,default=1000.d0,comment="Inverse temperature, at T=0 is used as a IR cut-off.")
  call parse_input_variable(xmu,"xmu",InputFile,default=0d0,comment="Chemical potential. If filling is 0d0 this sets the total density")
  call parse_input_variable(Lmats,"LMATS",InputFile,default=1000,comment="Number of Matsubara frequencies.")
  call parse_input_variable(Lreal,"LREAL",InputFile,default=1000,comment="Number of real-axis frequencies.")
  call parse_input_variable(wini,"WINI",InputFile,default=-5.d0,comment="Smallest real-axis frequency")
  call parse_input_variable(wfin,"WFIN",InputFile,default=5.d0,comment="Largest real-axis frequency")
  call parse_input_variable(eps,"EPS",InputFile,default=0.01d0,comment="Broadening on the real-axis.")
  !
  call parse_input_variable(w90file,"w90file",InputFile,default="hij.conf")
  call parse_input_variable(hkfile,"hkfile",InputFile,default="hk.conf")
  call parse_input_variable(latfile,"latfile",InputFile,default="lat.conf")
  call parse_input_variable(kpathfile,"kpathfile",InputFile,default="kpath.conf")
  call parse_input_variable(ineqfile,"ineqfile",InputFile,default="ineq.conf")
  call parse_input_variable(orderfile,"Orderfile",InputFile,default="order.conf")
  call parse_input_variable(Spinor,"Spinor",InputFile,default=.false.)
  call parse_input_variable(FSflag,"FSflag",InputFile,default=.false.)
  call parse_input_variable(EFflag,"EFflag",InputFile,default=.false.)
  call parse_input_variable(Nkvec,"NKVEC",InputFile,default=[10,10,10])
  call parse_input_variable(nkpath,"NKPATH",InputFile,default=500)
  call print_input()
  call save_input(InputFile)
  !
  call add_ctrl_var(Nlat,"Nlat")
  call add_ctrl_var(Norb,"Norb")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")

  Nlso = Nlat*Nspin*Norb

  inquire(file=reg(hkfile),exist=bool_hk)
  inquire(file=reg(latfile),exist=bool_lat)
  inquire(file=reg(kpathfile),exist=bool_kpath)
  inquire(file=reg(ineqfile),exist=bool_ineq)
  inquire(file=reg(orderfile),exist=bool_order)


  !Get/Set Wannier ordering:
  if(bool_order)then
     open(free_unit(unit),file=reg(orderfile))
     read(unit,*)Nin_w90(1),Nin_w90(2),Nin_w90(3)
     read(unit,*)OrderIn_w90(1),OrderIn_w90(2),OrderIn_w90(3)
     close(unit)
  else
     Nin_w90    =[Nspin,Norb,Nlat]
     OrderIn_w90=[character(len=5)::"Nspin","Norb","Nlat"]
  endif

  !Setup the path in the BZ.
  if(bool_kpath)then
     Npts = file_length(reg(kpathfile))
     allocate(kpath(Npts,3))
     allocate(points_name(Npts))
     open(free_unit(unit),file=reg(kpathfile))
     do i=1,Npts
        read(unit,*)points_name(i),kpath(i,:)
     enddo
     close(unit)
  else
     Npts = 9
     allocate(kpath(Npts,3),points_name(Npts))
     kpath(1,:)=[0.5d0,0.5d0,0d0]
     kpath(2,:)=[0.5d0,0.5d0,0.5d0]
     kpath(3,:)=[0d0,0d0,0d0]
     kpath(4,:)=[0.5d0,0d0,0d0]
     kpath(5,:)=[0.5d0,0.5d0,0d0]
     kpath(6,:)=[0d0,0d0,0d0]
     kpath(7,:)=[0d0,0d0,0.5d0]
     kpath(8,:)=[0.5d0,0d0,0.5d0]
     kpath(9,:)=[0.5d0,0.5d0,0.5d0]
     points_name=[character(len=40) ::'M', 'R', '{/Symbol} G', 'X', 'M', '{/Symbol} G', 'Z','A', 'R']
  endif


  !Setup inequivalent sites in the unit cell
  allocate(ineq_sites(Nlat));ineq_sites=1
  if(bool_ineq)then
     open(free_unit(unit),file=reg(ineqfile))
     do i=1,Nlat
        read(unit,*)ineq_sites(i)
        write(*,"(A,I5,A,I5)")"Site",i,"corresponds to ",ineq_sites(i)
     enddo
     close(unit)
  endif


  !Set basis vectors:
  if(bool_lat)then
     open(free_unit(unit),file=reg(latfile))
     read(unit,*)e1
     read(unit,*)e2
     read(unit,*)e3
     close(unit)
  else
     e1 = [1d0,0d0,0d0]
     e2 = [0d0,1d0,0d0]
     e3 = [0d0,0d0,1d0]
  endif
  call TB_set_ei(e1,e2,e3)
  call TB_build_bk(verbose=.true.)


  !Setup Wannier90 or read H(k) from file:
  call start_timer
  call TB_w90_setup(reg(w90file),nlat=Nlat,norb=Norb,nspin=Nspin,Spinor=spinor,verbose=.true.)
  call stop_timer("TB_w90_setup")

  if(bool_hk)then
     call TB_read_hk(Hk,reg(hkfile),Nkvec)
     if(size(Hk,1)/=Nlat*Nspin*Norb)stop "tb_DFT error: wrong size in Hk as read from file"
  else
     call start_timer  
     call TB_w90_FermiLevel(Nkvec,filling,Ef)
     call stop_timer("TB_w90_FermiLevel")
     !
     Nktot=product(Nkvec)
     allocate(Hk(Nlso,Nlso,Nktot))
     call start_timer
     call TB_build_model(Hk,Nlso,Nkvec)
     Hk = TB_reshape_array(Hk,Nin=Nin_w90,OrderIn=OrderIn_w90,&
          OrderOut=[character(len=5):: "Norb","Nspin","Nlat"])
     call TB_write_hk(reg(hkfile),Nkvec)
     call stop_timer("TB_build_model")
  endif


  write(*,*)"Using Nk_total="//str(size(Hk,3))
  allocate(Hloc(Nlso,Nlso))
  Hloc= sum(Hk(:,:,:),dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  if(master)call TB_write_Hloc(Hloc,"w90Hloc.dat")


  
  !Solve for the renormalized bands:
  call start_timer
  if(master)call TB_Solve_model(TB_w90_model,Nlso,kpath,Nkpath,&
       colors_name=[black,red,green,blue,magenta,black,red,green,blue,magenta],&
       points_name=points_name,& 
       file="Bands_tbDFT",iproject=.true.)
  call stop_timer("get Bands")


  !Build the local GF:
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  call dmft_gloc_matsubara(Hk,Gmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  call dmft_gloc_realaxis(Hk,Greal,zeros(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  call dmft_print_gf_matsubara(Gmats,"Gloc",4)
  call dmft_print_gf_realaxis(Greal,"Gloc",4)



  if(FSflag)then
     call TB_FSurface(Nlso,0d0,Nkvec(1:2),&
          colors_name=[black,red,red,green,blue],&
          file='FS_tbDFT',cutoff=1d-1,Niter=3,Nsize=2)
  endif

  call TB_w90_delete()



#ifdef _MPI
  call finalize_MPI()
#endif



end program tb_DFT
