program ed_W90
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  !
  !#############################################
  !#                                           #
  !#       THIS CODE IS MPI COMPILED ONLY      #
  !#                                           #
  !#############################################
  !
  !
  !#########   VARIABLEs DECLARATION   #########
  !
  integer                                             :: iloop,unit1
  integer                                             :: i,j,io,jo,ndx
  integer                                             :: iorb,jorb
  integer                                             :: ispin,jspin
  integer                                             :: ilat,jlat,Nlat
  integer                                             :: ifreq,Lfreq
  integer                                             :: ilayer,Nlayer
  integer                                             :: NlNsNo
  logical                                             :: converged,converged1,converged2,converged3
  real(8)                                             :: wmixing
  real(8)                                             :: xmu_tmp
  real(8),allocatable,dimension(:)                    :: dens
  character(len=60)                                   :: finput
  character(len=32)                                   :: geometry
  character(len=32)                                   :: z_symmetry
  !Mpi:
  integer                                             :: comm,rank,ier
  logical                                             :: master
  !Bath:
  integer                                             :: Nb
  real(8)   ,allocatable,dimension(:)                 :: Bath,Bath_
  !Hamiltoninas:
  integer                                             :: ik,Nk,Lk,Nkpath
  complex(8),allocatable,dimension(:,:,:)             :: Hk
  complex(8),allocatable,dimension(:,:,:,:)         :: Hloc
  real(8)   ,allocatable,dimension(:)                 :: Wtk
  !local dmft Weisss:
  complex(8),allocatable,dimension(:,:,:,:,:)       :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:)       :: Gmats,Greal,dos_Greal
  complex(8),allocatable,dimension(:,:,:,:,:)       :: Weiss,Weiss_old
  !Irreducible dmft Weisss:
  !meshes:
  real(8)                                             :: dw
  real(8)   ,allocatable,dimension(:)                 :: wr,wm
  !convergence test:
  integer                                             :: memory
  complex(8),allocatable,dimension(:)                 :: conv_funct
  !custom variables for chempot search:
  logical                                             :: converged_check
  integer                                             :: conv_n_loop=0
  real(8)                                             :: sumdens,xmu_old
  real(8)   ,allocatable,dimension(:,:)               :: orb_dens_lat,orb_mag_lat
  real(8)   ,allocatable,dimension(:)                 :: orb_dens_single,orb_mag_single
  logical                                             :: look4n=.true.
  !custom variables misc:
  logical                                             :: computeG0loc
  logical                                             :: bulk_magsym
  complex(8),allocatable,dimension(:,:)               :: zeta
  complex(8),allocatable,dimension(:,:,:,:,:)         :: Nmatrix_nn
  complex(8),allocatable,dimension(:,:,:)             :: Gloc
  !Ek calculation:
  logical                                             :: computeEk
  real(8)                                             :: Ek
  real(8)   ,allocatable,dimension(:)                 :: Ekm
  complex(8),allocatable,dimension(:,:,:,:,:,:,:)     :: Gkmats
  !fake flags
  !For printing DOS
  character(len=128)                            :: dos_file="dos_Hk"

  !#########   MPI INITIALIZATION   #########
  !
  call init_MPI(ier)
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  master = get_Master_MPI(comm)
  rank = get_Rank_MPI(comm)
  write(*,*) "ier=",ier
  write(*,*) "master=",master
  write(*,*) "rank=",rank
  !
  !
  !#########    VARIABLE PARSING    #########
  !
  !
  call parse_cmd_variable(finput,           "FINPUT",              default='inputED_W90.in')
  call parse_input_variable(nk,             "NK",finput,           default=10)
  call parse_input_variable(NLAT,           "NLAT",finput,         default=1)
  call parse_input_variable(nkpath,         "NKPATH",finput,       default=20)
  call parse_input_variable(wmixing,        "WMIXING",finput,      default=0.5d0)
  call parse_input_variable(computeG0loc,   "COMPUTEG0loc",finput, default=.false.)
  call parse_input_variable(geometry,       "GEOMETRY",finput,     default="bulk")
  call parse_input_variable(z_symmetry,     "ZSYMMETRY",finput,    default="FERRO")
  call parse_input_variable(bulk_magsym,    "BULKMAGSYM",finput,   default=.false.)
  call parse_input_variable(memory,         "MEMORY",finput,       default=3)
  call parse_input_variable(computeEk,      "COMPUTEEK",finput,    default=.false.)
  !
  call ed_read_input(trim(finput),comm)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !
  geometry=reg(geometry)
  z_symmetry=reg(z_symmetry)
  write(*,*) "Basic inputs prints start"
  write(*,*) "geometry=",geometry
  write(*,*) "bulk_magsym=",bulk_magsym
  write(*,*) "Nlat=",Nlat
  !
  !
  Nlayer=Nlat/2

  NlNsNo=Nspin*Norb
  write(*,*) "Nk=",Nk
  write(*,*) "Nmats=",Lmats
  write(*,*) "Nreal=",Lreal
  write(*,*) "Norb=",Norb
  write(*,*) "Nspin=",Nspin
  write(*,*) "Nlat=",Nlat
  write(*,*) "Nlayer=",Nlayer
  write(*,*) "NlNsNo=",NlNsNo
  write(*,*) "Basic inputs prints end"
  !
  !
  !##################       ALLOCATION       ##################
  !
  !
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));                 Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));                 Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));                 Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));                 Greal=zero
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));                 Weiss=zero
  allocate(Weiss_old(Nspin,Nspin,Norb,Norb,Lmats));             Weiss_old=zero
  !
  allocate(Hk(NlNsNo,NlNsNo,Nk*Nk*Nk));                              Hk=zero
  allocate(Hloc(Nspin,Nspin,Norb,Norb));                    Hloc=zero
  allocate(Wtk(Nk*Nk*Nk));                                           Wtk=1.d0/(Nk*Nk*Nk)
  !
  allocate(wr(Lreal));wr=0.0d0;                                      wr=linspace(wini,wfin,Lreal,mesh=dw)
  allocate(wm(Lmats));wm=0.0d0;                                      wm = pi/beta*real(2*arange(1,Lmats)-1,8)
  !
  write(*,*) "Basic allocation prints start"
  write(*,*) "Smats=",shape(Smats)
  write(*,*) "Gmats=",shape(Gmats)
  write(*,*) "Sreal=",shape(Sreal)
  write(*,*) "Greal=",shape(Greal)
  write(*,*) "feild=",shape(Weiss)
  write(*,*) "feild_old=",shape(Weiss_old)
  write(*,*) "Hk=",shape(Hk)
  write(*,*) "Hloc=",shape(Hloc)
  write(*,*) "Wtk=",shape(Wtk)
  write(*,*) "Basic allocation prints end"
  !
  !
  !
  !##################        BUILD Hk        ##################
  !
  call read_myhk("W90_hr_bulk.dat","Hk.dat","Hloc","Kpoints.dat")
  !
  write(*,*) "check DFT"

  allocate(dos_Greal(Nspin,Nspin,Norb,Norb,Lreal));  dos_Greal=zero
  call dmft_gloc_realaxis(Hk,Wtk,dos_Greal,zeros(1,1,Norb,Norb,Lreal))
  call dmft_print_gf_realaxis(dos_Greal(:,:,:,:,:),trim(dos_file),iprint=1)
  !call dmft_gloc_realaxis(comm,Hk,Wtk,Greal,Sreal,mpi_split='k')
  !if(master)call build_eigenbands("W90_hr_bulk.dat","Bands.dat","Hk_path.dat","Kpoints_path.dat",zeros(1,1,Norb,Norb,Lreal))
  deallocate(dos_Greal)
  !stop("test")
  !
  !##################          BATH          ##################
  Nb=get_bath_dimension(Hloc(:,:,:,:))
  if(master)write(LOGfile,*)"   Bath_size: ",Nb," layers: ",Nlayer
  allocate(Bath(Nb));    Bath=0.0d0
  allocate(Bath_(Nb));    Bath_=0.0d0
  write(*,*) "Nb=",Nb
  write(*,*) "Bath=",shape(Bath)
  write(*,*) "masterm=",master
  write(*,*) "cg_scheme=",cg_scheme
  write(*,*) "computeEk=",computeEk
  write(*,*) "ed_verbose=",ed_verbose
  !
  !
  !##################      INIT SOLVER       ##################
  write(*,*) "INIT SOLVER start by chhcl"
  !
  write(*,*) "Singlecell initialization by chhcl"
  call ed_init_solver(comm,Bath,Hloc(:,:,:,:))
  !
  write(*,*) "INIT SOLVER end by chhcl"
  !
  !
  !##################          DMFT          ##################
  !
  iloop=0 ; converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     if(master)call start_loop(iloop,nloop,"DMFT-loop")
     !
     !--------  solve impurity (diagonal basis)  --------
     call ed_solve(comm,Bath,Hloc(:,:,:,:))
     !
     !--------    get sigmas   (diagonal basis)  --------
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_real(Sreal)

     !------  get local Gf's  (Wannier90 basis)  ------
     call dmft_gloc_matsubara(comm,Hk,Wtk,Gmats,Smats,mpi_split='k')
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=1)

     !
     !------    get Weiss     (Wannier90 basis)  --------
     write(*,*) "get Weiss by chhcl"
     if(master)then
        if(cg_scheme=='weiss')then
           call dmft_weiss(Gmats,Smats,Weiss,Hloc)
           call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1)
        elseif(cg_scheme=='delta')then
           call dmft_delta(Gmats,Smats,Weiss,Hloc)
           call dmft_print_gf_matsubara(Weiss,"Delta",iprint=1)
        endif
     endif
     call Bcast_MPI(comm,Weiss)
     !------    fit Weiss     ------
     call ed_chi2_fitgf(comm,Weiss(:,:,:,:,:),Bath,ispin=1)
     !
     !
     !------    mix Weiss     ------
     if(master)then
        write(*,*) "wmixing=",wmixing
        write(*,*) "Weiss_old=",shape(Weiss_old)
        write(*,*) "Weiss=",shape(Weiss)
        write(*,*) "mix Weiss by chhcl"
        write(*,*) "Bath old=",Bath_
     endif
     if(iloop>1)then
        Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_old
        Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
     endif
     Weiss_old=Weiss
     Bath_=Bath
     !
     !
     if(master)then
        write(*,*) "Bath new=",Bath
        write(*,*)"Norb=",Norb
        write(*,*)"Nlat=",Nlat
        write(*,*)"Nlayer=",Nlayer
        write(*,*) "nread=",nread

        if(nread/=0d0)then
           allocate(dens(Norb));dens=0d0
           call ed_get_dens(dens)
           sumdens=sum(dens)
           write(*,*) "Before xmu=",xmu
           write(*,*) "Before sumdens=",sumdens
           call search_chemical_potential(xmu,sumdens,converged)
           write(*,*) "After xmu=",xmu
           unit1=12
           open(unit=unit1,file="xmu.restart",form="FORMATTED",status="REPLACE",action="WRITE")
           write(unit=unit1,fmt=*) xmu
           close(unit=unit1)    
           deallocate(dens)
        endif
        !
        !
        write(LOGfile,*) "   ------------------- convergence --------------------"
        converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.,index=1,total=3)
        write(*,*) "Inside  converged=",converged
        !converged = converged .AND. check_convergence(Weiss(1,1,2,2,:),dmft_error,nsuccess,nloop,reset=.false.,index=2,total=3)
        !converged = converged .AND. check_convergence(Weiss(1,1,3,3,:),dmft_error,nsuccess,nloop,reset=.false.,index=3,total=3)
        !converged1 = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.,index=1,total=3)
        !converged2 = check_convergence(Weiss(1,1,2,2,:),dmft_error,nsuccess,nloop,reset=.false.,index=2,total=3)
        !converged3 = check_convergence(Weiss(1,1,3,3,:),dmft_error,nsuccess,nloop,reset=.false.,index=3,total=3)
        !converged  = converged1.AND.converged2
        !converged  = converged1.AND.converged2 .AND. converged3
     endif
     !
     !---- each loop operations -----
     !call Bcast_MPI(comm,Bath)
     call MPI_Barrier(comm,ier)
     call Bcast_MPI(comm,xmu)
     !call Bcast_MPI(comm,converged1)
     !call Bcast_MPI(comm,converged2)
     !call Bcast_MPI(comm,converged3)
     call Bcast_MPI(comm,converged)
     !write(*,*) "rank=",rank," iloop=",iloop
     !write(*,*) "rank=",rank," converged=",converged
     !write(*,*) "rank=",rank," converged1=",converged1
     !write(*,*) "rank=",rank," converged2=",converged2
     !write(*,*) "rank=",rank," converged3=",converged3
     !write(*,*)"ier=",ier
     !
     if(master)call end_loop
     !
  enddo
  deallocate(Bath,Bath_)
  !
  !
  !##################    POST-PROCESSING     ##################

  !call dmft_kinetic_energy(comm,Hk,Wtk,Smats)
  !write(*,*) "sigma_matsubara"
  !call ed_get_sigma_matsubara(Smats)
  !write(*,*) "sigma_real"
  !call ed_get_sigma_real(Sreal)
  !write(*,*) "dmft_gloc"
  call dmft_gloc_realaxis(comm,Hk,Wtk,Greal,Sreal,mpi_split='k')
  !if(master)call dmft_print_gf_realaxis(Greal,"Gloc",iprint=1)
  !call dmft_gloc_realaxis(comm,Hk,Wtk,Greal,Sreal,mpi_split='k')
  if(master)call dmft_print_gf_realaxis(Greal,"Greal",iprint=1)
  !call dmft_kinetic_energy(comm,Hk,Wtk,Smats)
  !
  !write(*,*) "Final"
  !call MPI_Barrier(comm,ier)
  !write(*,*)"ier=",ier
  !------   compute Bands  ------
  if(master.and.geometry=="bulk")call build_eigenbands("W90_hr_bulk.dat","Bands.dat","Hk_path.dat","Kpoints_path.dat",Sreal)
  !
  !master = get_Master_MPI(comm)
  !rank = get_Rank_MPI(comm)
  !write(*,*) "master=",master
  !write(*,*) "rank=",rank
  call finalize_MPI()
  !write(*,*) "After Final"
  !
  !
contains
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Read the Non interacting Hamiltonian from  file
  !         Also this is just for testing the correct interface with the
  !         translator of the W90 output
  !         The re-ordering part can be used or not, depending on what the user of W90 did.
  !+------------------------------------------------------------------------------------------+!
  subroutine read_myhk(fileHR,fileHk,fileHloc,fileKpoints,Nlat_notfake_)
    implicit none
    character(len=*)            ,intent(in)           ::   fileHR
    character(len=*)            ,intent(in)           ::   fileHk
    character(len=*)            ,intent(in)           ::   fileHloc
    character(len=*)            ,intent(in)           ::   fileKpoints
    integer                     ,intent(in),optional  ::   Nlat_notfake_
    integer                                           ::   Nlat_notfake,dim_notfake,dim_plane
    logical                                           ::   IOfile
    real(8)                                           ::   mu,bk_x(3),bk_y(3),bk_z(3)
    real(8)         ,allocatable                      ::   Aw(:,:)
    integer         ,allocatable                      ::   Nkvec(:)
    real(8)         ,allocatable                      ::   Kvec(:,:)
    complex(8)      ,allocatable                      ::   Hloc_w90(:,:)
    complex(8)      ,allocatable                      ::   Hk_tmp(:,:,:)
    complex(8)      ,allocatable                      ::   Potential_so(:,:)
    complex(8)      ,allocatable                      ::   Potential_nn(:,:,:,:,:)
    complex(8)      ,allocatable                      ::   Gmats(:,:,:),Greal(:,:,:)
    complex(8)      ,allocatable                      ::   Gso(:,:,:,:,:,:)
    !
    !
    write(*,*) "Read my HK"
    bk_x = [1.d0,0.d0,0.d0]*2*pi
    bk_y = [0.d0,1.d0,0.d0]*2*pi
    bk_z = [0.d0,0.d0,1.d0]*2*pi
    call TB_set_bk(bk_x,bk_y,bk_z)
    !
    !
    Lk=Nk*Nk*Nk
    if(master)write(LOGfile,*)" Bulk tot k-points:",Lk
    !
    !DEBUG>>
    Nlat_notfake=Nlat
    if(present(Nlat_notfake_))then
       Nlat_notfake=Nlat_notfake_
       dim_plane=2*Nspin*Norb
       dim_notfake=Nlat_notfake*Nspin*Norb
       deallocate(Hk)
       allocate(Hk(dim_notfake,dim_notfake,Nk*Nk*Nk));Hk=zero
    endif
    if(master) write(*,*)"dim_notfake=",dim_notfake," shape(Hk)=",shape(Hk)
    !>>DEBUG
    !
    allocate(Hloc_w90(NlNsNo,NlNsNo));Hloc_w90=zero
    allocate(Kvec(Lk,3));Kvec=0d0
    allocate(Nkvec(3));Nkvec=0
    Nkvec=[Nk,Nk,Nk]
    !
    inquire(file=fileHk,exist=IOfile)
    !
    write(*,*) "Hk=",fileHK
    write(*,*) "inquire for Hk=",IOfile
    if(IOfile)then
       write(LOGfile,*) " Reading existing Hk from:  ",fileHk
       call TB_read_hk(Hk,fileHk,Nspin*Norb*Nlat_notfake,1,1,Nlat,Nkvec,Kvec)
    else
       write(LOGfile,*) " Transforming HR from:  ",fileHR
       !call TB_hr_to_hk(Hk,fileHR,Nspin,Norb,Nlat_notfake,Nkvec,Kvec,fileHk,fileKpoints)
       call TB_hr_to_hk([1.d0,0.d0,0.d0],[0.d0,1.d0,0.d0],[0.d0,0.d0,1.d0],Hk,Hloc_w90,fileHR,Nspin,Norb,  &
    Nlat_notfake,Nkvec,Kvec,fileHk,fileKpoints)
    endif

    Hloc_w90=sum(Hk,dim=3)/Lk
    Hloc(1,1,1,1)=Hloc_w90(1,1)
    Hloc(1,1,2,2)=Hloc_w90(2,2)
    Hloc(1,1,3,3)=Hloc_w90(3,3)
    !
    !DEBUG>>
    !
  end subroutine read_myhk
  !
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: solve H(k) along path in the BZ.
  !+------------------------------------------------------------------------------------------+!
  subroutine build_eigenbands(fileHR,fileband,fileHk_path,fileKpoints_path,Sreal_)
    implicit none
    character(len=*)            ,intent(in)           ::   fileHR
    character(len=*)            ,intent(in),optional  ::   fileband,fileHk_path,fileKpoints_path
    complex(8)      ,intent(in),optional  ::   Sreal_(:,:,:,:,:)
    integer                                           ::   Npts,Nkpathread
    real(8)         ,allocatable                      ::   kpath(:,:),kgrid(:,:),Scorr(:,:)
    complex(8)      ,allocatable                      ::   Hkpath(:,:,:)
    complex(8)      ,allocatable                      ::   Gkreal(:,:,:,:,:,:)
    complex(8)      ,allocatable                      ::   Gkr(:,:,:,:)
    type(rgb_color) ,allocatable                      ::   colors(:),colors_orb(:)
    logical                                           ::   IOfile
    !
    !
    allocate(colors(NlNsNo))
    allocate(colors_orb(Norb))
    colors_orb=[red1,green1,blue1]
    do i=1,Nspin
       colors(1+(i-1)*Norb:Norb+(i-1)*Norb)=colors_orb
    enddo
    !
    write(LOGfile,*)
    write(LOGfile,*)"Build bulk H(k) along the path M-R-G-M-X-G-X"
    write(LOGfile,*)
    Npts = 8
    Lk=(Npts-1)*Nkpath
    allocate(kpath(Npts,3))
    kpath(1,:)=[0,0,0]*pi
    kpath(2,:)=[1,0,0]*pi
    kpath(3,:)=[1,1,0]*pi
    kpath(4,:)=[0,0,0]*pi
    kpath(5,:)=[1,1,1]*pi
    kpath(6,:)=[1,0,1]*pi
    kpath(7,:)=[0,0,1]*pi
    kpath(8,:)=[0,0,0]*pi
    !
    allocate(kgrid(Lk,3))  ;kgrid=0d0
    allocate(Hkpath(NlNsNo,NlNsNo,Lk));Hkpath=zero
    allocate(Gkreal(Lk,Nspin,Nspin,Norb,Norb,Lreal));Gkreal=zero
    allocate(Scorr(NlNsNo,NlNsNo));Scorr=0d0
    if(present(Sreal_))Scorr=real(nn2so_reshape(Sreal_(:,:,:,:,1),Nspin,Norb),8)
    if(master)write(*,*)"Scorr=",Scorr
    !
    inquire(file=fileHk_path,exist=IOfile)
    !
    write(*,*)"fileHk_path=",fileHk_path
    write(*,*)"inquire for fileHk_path=",IOfile
    if(IOfile)then
       write(LOGfile,*) "   Reading existing Hkpath on: ",fileHk_path

       call TB_read_hk(Hkpath,fileHk_path,Nspin*Norb,Nkpathread,kpath,kgrid)
       if(Nkpathread.ne.Nkpath) stop "Eigenbands wrong Nkpath readed"
    else
       write(LOGfile,*) "   Solving model on path in main"
       call TB_solve_model(   fileHR,Nspin,Norb,Nlat,kpath,Nkpath,colors               &
            ,   [character(len=20) ::'G', 'X', 'M', 'G', 'R', 'A', 'Z','G']  &
            ,   fileband                                                 &
            ,   fileHk_path                                              &
            ,   fileKpoints_path                                         &
            ,   Scorr                                                    &
            ,   Hkpath                                                   &
            ,   kgrid                                                    )
    endif
    !
    allocate(Gkr(Lk,Nspin*Norb,Nspin*Norb,Lreal));Gkr=zero
    !do ik=1,Lk
    !   call dmft_gk_realaxis(Hkpath(:,:,ik),1.d0/Lk,Gkreal(ik,:,:,:,:,:,:),Sreal)
    !   !faccio questa cosa qui sotto per separare per bene i due blocchi di spin
    !   do ifreq=1,Lreal
    !      Gkr(ik,:,:,ifreq)=matmul(P,matmul(nn2so_reshape(Gkreal(ik,:,:,:,:,:,ifreq),Nspin,Norb),transpose(P)))
    !   enddo
    !enddo
    !!
    !open(unit=106,file='Akw_s1.dat',status='unknown',action='write',position='rewind')
    !open(unit=107,file='Akw_s2.dat',status='unknown',action='write',position='rewind')
    !do ifreq=1,Lreal
    !   write(106,'(9000F18.12)')wr(ifreq),(30.*trace(-aimag(Gkr(ik,1:Nlat*Norb,1:Nlat*Norb,ifreq))/pi)+10.*ik,ik=1,Lk)
    !   write(107,'(9000F18.12)')wr(ifreq),(30.*trace(-aimag(Gkr(ik,1+Nlat*Norb:Nlat*Nspin*Norb,1+Nlat*Norb:Nlat*Nspin*Norb,ifreq))/pi)+10.*ik,ik=1,Lk)
    !enddo
    !close(106)
    !close(107)
    write(LOGfile,*)"Im done on the path in main"
    !
    deallocate(kgrid,Hkpath,Gkreal,Scorr,Gkr,colors,colors_orb)
    !
  end subroutine build_eigenbands
  !
  !
  !
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Put in proper order the input coming from w90. Depends on the w90 users.
  !+------------------------------------------------------------------------------------------+!
  !
  !
  !
  !____________________________________________________________________________________________!
  !                                       Gfs
  !____________________________________________________________________________________________!
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: G0_loc functions
  !+------------------------------------------------------------------------------------------+!
  function inverse_g0k(iw,hk_,Nlat,mu_) result(g0k)
    implicit none
    complex(8),intent(in)                                  :: iw
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)  :: hk_
    real(8),intent(in),optional                            :: mu_
    real(8)                                                :: mu
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)  :: g0k,g0k_tmp
    integer                                                :: i,ndx,Nlat
    integer (kind=4), dimension(6)                         :: ipiv
    integer (kind=1)                                       :: ok
    integer (kind=4), parameter                            :: lwork=2000
    complex (kind=8), dimension(lwork)                     :: work
    real    (kind=8), dimension(lwork)                     :: rwork
    !
    mu=0.d0
    if(present(mu_))mu=mu_
    g0k=zero;g0k_tmp=zero
    !
    g0k=(iw+mu)*eye(Nlat*Nspin*Norb)-hk_
    g0k_tmp=g0k
    !
    call inv(g0k)
    call inversion_test(g0k,g0k_tmp,1.e-6,Nlat)
  end function inverse_g0k
  !
  !
  !
  !____________________________________________________________________________________________!
  !                                     utilities
  !____________________________________________________________________________________________!
  !+------------------------------------------------------------------------------------------+!
  !PURPOSE: Inversion test
  !+------------------------------------------------------------------------------------------+!
  subroutine inversion_test(A,B,tol,Nlat)
    implicit none
    integer (kind=4), intent(in)   ::   Nlat
    complex (kind=8), intent(in)   ::   A(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    complex (kind=8), intent(in)   ::   B(Nlat*Nspin*Norb,Nlat*Nspin*Norb)
    real    (kind=4), intent(in)   ::   tol
    real    (kind=4)               ::   error
    integer (kind=2)               ::   dime

    if (size(A).ne.size(B)) then
       write(LOGfile,*) "Matrices not equal cannot perform inversion test"
       stop
    endif
    dime=maxval(shape(A))
    error=abs(float(dime)-real(sum(matmul(A,B))))
    if (error.gt.tol) write(LOGfile,*) "inversion test fail",error
  end subroutine inversion_test
  !
  !
  !
end program ed_W90








!    kpath( 1,:)=kpoint_Gamma
!    kpath( 2,:)=kpoint_X1
!    kpath( 3,:)=kpoint_M1
!    kpath( 4,:)=kpoint_X2
!    kpath( 5,:)=kpoint_Gamma
!    kpath( 6,:)=kpoint_X3
!    kpath( 7,:)=kpoint_M3
!    kpath( 8,:)=kpoint_R
!    kpath( 9,:)=kpoint_M2
!    kpath(10,:)=kpoint_X1
!    kpath(11,:)=kpoint_M1
!    kpath(12,:)=kpoint_X2
!    kpath(13,:)=kpoint_Gamma
!    kpath(14,:)=kpoint_X3
!    kpath(15,:)=kpoint_M3
!    kpath(16,:)=kpoint_R
!    kpath(17,:)=kpoint_M2
!    kpath(18,:)=kpoint_X3
