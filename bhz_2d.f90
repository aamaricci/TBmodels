! DESCRIPTION
!   Solve the non-interacting BHZ model and generate the hamiltoniana matrix H(k), 
!   Model parameters are: t_dd,v0/tpd,ep0,ed0.
!   Using the 2d square lattice dispersion. 
!   The matrix is written in the Wannier90 form, as expected by w2CTQMC & ggED code.
!   Output is on file *hkfile.in default.
!
! MODEL Hamiltonian is:
!
! |     h^{2x2}(k)              &         hso^{2x2}(k)        |
! |      [hso^{2x2}]*(k)        &        [h^{2x2}]*(-k)       |
!
!
! h^{2x2}(k):=
!
! | m-(Cos{kx}+Cos{ky})         & \lambda*(Sin{kx}-i*Sin{ky}) |
! | \lambda*(Sin{kx}+i*Sin{ky}) & -m+(Cos{kx}+Cos{ky})        |
!
! hso^{2x2}(k):=
! | xi*rh*(sin(kx)-xi*sin(ky))  &         \delta              |
! |         -\delta             &             0               |
!
program bhz_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                           :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                     :: Nk,Nktot,Nkpath,Nkx,Npts,L
  integer                                     :: Nky,Nlat,Nx,Ny
  integer                                     :: i,j,k,ik,iorb,jorb,ispin,io,jo
  integer                                     :: ilat,jlat
  integer                                     :: ix,iy,iz
  real(8)                                     :: kx,ky,kz
  real(8),dimension(:,:),allocatable          :: kgrid,kpath,ktrims,Rgrid
  integer,dimension(:,:),allocatable          :: Links
  complex(8),dimension(:,:,:),allocatable     :: Hk
  complex(8),dimension(:,:,:,:),allocatable   :: Hlat
  real(8),dimension(:),allocatable            :: Wtk

  real(8)                                     :: chern,z2
  real(8)                                     :: mh,rh,lambda,delta
  real(8)                                     :: xmu,beta,eps,Eout(2)
  real(8)                                     :: dens(Nso)
  complex(8)                                  :: Hloc(Nso,Nso),arg
  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats,Greal
  complex(8),dimension(:,:,:,:,:,:),allocatable :: iGmats,iGreal
  character(len=20)                           :: file
  logical                                     :: iexist
  complex(8),dimension(Nso,Nso)               :: Gamma1,Gamma2,Gamma5
  complex(8),dimension(:,:,:),allocatable     :: ftHk
  complex(8),dimension(:,:,:,:),allocatable   :: ftHlat
  real(8),dimension(2)                        :: vecK,vecRi,vecRj

  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(L,"L","inputBHZ.conf",default=2048)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=1d0)
  call parse_input_variable(rh,"RH","inputBHZ.conf",default=0d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(delta,"DELTA","inputBHZ.conf",default=0d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=1000.d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz.in")
  call save_input_file("inputBHZ.conf")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(-10d0,"wini")
  call add_ctrl_var(10d0,"wfin")
  call add_ctrl_var(eps,"eps")


  Nky  = Nkx
  Nktot= Nkx*Nky
  !
  Nx   = Nkx
  Ny   = Nkx
  Nlat = Nx*Ny


  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)


  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])


  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:  
  write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  allocate(Wtk(Nktot))
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])
  Wtk = 1d0/Nktot

  call TB_write_hk(Hk,trim(file),Nso,&
       Nd=Norb,Np=1,Nineq=1,&
       Nkvec=[Nkx,Nkx])


  !GET LOCAL PART OF THE HAMILTONIAN
  Hloc=sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc)


  !SOLVE ALONG A PATH IN THE BZ.
  Npts=5
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_X1
  kpath(2,:)=kpoint_Gamma
  kpath(3,:)=kpoint_M1
  kpath(4,:)=kpoint_X1
  kpath(5,:)=kpoint_Gamma
  call TB_Solve_model(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[red1,blue1,red1,blue1],&
       points_name=[character(len=20) :: 'X', 'G', 'M', 'X', 'G'],&
       file="Eigenband.nint")




  !Build the local GF:
  allocate(Gmats(Nspin,Nspin,Norb,Norb,L))
  allocate(Greal(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_matsubara(Hk,Wtk,Gmats,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(Hk,Wtk,Greal,zeros(Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_matsubara(Gmats,"Gloc",1)
  call dmft_print_gf_realaxis(Greal,"Gloc",1)



  !Occupation:
  do ispin=1,Nspin
     do iorb=1,Norb
        io = iorb+(ispin-1)*Norb
        dens(io) = fft_get_density(Gmats(ispin,ispin,iorb,iorb,:),beta)
     enddo
  enddo
  open(10,file="observables.nint")
  write(10,"(20F20.12)")(dens(iorb),iorb=1,Nso),sum(dens)
  close(10)
  write(*,"(A,20F14.9)")"Occupations =",(dens(iorb),iorb=1,Nso),sum(dens)



  !Kinetic Energy
  call dmft_kinetic_energy(Hk,Wtk,zeros(Nspin,Nspin,Norb,Norb,L))



  !GET Z2 INVARIANT:
  allocate(ktrims(2,4))
  ktrims=reshape( [ [0d0,0d0] , [0d0,pi] , [pi,0d0] , [pi,pi] ] , shape(ktrims))
  call get_z2_number(ktrims,[2,4],z2)
  print*,z2


  allocate(kgrid(Nktot,2))
  call TB_build_kgrid([Nkx,Nkx],kgrid)
  z2 = 0d0
  do ik=1,Nktot
     z2 = z2 - chern_nk(kgrid(ik,:))/4d0/pi
  enddo
  print*,z2*(2*pi/Nkx)*(2*pi/Nkx)  
  z2=-simps2d(chern_nk,[-pi,pi],[-pi,pi],N0=200,iterative=.false.)
  print*,z2/4d0/pi
  deallocate(kgrid)

  call get_Chern_Number(Hk,[Nkx,Nkx],2,Nkx/pi2*Nkx/pi2,z2)
  print*,z2

  ! call get_shcond()

  deallocate(Hk)


  !##################################################################


  allocate(Hlat(Nso,Nso,Nlat,Nlat))
  allocate(Hk(Nso,Nso,Nktot))


  !>Build direct Lattice Hamiltonian
  allocate(Links(4,2))          !Links: right,up,left,down
  Links(1,:) = [1,0]
  Links(2,:) = [0,1]
  Links(3,:) = -Links(1,:)
  Links(4,:) = -Links(2,:)
  call TB_build_model(Hlat,ts_model,Nso,[Nx,Ny],Links,pbc=.false.)
  !
  !>Build reciprocal Lattice Hamiltonian
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx],wdos=.false.)


  allocate(iGmats(Nlat,Nspin,Nspin,Norb,Norb,L))
  allocate(iGreal(Nlat,Nspin,Nspin,Norb,Norb,L))
  allocate(ftHk(Nlat*Nso,Nlat*Nso,1))

  do ilat=1,Nlat
     do jlat=1,Nlat
        do io=1,Nso
           do jo=1,Nso
              i = io + (ilat-1)*Nso
              j = jo + (jlat-1)*Nso
              ftHk(i,j,1) = Hlat(io,jo,ilat,jlat)
           enddo
        enddo
     enddo
  enddo
  call dmft_gloc_matsubara(ftHk,[1d0],iGmats,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
  call dmft_gloc_realaxis(ftHk,[1d0],iGreal,zeros(Nlat,Nspin,Nspin,Norb,Norb,L))
  call dmft_print_gf_matsubara(iGmats,"ftGij",4)
  call dmft_print_gf_realaxis(iGreal,"ftGij",4)
  deallocate(iGmats,iGreal,ftHk)
  print*,""
  print*,""
  print*,""



  allocate(Kgrid(Nktot,2))
  call TB_build_kgrid([Nkx,Nkx],Kgrid)
  allocate(Rgrid(Nlat,2))
  call TB_build_Rgrid([Nx,Ny],Rgrid)


  !Build up the K-space Hamiltonian thru FT
  print*,"RECIPROCAL SPACE: H_bhz(kx,ky)=FT[H_bhz(i,j)]"
  allocate(ftHk(Nso,Nso,Nktot))
  allocate(ftHlat(Nso,Nso,Nlat,Nlat))


  !> from Hlat --> Hk = FT(Hlat)
  ftHk=zero
  call start_timer
  do ik=1,Nktot
     vecK = Kgrid(ik,:)
     do ilat=1,Nlat
        vecRi = Rgrid(ilat,:)
        do jlat=1,Nlat
           vecRj = Rgrid(jlat,:)
           !
           arg=dot_product(vecK,vecRj-vecRi)
           !
           ftHk(:,:,ik)=ftHk(:,:,ik) + exp(-xi*arg)*Hlat(:,:,ilat,jlat)/Nlat
           !
        enddo
     enddo
     call eta(ik,Nktot)
  enddo
  where(abs(ftHk)<1.d-6)ftHk=zero
  call stop_timer


  ftHk = ftHk-Hk
  print*,"Identity test for Hk:"
  if(sum(abs(ftHk))>1d-6)then
     print*,"Failed! :-("
  else
     print*,"Passed! :-)"
  endif

  print*,""
  print*,""
  print*,""

  !Build up the real-space Hamiltonian thru FT:"
  print*,"REAL SPACE:       H_bhz(i,j)=FT^-1(H_bhz(kx,ky))"
  ftHlat=zero
  call start_timer
  do ik=1,Nktot
     vecK = Kgrid(ik,:)
     !
     do ilat=1,Nlat
        vecRi = Rgrid(ilat,:)
        do jlat=1,Nlat
           vecRj = Rgrid(jlat,:)
           !
           arg=dot_product(vecK,vecRj-vecRi)
           !
           ftHlat(:,:,ilat,jlat)= ftHlat(:,:,ilat,jlat) + exp(xi*arg)*Hk(:,:,ik)/Nktot!hk_model(vecK,Nso)/Nktot
        enddo
     enddo
     call eta(ik,Nktot)
  enddo
  where(abs(Hlat)<1.d-6)Hlat=zero
  call stop_timer



  ftHlat = ftHlat-Hlat
  print*,"Identity test for Hlat:"
  if(sum(abs(ftHlat))>1d-6)then
     print*,"Failed! :-("
  else
     print*,"Passed! :-)"
  endif



contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8) :: ek
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
  end function hk_model



  function ts_model(link,Nso) result(Hts)
    integer                       :: link
    integer                       :: Nso
    complex(8),dimension(Nso,Nso) :: Hts
    select case(link)
    case (0) !LOCAL PART
       Hts =  Mh*Gamma5
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma1
    case (2) !UP HOPPING
       Hts = -0.5d0*Gamma5 + xi*0.5d0*lambda*Gamma2
    case (3) !LEFT HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma1
    case (4) !DOWN HOPPING
       Hts = -0.5d0*Gamma5 - xi*0.5d0*lambda*Gamma2
    case default 
       stop "ts_model ERROR: link != {0,...,4}"
    end select
  end function ts_model




  ! !+----------------------------------------------------------------+
  ! !PURPOSE  : Get the BHZ lattice Hamiltonian, using Nrow,Ncol
  ! !+----------------------------------------------------------------+
  ! subroutine get_Hij_model(Nx,Ny,Nso,pbc,Hij)
  !   integer                                     :: Nx
  !   integer                                     :: Ny
  !   integer                                     :: Nso
  !   complex(8),dimension(Nso,Nso,Nx*Ny,Nx*Ny)   :: Hij
  !   logical                                     :: pbc
  !   integer                                     :: Nlat,Nlso
  !   integer                                     :: ix,iy
  !   integer                                     :: i,jj,link(4)
  !   !
  !   Nlat=Nx*Ny
  !   Nlso=Nlat*Nso
  !   Hij = 0d0
  !   !
  !   do ix=0,Nx-1
  !      do iy=0,Ny-1
  !         i= 1 + ix + iy*Nx
  !         !
  !         !RIGHT hop 
  !         link(1)= i + 1
  !         if((ix+1)==Nx) then
  !            link(1)=0
  !            if(pbc)link(1)=1+iy*Nx
  !         end if
  !         !UP  hop
  !         link(2)= i + Nx 
  !         if((iy+1)==Ny) then
  !            link(2)=0  
  !            if(pbc)link(2)=ix+1
  !         end if
  !         !LEFT  hop
  !         link(3)= i - 1    
  !         if((ix-1)<0)     then
  !            link(3)=0  
  !            if(pbc)link(3)=Nx+iy*Nx
  !         end if
  !         !DOWN  hop
  !         link(4)= i - Nx 
  !         if((iy-1)<0)     then
  !            link(4)=0  
  !            if(pbc)link(4)=ix+1+(Ny-1)*Nx
  !         end if
  !         !
  !         Hij(:,:,i,i) = Mh*Gamma5
  !         do jj=1,4
  !            if(link(jj)==0)cycle
  !            Hij(:,:,i,link(jj)) = ts_model(jj,Nso)
  !         enddo
  !         !
  !      enddo
  !   enddo
  ! end subroutine get_Hij_model











  subroutine get_z2_number(ktrims,band_indices,z2)
    real(8),dimension(:,:),intent(in)       :: ktrims
    integer,dimension(:),intent(in)         :: band_indices
    complex(8),dimension(:,:,:),allocatable :: Htrims
    real(8),dimension(:,:),allocatable      :: Etrims
    complex(8),dimension(:),allocatable     :: Delta
    real(8)                                 :: z2
    integer                                 :: i,j,Ntrim,itrim,Nocc,unit
    !
    Ntrim=size(Ktrims,2)
    Nocc = size(band_indices)
    allocate(Htrims(Nso,Nso,Ntrim),Etrims(Nocc,Ntrim))
    allocate(Delta(Ntrim))
    !
    do itrim=1,Ntrim
       Htrims(:,:,itrim) = hk_model(Ktrims(:,itrim),Nso)
       Delta(itrim)=-sign(1d0,dreal(Htrims(1,1,itrim)))
    enddo
    z2=product(Delta(:))
    if(z2>0)then
       z2=0d0
    else
       z2=1d0
    end if
    open(free_unit(unit),file="z2_invariant.dat")
    write(unit,*) z2
    close(unit)
  end subroutine get_z2_number


  function nd_model(kpoint,M) result(dk)
    real(8),dimension(:),intent(in) :: kpoint
    integer                         :: M
    real(8),dimension(M)            :: dk
    real(8)                         :: kx,ky,norm
    kx=kpoint(1)
    ky=kpoint(2)
    dk=[lambda*sin(kx),lambda*sin(ky),(mh-cos(kx)-cos(ky))]
    norm = dot_product(dk,dk)
    dk = dk/sqrt(norm)
    where(abs(dk)<1.d-12)dk=0d0
  end function nd_model


  subroutine djac_dk(kpoint,M,ddk)
    real(8),dimension(:)            :: kpoint
    real(8),dimension(size(kpoint)) :: k_
    integer                         :: M
    real(8),dimension(M)            :: fvec,wa1
    real(8)                         :: ddk(M,size(kpoint))
    call djacobian(nd_model,kpoint,M,ddk)
  end subroutine djac_dk

  function chern_nk(kpoint) result(ck)
    real(8),dimension(:) :: kpoint
    real(8) :: dk(3),dk_(3)
    real(8) :: ddk(3,2)
    real(8) :: ck,norm
    dk  = nd_model(kpoint,3)
    call djac_dk(kpoint,3,ddk)
    ck  = s3_product(dk,ddk(:,1),ddk(:,2))
  end function chern_nk




  subroutine get_Chern_Number(Hk,Nkvec,Noccupied,one_over_area,Chern)
    complex(8),intent(in),dimension(:,:,:)    :: Hk    ![Nlso][Nlso][Nktot]
    integer,intent(in),dimension(2)           :: Nkvec ![Nk1][Nk2]: prod(Nkvec)=Nktot
    integer,intent(in)                        :: Noccupied
    real(8),intent(in)                        :: one_over_area
    real(8),intent(out)                       :: Chern
    !
    integer                                   :: Nlso
    integer                                   :: Nktot
    integer                                   :: Nkx,Nky
    integer                                   :: ikx,iky
    integer                                   :: ikxP,ikyP
    integer                                   :: ik,iocc
    complex(8),dimension(:,:),allocatable     :: Eigvec ![Nlso][Nlso]
    real(8),dimension(:),allocatable          :: Eigval ![Nlso]
    complex(8),dimension(:,:),allocatable     :: Gmat
    complex(8),dimension(:,:,:,:),allocatable :: BlochStates ![Nkx][Nky][Noccupied][Nlso]
    complex(8),dimension(4)                   :: Ulink
    real(8),dimension(:,:),allocatable        :: BerryCurvature
    real(8)                                   :: berry_phase
    integer                                   :: unit
    !
    Nlso  = size(Hk,1)
    Nktot = size(Hk,3)
    Nkx   = Nkvec(1)
    Nky   = Nkvec(2)
    call assert_shape(Hk,[Nlso,Nlso,Nktot],"Get_Chern_NUmber","Hk")
    if(Nkx*Nky/=Nktot)stop "ERROR Get_Chern_Number: Nktot = prod(Nkvec)"
    !
    !
    !1. Get the Bloch states from H(:,:,k)
    allocate(Eigvec(Nlso,Nlso))
    allocate(Eigval(Nlso))
    allocate(BlochStates(Nkx,Nky,Noccupied,Nlso))
    allocate(BerryCurvature(Nkx,Nky))
    allocate(Gmat(Noccupied,Noccupied))
    ik=0
    do ikx=1,Nkx
       do iky=1,Nky
          ik=ik+1
          Eigvec = Hk(:,:,ik)
          call eigh(Eigvec,Eigval)
          do iocc=1,Noccupied
             BlochStates(ikx,iky,iocc,:) = Eigvec(:,iocc)
          enddo
       enddo
    enddo
    deallocate(Eigvec,Eigval)
    !
    !
    !2. Evaluate the Berry Curvature
    chern=0d0
    do ikx= 1, Nkx
       ikxP = modulo(ikx,Nkx) + 1
       !ikxM = modulo(ikx-2,Nkx) + 1
       do iky= 1, Nky
          ikyP = modulo(iky,Nky) + 1
          !ikyM = modulo(iky-2,Nky) + 1
          !
          if(Noccupied==1)then
             Ulink(1) = dot_product(BlochStates(ikx,iky,1,:)  , BlochStates(ikx,ikyP,1,:))
             Ulink(2) = dot_product(BlochStates(ikx,ikyP,1,:) , BlochStates(ikxP,ikyP,1,:))
             Ulink(3) = dot_product(BlochStates(ikxP,ikyP,1,:), BlochStates(ikxP,iky,1,:))
             Ulink(4) = dot_product(BlochStates(ikxP,iky,1,:) , BlochStates(ikx,iky,1,:))
             !
          else
             !
             do i=1,Noccupied
                do j=1,Noccupied
                   gmat(i,j) = dot_product(BlochStates(ikx,iky,i,:)  , BlochStates(ikx,ikyP,j,:))
                enddo
             enddo
             Ulink(1) = det(gmat)
             !
             do i=1,Noccupied
                do j=1,Noccupied
                   gmat(i,j) = dot_product(BlochStates(ikx,ikyP,i,:) , BlochStates(ikxP,ikyP,j,:))
                enddo
             enddo
             Ulink(2) = det(gmat)
             !
             do i=1,Noccupied
                do j=1,Noccupied
                   gmat(i,j) = dot_product(BlochStates(ikxP,ikyP,i,:), BlochStates(ikxP,iky,j,:))
                enddo
             enddo
             Ulink(3) = det(gmat)
             !
             do i=1,Noccupied
                do j=1,Noccupied
                   gmat(i,j) = dot_product(BlochStates(ikxP,iky,i,:) , BlochStates(ikx,iky,j,:))
                enddo
             enddo
             Ulink(4) = det(gmat)
             !
          endif
          !
          berry_phase = dimag(zlog( product(Ulink(:))  ))
          chern = chern + berry_phase
          BerryCurvature(ikx,iky) = berry_phase*one_over_area
          !
       enddo
    enddo
    !
    chern=chern/pi2
    !
    open(unit=free_unit(unit),file="Chern_Number.dat")
    write(unit,*)chern
    close(unit)
    !
    call splot3d("Berry_Curvature.dat",&
         linspace(0d0,pi2,Nkx,iend=.false.),&
         linspace(0d0,pi2,Nky,iend=.false.),&
         BerryCurvature(:Nkx,:Nky))
    !
  end subroutine Get_Chern_Number






  subroutine print_Hlat(Hij,file)
    complex(8),dimension(Nso,Nso,Nlat,Nlat) :: Hij
    integer                                 :: ilat,jlat,iso,jso,unit
    character(len=*)                        :: file
    open(free_unit(unit),file=trim(file))
    do ilat=1,Nlat
       do iso=1,Nso
          do jlat=1,Nlat
             do jso=1,Nso
                write(unit,"(F5.2,1x)",advance="no")dreal(Hij(iso,jso,ilat,jlat))
             enddo
             write(unit,"(A2)",advance="no")"  "
          enddo
          write(unit,*)
       enddo
       write(unit,*)""
    enddo
    write(unit,*)
    close(unit)
  end subroutine print_Hlat


  subroutine print_Hk(Hk,file)
    complex(8),dimension(:,:,:)     :: Hk ![Nso][Nso][Nk]
    integer                         :: Nso,Nk,unit
    integer                         :: ik,iso,jso
    character(len=*)                :: file
    open(free_unit(unit),file=trim(file))
    Nk = size(Hk,3)
    Nso= size(Hk,1)
    call assert_shape(Hk,[Nso,Nso,Nk],"print_Hk","Hk")
    do ik=1,Nk
       do iso=1,Nso
          write(unit,"(1000(A1,F6.3,A1,F6.3,A5))")&
               ("(",dreal(Hk(iso,jso,ik)),",",dimag(Hk(iso,jso,ik)),")    ",jso=1,Nso)
       enddo
       write(unit,*)
    enddo
    write(unit,*)
    close(unit)
  end subroutine print_Hk





  ! function get_kinetic_energy(Hk,Liw) result(ed_Ekin)
  !   integer                                  :: Lk,No,Liw
  !   integer                                  :: i,ik,iorb
  !   complex(8),dimension(:,:,:)              :: Hk ![Nso][Nso][Nk]
  !   real(8),dimension(size(Hk,3))            :: Wtk
  !   !
  !   real(8),dimension(:),allocatable         :: wm
  !   complex(8),dimension(:,:),allocatable    :: Ak,Bk
  !   complex(8),dimension(:,:),allocatable    :: Ck,Zk
  !   complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk
  !   real(8)                                  :: Tail0,Tail1,spin_degeneracy
  !   !
  !   real(8)                                  :: H0,ed_Ekin
  !   !
  !   No = size(Hk,1)
  !   Lk = size(Hk,3)
  !   if(No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb_total]"
  !   if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"
  !   !
  !   allocate(wm(Liw))
  !   allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Gk(No,No),Tk(No,No))
  !   !
  !   wm = pi/beta*dble(2*arange(1,Liw)-1)
  !   Wtk=1d0/Lk
  !   !
  !   H0=0d0
  !   Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
  !   do ik=1,Lk
  !      Bk=-Hk(:,:,ik)
  !      do i=1,Liw
  !         Gk = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik)
  !         select case(No)
  !         case default
  !            call matrix_inverse(Gk)
  !         case(1)
  !            Gk = 1d0/Gk
  !         end select
  !         Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
  !         Ck = matmul(Hk(:,:,ik),Gk - Tk)
  !         H0 = H0 + Wtk(ik)*trace_matrix(Ck,No)
  !      enddo
  !   enddo
  !   spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
  !   H0=H0/beta*2.d0*spin_degeneracy
  !   !
  !   Tail0=0d0
  !   Tail1=0d0
  !   do ik=1,Lk
  !      Ak= Hk(:,:,ik)
  !      Bk=-Hk(:,:,ik)
  !      Ck= matmul(Ak,Bk)
  !      Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,No)
  !      Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,No)
  !   enddo
  !   Tail0=spin_degeneracy*Tail0
  !   Tail1=spin_degeneracy*Tail1*beta
  !   ed_Ekin=H0+Tail0+Tail1
  !   deallocate(wm,Ak,Bk,Ck,Zk,Zeta,Gk,Tk)
  ! end function get_kinetic_energy


  ! function get_local_energy(Hk,Liw) result(ed_Eloc)
  !   integer                                  :: Lk,No,Liw
  !   integer                                  :: i,ik,iorb
  !   complex(8),dimension(:,:,:)              :: Hk ![Nso][Nso][Nk]
  !   real(8),dimension(size(Hk,3))            :: Wtk
  !   !
  !   real(8),dimension(:),allocatable         :: wm
  !   complex(8),dimension(:,:),allocatable    :: Ak,Bk,Hloc
  !   complex(8),dimension(:,:),allocatable    :: Ck,Zk
  !   complex(8),dimension(:,:),allocatable    :: Zeta,Gk,Tk
  !   real(8)                                  :: Tail0,Tail1,spin_degeneracy
  !   !
  !   real(8)                                  :: H0,ed_Eloc
  !   !
  !   No = size(Hk,1)
  !   Lk = size(Hk,3)
  !   if(No/=size(Hk,2))stop "get_kinetic_energy: size(Hk,1)!=size(Hk,2) [Norb_total]"
  !   if(Lk/=size(Wtk))stop "get_kinetic_energy: size(Wtk)!=size(Hk,3) [L_k]"
  !   !
  !   allocate(wm(Liw))
  !   allocate(Ak(No,No),Bk(No,No),Ck(No,No),Zk(No,No),Zeta(No,No),Gk(No,No),Tk(No,No),Hloc(No,No))
  !   !
  !   wm = pi/beta*dble(2*arange(1,Liw)-1)
  !   Wtk=1d0/Lk
  !   !
  !   Hloc=sum(Hk,3)/Lk
  !   H0=0d0
  !   Zk=0d0 ; forall(i=1:No)Zk(i,i)=1d0
  !   do ik=1,Lk
  !      Bk=-Hk(:,:,ik)
  !      do i=1,Liw
  !         Gk = (xi*wm(i)+xmu)*Zk(:,:) - Hk(:,:,ik)
  !         select case(No)
  !         case default
  !            call matrix_inverse(Gk)
  !         case(1)
  !            Gk = 1d0/Gk
  !         end select
  !         Tk = Zk(:,:)/(xi*wm(i)) - Bk(:,:)/(xi*wm(i))**2
  !         Ck = matmul(Hloc,Gk - Tk)
  !         H0 = H0 + Wtk(ik)*trace_matrix(Ck,No)
  !      enddo
  !   enddo
  !   spin_degeneracy=3.d0-Nspin !2 if Nspin=1, 1 if Nspin=2
  !   H0=H0/beta*2.d0*spin_degeneracy
  !   !
  !   Tail0=0d0
  !   Tail1=0d0
  !   do ik=1,Lk
  !      Ak= Hloc
  !      Bk=-Hk(:,:,ik)
  !      Ck= matmul(Ak,Bk)
  !      Tail0 = Tail0 + 0.5d0*Wtk(ik)*trace_matrix(Ak,No)
  !      Tail1 = Tail1 + 0.25d0*Wtk(ik)*trace_matrix(Ck,No)
  !   enddo
  !   Tail0=spin_degeneracy*Tail0
  !   Tail1=spin_degeneracy*Tail1*beta
  !   ed_Eloc=H0+Tail0+Tail1
  !   deallocate(wm,Ak,Bk,Ck,Zk,Zeta,Gk,Tk,Hloc)
  ! end function get_local_energy


  ! function trace_matrix(M,dim) result(tr)
  !   integer                       :: dim
  !   complex(8),dimension(dim,dim) :: M
  !   complex(8) :: tr
  !   integer                       :: i
  !   tr=dcmplx(0d0,0d0)
  !   do i=1,dim
  !      tr=tr+M(i,i)
  !   enddo
  ! end function trace_matrix

  ! subroutine get_shcond()
  !   real(8),dimension(L)                :: wm,vm
  !   complex(8),dimension(2,2)           :: g0k,KerU,KerD
  !   complex(8),dimension(2,2,2,Nktot,L) :: Gk
  !   complex(8),dimension(L)             :: Kmats
  !   complex(8),dimension(2,2,2,Nktot)   :: Vkx,Vky
  !   complex(8)                          :: Ksum
  !   real(8)                             :: kx,ky,C_qsh
  !   integer                             :: iw,iv,ik,i,j
  !   wm = pi/beta*real(2*arange(1,L)-1,8)
  !   vm = pi/beta*real(2*arange(1,L)-2,8)
  !   Kmats=zero
  !   ik=0
  !   do i=1,Nkx
  !      kx = kxgrid(i)
  !      do j=1,Nkx
  !         ky = kxgrid(j)
  !         ik=ik+1
  !         Vkx(1,:,:,ik) = sin(kx)*pauli_tau_z + lambda*cos(kx)*pauli_tau_x
  !         Vky(1,:,:,ik) = sin(ky)*pauli_tau_z + lambda*cos(ky)*pauli_tau_y
  !         Vkx(2,:,:,ik) = sin(kx)*pauli_tau_z - lambda*cos(kx)*pauli_tau_x
  !         Vky(2,:,:,ik) = sin(ky)*pauli_tau_z + lambda*cos(ky)*pauli_tau_y
  !         do iw=1,L
  !            g0k = (xi*wm(iw)+xmu)*eye(2)-hk_bhz2x2(kx,ky)
  !            call inv(g0k)
  !            Gk(1,:,:,ik,iw) = g0k
  !            g0k = (xi*wm(iw)+xmu)*eye(2)-conjg(hk_bhz2x2(-kx,-ky))
  !            call inv(g0k)
  !            Gk(2,:,:,ik,iw) = g0k
  !         enddo
  !      enddo
  !   enddo
  !   call start_timer
  !   do iv=1,L
  !      Ksum=zero
  !      do iw=1,L-iv+1
  !         do ik=1,Nktot
  !            KerU = matmul( Vky(1,:,:,ik) , Gk(1,:,:,ik,iw+iv-1) )
  !            KerU = matmul( KerU          , Vkx(1,:,:,ik)      )
  !            KerU = matmul( KerU          , Gk(1,:,:,ik,iw)    )
  !            KerD = matmul( Vky(2,:,:,ik) , Gk(2,:,:,ik,iw+iv-1) )
  !            KerD = matmul( KerD          , Vkx(2,:,:,ik)      )
  !            KerD = matmul( KerD          , Gk(2,:,:,ik,iw)    )
  !            Ksum = Ksum + trace_matrix(KerU,2)-trace_matrix(KerD,2)
  !         enddo
  !      enddo
  !      Kmats(iv) = -Ksum/beta*2*pi/Nktot
  !      call eta(iv,L)
  !   enddo
  !   call stop_timer
  !   C_qsh = dreal(Kmats(2))/vm(2)
  !   open(100,file="qsh_conductance.nint")
  !   write(100,*) C_qsh
  !   close(100)
  !   print*,C_qsh
  ! end subroutine get_shcond


  ! function vkx_model(kpoint,N) result(vkx)
  !   real(8),dimension(:)      :: kpoint
  !   real(8)                   :: kx
  !   complex(8),dimension(N,N) :: vkx
  !   integer                   :: N
  !   kx=kpoint(1)
  !   vkx = zero
  !   vkx(1:2,1:2) = sin(kx)*pauli_tau_z + lambda*cos(kx)*pauli_tau_x
  !   vkx(3:4,3:4) = conjg(sin(-kx)*pauli_tau_z + lambda*cos(-kx)*pauli_tau_x) 
  ! end function vkx_model

  ! function vky_model(kpoint,N) result(vky)
  !   real(8),dimension(:)      :: kpoint
  !   real(8)                   :: ky
  !   complex(8),dimension(N,N) :: vky
  !   integer                   :: N
  !   ky=kpoint(2)
  !   vky = zero
  !   vky(1:2,1:2) = sin(ky)*pauli_tau_z + lambda*cos(ky)*pauli_tau_y
  !   vky(3:4,3:4) = conjg(sin(-ky)*pauli_tau_z + lambda*cos(-ky)*pauli_tau_y) 
  ! end function vky_model

end program bhz_2d


