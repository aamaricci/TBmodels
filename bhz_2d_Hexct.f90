! \left(
! \begin{array}{cccc}
!  -\text{Cx}-\text{Cy}+M+\Sigma  & \text{$\Delta $0}+\text{$\Delta $z}+\lambda  \text{Sx}+i \lambda  \text{Sy} & 0 & \text{$\Delta $x} \\
!  \text{$\Delta $0}+\text{$\Delta $z}+\lambda  \text{Sx}-i \lambda  \text{Sy} & \text{Cx}+\text{Cy}-M-\Sigma  & \text{$\Delta $x} & 0 \\
!  0 & \text{$\Delta $x} & -\text{Cx}-\text{Cy}+M+\Sigma  & \text{$\Delta $0}-\text{$\Delta $z}-\lambda  \text{Sx}+i \lambda  \text{Sy} \\
!  \text{$\Delta $x} & 0 & \text{$\Delta $0}-\text{$\Delta $z}-\lambda  \text{Sx}-i \lambda  \text{Sy} & \text{Cx}+\text{Cy}-M-\Sigma  \\
! \end{array}
! \right)


program bhz_2d
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                           :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                     :: Nk,Nktot,Nkpath,Nkx,Npts,Lmats,Lreal
  integer                                     :: Nky,Nlat,Nx,Ny
  integer                                     :: i,j,k,ik,iorb,jorb,ispin,io,jo
  integer                                     :: ilat,jlat
  integer                                     :: ix,iy,iz
  real(8)                                     :: kx,ky,kz
  real(8),dimension(:,:),allocatable          :: kgrid,kpath,ktrims,Rgrid
  complex(8),dimension(:,:,:),allocatable     :: Hk

  real(8)                                     :: chern,z2
  real(8)                                     :: mh,rh,lambda,delta5,delta0,deltaX,deltaZ
  real(8)                                     :: xmu,beta,eps,Eout(2)
  real(8)                                     :: dens(Nso)
  complex(8)                                  :: Hloc(Nso,Nso),arg
  complex(8),dimension(:,:,:,:,:),allocatable :: Gmats,Greal
  complex(8),dimension(:,:,:,:,:,:),allocatable :: iGmats,iGreal
  character(len=20)                           :: infile
  logical                                     :: iexist
  complex(8),dimension(Nso,Nso)               :: Gamma1,Gamma2,Gamma5
  complex(8),dimension(Nso,Nso)               :: Gamma0,GammaX,GammaZ
  complex(8),dimension(:,:,:),allocatable     :: ftHk
  complex(8),dimension(:,:,:,:),allocatable   :: ftHlat
  real(8),dimension(2)                        :: vecK,vecRi,call

  call parse_cmd_variable(infile,"INFILE",default="used.inputED_BHZ.in")
  call parse_input_variable(nkx,"NKX",infile,default=25)
  call parse_input_variable(nkpath,"NKPATH",infile,default=500)
  call parse_input_variable(Lmats,"Lmats",infile,default=2048)
  call parse_input_variable(Lreal,"Lreal",infile,default=2048)
  call parse_input_variable(mh,"MH",infile,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",infile,default=0.3d0)
  call parse_input_variable(delta5,"DELTA5",infile,default=0d0)
  call parse_input_variable(delta0,"DELTA0",infile,default=0d0)
  call parse_input_variable(deltaX,"DELTAX",infile,default=0d0)
  call parse_input_variable(deltaZ,"DELTAZ",infile,default=0d0)
  call parse_input_variable(xmu,"XMU",infile,default=0.d0)
  call parse_input_variable(eps,"EPS",infile,default=4.d-2)
  call parse_input_variable(beta,"BETA",infile,default=1000.d0)
  call save_input_file(infile)
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
  gamma1=kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2=kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z)

  gamma0=kron_pauli( pauli_sigma_0, pauli_tau_x)
  gammaX=kron_pauli( pauli_sigma_x, pauli_tau_x)
  gammaZ=kron_pauli( pauli_sigma_z, pauli_tau_x)

  call TB_set_ei([1d0,0d0],[0d0,1d0])
  call TB_set_bk([pi2,0d0],[0d0,pi2])


  !SOLVE AND PLOT THE FULLY HOMOGENOUS PROBLEM:  
  write(*,*) "Using Nk_total="//txtfy(Nktot)
  allocate(Hk(Nso,Nso,Nktot))
  call TB_build_model(Hk,hk_model,Nso,[Nkx,Nkx])


  !GET LOCAL PART OF THE HAMILTONIAN
  Hloc=sum(Hk,dim=3)/Nktot
  where(abs(Hloc)<1d-6)Hloc=zero
  call TB_write_Hloc(Hloc)


  !SOLVE ALONG A PATH IN THE BZ.
  Npts=4
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_gamma
  kpath(2,:)=kpoint_x1
  kpath(3,:)=kpoint_m1
  kpath(4,:)=kpoint_gamma
  call TB_Solve_model(Hk_model,Nso,kpath,Nkpath,&
       colors_name=[red1,blue1,red3,blue3],&
       points_name=[character(len=20) :: "{/Symbol G}","X","M","{/Symbol G}"],&
       file="Eigenband.nint")
  ! call TB_Solve_model(Hk_diag,Nso,kpath,Nkpath,&
  !      colors_name=[gray0,gray0,gray0,gray0],&
  !      points_name=[character(len=20) :: "{/Symbol G}","X","M","{/Symbol G}"],&
  !      file="DiagEigenband.nint")


  ! !Build the local GF:
  ! allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  ! allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  ! call dmft_gloc_matsubara(Hk,Gmats,zeros(Nspin,Nspin,Norb,Norb,Lmats))
  ! call dmft_gloc_realaxis(Hk,Greal,zeros(Nspin,Nspin,Norb,Norb,Lreal))
  ! call dmft_print_gf_matsubara(Gmats,"G0loc",3)
  ! call dmft_print_gf_realaxis(Greal,"G0loc",3)





  !GET Z2 INVARIANT:
  allocate(ktrims(2,4))
  ktrims=reshape( [ [0d0,0d0] , [0d0,pi] , [pi,0d0] , [pi,pi] ] , shape(ktrims))
  call get_z2_number(ktrims,[2,4],z2)

  deallocate(Hk)


  !##################################################################



contains


  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8) :: ek
    real(8)                   :: kx,ky
    complex(8),dimension(N,N) :: hk,t0,tx,tz,t5
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*Gamma5 + lambda*sin(kx)*Gamma1 + lambda*sin(ky)*Gamma2
    t5 = delta5*Gamma5
    t0 = delta0*Gamma0
    tx = deltaX*GammaX
    tz = deltaZ*GammaZ
    Hk = Hk + t5 + tx + t0  + tz
  end function hk_model


  function hk_diag(kpoint,N) result(hk)
    real(8),dimension(:)      :: kpoint
    integer                   :: N
    real(8)                   :: ek,cx,cy,sx,sy,sk,dd,d1,sq,ep,em
    real(8)                   :: kx,ky
    real(8),dimension(N)      :: hkdiag
    complex(8),dimension(N,N) :: hk,t0,tx,tz,t5
    if(N/=4)stop "hk_model: error in N dimensions"
    kx=kpoint(1)
    ky=kpoint(2)
    cx = cos(kx)
    cy = cos(ky)
    sx = sin(kx)
    sy = sin(ky)
    ek = (cx + cy - Mh - delta5)**2
    sk = sx**2 + sy**2
    dd = delta0**2 + deltaX**2 + deltaZ**2
    d1 = deltaX**2 + (DeltaZ+Sx*lambda)**2
    sq = sqrt(d1*delta0**2)
    em = sqrt(dd + 2*Sx*deltaZ*lambda + sk*lambda**2 + ek - 2*sq)
    ep = sqrt(dd + 2*Sx*deltaZ*lambda + sk*lambda**2 + ek + 2*sq)
    Hk = diag([-em,em,-ep,ep])
  end function hk_diag













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

end program bhz_2d


