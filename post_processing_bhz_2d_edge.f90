program bhz_fcc_stripe
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                         :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                   :: Lmats,Lreal
  integer                                   :: Nk,Nktot,Nkpath,Nkx,Ly,Npts,Nslat
  integer                                   :: Nx,Ny
  integer                                   :: i,j,k,ik,iorb,jorb,ilat,jlat
  integer                                   :: iso,jso,io,jo,ii,jj
  integer,allocatable                       :: ixRvec(:),iyRvec(:)
  integer                                   :: ix,iy,iz
  real(8)                                   :: kx,ky,kz
  real(8),dimension(:),allocatable          :: kxgrid
  real(8),dimension(:,:),allocatable        :: kpath,dens,rho,zed
  complex(8),dimension(:,:,:),allocatable   :: Hk,Hkr,Gk
  real(8),dimension(:),allocatable          :: wm,wr
  complex(8),dimension(Nso,Nso)             :: Gamma1,Gamma2,Gamma5

  real(8)                                   :: mh,lambda,e0
  real(8)                                   :: xmu,beta,eps,wini,wfin
  complex(8),dimension(:,:,:,:),allocatable :: Gmats,Greal,Hlat
  real(8),dimension(:,:,:,:,:),allocatable  :: Akreal
  character(len=20)                         :: file,nkstring
  logical                                   :: iexist,nogf
  complex(8),allocatable,dimension(:,:,:,:) :: Zeta_mats,Zeta_real,Smats,Sreal
  real(8),dimension(:,:),allocatable        :: Zmats

  call parse_input_variable(nkx,"NKX","inputBHZ.conf",default=25)
  call parse_input_variable(Ly,"Ly","inputBHZ.conf",default=50)
  call parse_input_variable(nkpath,"NKPATH","inputBHZ.conf",default=500)
  call parse_input_variable(e0,"e0","inputBHZ.conf",default=1d0)
  call parse_input_variable(mh,"MH","inputBHZ.conf",default=3d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.conf",default=0.3d0)
  call parse_input_variable(Lmats,"Lmats","inputBHZ.conf",default=4192)
  call parse_input_variable(Lreal,"Lreal","inputBHZ.conf",default=1024)
  call parse_input_variable(xmu,"XMU","inputBHZ.conf",default=0.d0)
  call parse_input_variable(wini,"wini","inputBHZ.conf",default=-12d0)
  call parse_input_variable(wfin,"wfin","inputBHZ.conf",default=12d0)
  call parse_input_variable(eps,"EPS","inputBHZ.conf",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.conf",default=100.d0)
  call parse_input_variable(file,"FILE","inputBHZ.conf",default="hkfile_bhz.in")
  call parse_input_variable(nogf,"NOGF","inputBHZ.conf",default=.false.)
  call save_input_file("inputBHZ.conf")

  if(mod(Ly,2)/=0)stop "Error. Please choose use Ly%2=0"


  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)


  !START SOLVING THE INHOMO PROBLEM: EDGE STATES IN THE STRIPE GEOMETRY.
  !Generate the Hamiltonian for the full BZ (a set of 1D lines)
  Nslat=Ly*Nso
  allocate(kxgrid(Nkx))
  allocate(Hkr(Nslat,Nslat,Nkx))

  allocate(Smats(Ly,Nso,Nso,Lmats))
  allocate(Sreal(Ly,Nso,Nso,Lreal))
  allocate(Zmats(Nslat,Nslat))
  Sreal=zero
  Smats=zero
  Zmats=eye(Nslat)


  print*,"Build H(k,R) model"
  kxgrid = kgrid(Nkx)
  Hkr = TB_build_model(bhz_edge_model,Ly,Nso,kxgrid,[0d0],[0d0],pbc=.false.)
  call write_hk_w90("Hkrfile_BHZ.data",&
       No=Nslat,&
       Nd=Norb,&
       Np=0,&
       Nineq=Ly,&
       Hk=Hkr,&
       kxgrid=kxgrid,kygrid=[0d0],kzgrid=[0d0])


  print*,"Solve H(k,R) along -pi:pi"
  Npts=3
  allocate(Kpath(Npts,1))
  kpath(1,:)=[-1]*pi
  kpath(2,:)=[ 0]*pi
  kpath(3,:)=[ 1]*pi
  call TB_solve_path(bhz_edge_model,Ly,Nso,kpath,Nkpath,&
       colors_name=[red1,gray88,blue1,gray88,blue1,gray88,red1,gray88],&
       points_name=[character(len=10) :: "-pi","0","pi"],&
       file="Eigenbands.nint",pbc=.false.)

  if(nogf) stop "Called with NOGF=TRUE! Got bands and exit."



  !Build the local GF:
  call read_sigma_matsubara(Smats,Ly/2)
  call read_sigma_real(Sreal,Ly/2)


  allocate(wm(Lmats),wr(Lreal))
  allocate(Gmats(Ly,Nso,Nso,Lmats),Zeta_mats(Ly,Nso,Nso,Lmats))
  allocate(Greal(Ly,Nso,Nso,Lreal),Zeta_real(Ly,Nso,Nso,Lreal))
  allocate(Gk(Ly,Nso,Nso))
  allocate(Akreal(Ly,Nso,Nso,Nkx,Lreal))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  wr = linspace(wini,wfin,Lreal)
  do ilat=1,Ly
     do i=1,Lmats
        zeta_mats(ilat,:,:,i) = (xi*wm(i)+xmu)*eye(Nso)     - Smats(ilat,:,:,i)
     enddo
     do i=1,Lreal
        zeta_real(ilat,:,:,i) = (wr(i)+xi*eps+xmu)*eye(Nso) - Sreal(ilat,:,:,i)
     enddo
  enddo
  call start_timer()
  do ik=1,Nkx
     do i=1,Lmats
        Gk = inverse_g0k(Zeta_mats(:,:,:,i),Ly,Nso,Hkr(:,:,ik))
        Gmats(:,:,:,i)=Gmats(:,:,:,i) + Gk
     enddo
     do i=1,Lreal
        Gk = inverse_g0k(Zeta_real(:,:,:,i),Ly,Nso,Hkr(:,:,ik))
        Greal(:,:,:,i)=Greal(:,:,:,i) + Gk
        Akreal(:,:,:,ik,i)=-dimag(Gk)/pi
     enddo
     call eta(ik,Nkx)
  enddo
  call stop_timer()
  Gmats=Gmats/Nkx
  Greal=Greal/Nkx

  call splot3d("Akw_11.nint",kxgrid,wr,Akreal(1,1,1,:,:))
  call splot3d("Akw_22.nint",kxgrid,wr,Akreal(2,1,1,:,:))
  call splot3d("Akw_33.nint",kxgrid,wr,Akreal(3,1,1,:,:))
  call splot3d("Akw_44.nint",kxgrid,wr,Akreal(4,1,1,:,:))
  call splot3d("Akw_55.nint",kxgrid,wr,Akreal(5,1,1,:,:))

  !Build local observables:
  allocate(dens(Ly,Nso),rho(Ly,Nso),zed(Ly,Nso))
  open(10,file="density.nint")
  open(11,file="rho.nint")
  open(12,file="zeta.nint")
  do iy=1,Ly
     do iorb=1,Nso
        call splot("Gloc_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_site"//reg(txtfy(iy,4))//"_iw.nint",wm,Gmats(iy,iorb,iorb,:))
        call splot("Gloc_l"//reg(txtfy(iorb))//reg(txtfy(iorb))//"_site"//reg(txtfy(iy,4))//"_realw.nint",wr,&
             -dimag(Greal(iy,iorb,iorb,:))/pi,dreal(Greal(iy,iorb,iorb,:)))
        dens(iy,iorb) = fft_get_density(Gmats(iy,iorb,iorb,:),beta)
        rho(iy,iorb)  = -dimag(Gmats(iy,iorb,iorb,1))/pi
        zed(iy,iorb)  = 1.d0/( 1.d0 + abs( dimag(Smats(iy,iorb,iorb,1))/(pi/beta) ))
     enddo
     write(10,"(I4,1000F20.12)")iy,(dens(iy,iorb),iorb=1,Nso)
     write(11,"(I4,1000F20.12)")iy,(rho(iy,iorb),iorb=1,Nso)
     write(12,"(I4,1000F20.12)")iy,(zed(iy,iorb),iorb=1,Nso)
  enddo
  close(10)
  close(11)
  close(12)



  print*,"Solve z*H(k,R) along -pi:pi"
  do ilat=1,Ly
     do iorb=1,Nso
        i = iorb + (ilat-1)*Nso
        Zmats(i,i) = zed(ilat,iorb)
     enddo
  enddo
  call TB_solve_path(bhz_edge_model,Ly,Nso,kpath,Nkpath,&
       colors_name=[red1,gray88,blue1,gray88,blue1,gray88,red1,gray88],&
       points_name=[character(len=10) :: "-pi","0","pi"],&
       file="Ren_Eigenbands.nint",pbc=.false.)




contains


  !BHZ on a stripe geometry;
  function bhz_edge_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx
    integer                             :: Nlat,N
    complex(8),dimension(N,N)           :: Hmat,Tmat,TmatH
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    kx=kpoint(1)
    Hrk=zero
    Hmat=h0_rk_bhz(kx,N)
    Tmat=t0_rk_bhz(N)
    TmatH=conjg(transpose(Tmat))
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat + dreal(Smats(i,:,:,1)) !< H(k) + Re(Sigma_iy(:Nso,:Nso;omega=0))
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
    Hrk = matmul(Zmats,Hrk)
  end function bhz_edge_model

  function h0_rk_bhz(kx,N) result(H)
    real(8)                    :: kx
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    H = (mh-e0*cos(kx))*gamma5 + lambda*sin(kx)*gamma1
  end function h0_rk_bhz

  function t0_rk_bhz(N) result(H)
    integer                    :: N
    complex(8),dimension(N,N)  :: H
    H = -0.5d0*e0*gamma5 + xi*0.5d0*lambda*gamma2
  end function T0_rk_bhz





  function inverse_g0k(zeta,Nlat,N,hk) result(g0k)
    complex(8),dimension(Nlat,N,N)               :: zeta
    integer                                      :: Nlat,N,i,j
    complex(8),dimension(Nlat*N,Nlat*N)          :: Hk
    complex(8),dimension(Nlat,N,N)               :: Diag
    complex(8),dimension(Nlat-1,N,N)             :: Sub,Over
    complex(8),dimension(Nlat,N,N)               :: g0k
    call get_tridiag(Nlat,N,Hk,Sub,Diag,Over)
    Diag = zeta - Diag
    call inv_tridiag(Nlat,N,-Sub,Diag,-Over,G0k)
  end function inverse_g0k







  !+----------------------------------------------------------------+
  !PURPOSE  : Get the BHZ lattice Hamiltonian, using Nrow,Ncol
  !+----------------------------------------------------------------+
  function get_lattice_hamiltonian(Nx,Ny,Nso,xpbc,ypbc) result(Hij)
    integer                                     :: Nx
    integer                                     :: Ny
    integer                                     :: Nso
    complex(8),dimension(Nso,Nso,Nx*Ny,Nx*Ny)   :: Hij
    logical,optional                            :: xpbc,ypbc
    logical                                     :: xpbc_,ypbc_
    integer                                     :: Nlat,Nlso
    integer                                     :: ix,iy
    integer                                     :: i,jj,link(4)
    !
    xpbc_=.true.;if(present(xpbc))xpbc_=xpbc
    ypbc_=.true.;if(present(ypbc))ypbc_=ypbc
    !
    Nlat=Nx*Ny
    Nlso=Nlat*Nso
    Hij = 0d0
    !
    !+- 2D LATTICE (Nx x Ny) -+!
    allocate(ixRvec(Nlat),iyRvec(Nlat))
    do ix=0,Nx-1
       do iy=0,Ny-1
          i= 1 + ix + iy*Nx
          !
          ixRvec(i)=ix+1
          iyRvec(i)=iy+1
          !
          !RIGHT hop 
          link(1)= i + 1
          if((ix+1)==Nx) then
             link(1)=0
             if(xpbc_)link(1)=1+iy*Nx
          end if
          !UP  hop
          link(2)= i + Nx 
          if((iy+1)==Ny) then
             link(2)=0  
             if(ypbc_)link(2)=ix+1
          end if
          !LEFT  hop
          link(3)= i - 1    
          if((ix-1)<0)     then
             link(3)=0  
             if(xpbc_)link(3)=Nx+iy*Nx
          end if
          !DOWN  hop
          link(4)= i - Nx 
          if((iy-1)<0)     then
             link(4)=0  
             if(ypbc_)link(4)=ix+1+(Ny-1)*Nx
          end if
          !
          Hij(:,:,i,i) = Mh*Gamma5
          do jj=1,4
             if(link(jj)==0)cycle
             Hij(:,:,i,link(jj)) = ts_model(jj)
          enddo
          !
       enddo
    enddo
  end function get_lattice_hamiltonian

  function ts_model(link) result(Hts)
    integer                       :: link
    complex(8),dimension(Nso,Nso) :: Hts
    select case(link)
    case default 
       stop "ERROR ts_model: link != 4" 
    case (1) !RIGHT HOPPING
       Hts = -0.5d0*e0*Gamma5 + xi*0.5d0*lambda*Gamma1
    case (2) !UP HOPPING
       Hts = -0.5d0*e0*Gamma5 + xi*0.5d0*lambda*Gamma2
    case (3) !LEFT HOPPING
       Hts = -0.5d0*e0*Gamma5 - xi*0.5d0*lambda*Gamma1
    case (4) !DOWN HOPPING
       Hts = -0.5d0*e0*Gamma5 - xi*0.5d0*lambda*Gamma2
    end select
  end function ts_model







  !----------------------------------------------------------------------------------------!
  ! purpose: read the local self-energy from disk
  !----------------------------------------------------------------------------------------!
  subroutine read_sigma_matsubara(Self,Nineq)
    complex(8),dimension(:,:,:,:),intent(inout)     :: Self
    integer,optional                                :: Nineq
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Self_ineq
    character(len=30)                               :: suffix
    integer                                         :: ilat,ispin,iorb,ineq
    real(8),dimension(:),allocatable                :: ome
    call assert_shape(Self,[Ly,Nspin*Norb,Nspin*Norb,Lmats],"read_sigma_matsubara","Self")
    allocate(ome(Lmats))
    if(present(Nineq))then
       allocate(Self_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
    else
       allocate(Self_ineq(Ly,Nspin,Nspin,Norb,Norb,Lmats))
    endif
    Self_ineq=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_iw.ed"
          call read_data("LSigma"//trim(suffix),Self_ineq(:,ispin,ispin,iorb,iorb,:),ome)
       enddo
    enddo
    if(present(Nineq))then
       do ilat=1,Ly
          do i=1,Lmats
             Self(ilat,:,:,i) = nn2nso_reshape( Self_ineq(ilat2ineq(ilat,Nineq),:,:,:,:,i) , Nspin,Norb)
          enddo
       enddo
    else
       do ilat=1,Ly
          do i=1,Lmats
             Self(ilat,:,:,i) = nn2nso_reshape( Self_ineq(ilat,:,:,:,:,i)                  , Nspin,Norb)
          enddo
       enddo
    endif
  end subroutine read_sigma_matsubara

  subroutine read_sigma_real(Self,Nineq)
    complex(8),dimension(:,:,:,:),intent(inout)     :: Self
    integer,optional                                :: Nineq
    complex(8),dimension(:,:,:,:,:,:),allocatable   :: Self_ineq
    character(len=30)                               :: suffix
    integer                                         :: ilat,ispin,iorb,ineq
    real(8),dimension(:),allocatable                :: ome
    call assert_shape(Self,[Ly,Nspin*Norb,Nspin*Norb,Lreal],"read_sigma_real","Self")
    allocate(ome(Lreal))
    if(present(Nineq))then
       allocate(Self_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
    else
       allocate(Self_ineq(Ly,Nspin,Nspin,Norb,Norb,Lreal))
    endif
    Self_ineq=zero
    do ispin=1,Nspin
       do iorb=1,Norb
          suffix="_l"//reg(txtfy(iorb))//"_s"//reg(txtfy(ispin))//"_realw.ed"
          call read_data("LSigma"//trim(suffix),Self_ineq(:,ispin,ispin,iorb,iorb,:),ome)
       enddo
    enddo
    do ilat=1,Ly
       ineq=ilat
       if(present(Nineq))ineq=ilat2ineq(ilat,Nineq)
       do i=1,Lreal
          Self(ilat,:,:,i) = nn2nso_reshape( Self_ineq(ineq,:,:,:,:,i) , Nspin,Norb)
       enddo
    enddo
  end subroutine read_sigma_real






  function ilat2ineq(ilat,Nineq) result(ineq)
    integer,intent(in) :: ilat,Nineq
    integer            :: ineq
    ineq=ilat
    if( ilat > Nineq )ineq=Ly-ilat+1
  end function ilat2ineq
  function ineq2ilat(ineq,Nineq) result(ilat)
    integer,intent(in) :: ineq,Nineq
    integer            :: ilat
    ilat=ineq
    if(ineq>Nineq)stop "ineq2ilat error: called with ineq > Nineq"
  end function ineq2ilat




  function nn2nso_reshape(Hnn,Nspin,Norb) result(Hso)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnn
    integer                                     :: Nspin,Norb
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hso
    integer                                     :: iorb,ispin,is
    integer                                     :: jorb,jspin,js
    Hso=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb  !spin-orbit stride
                js = jorb + (jspin-1)*Norb  !spin-orbit stride
                Hso(is,js) = Hnn(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function nn2nso_reshape



end program bhz_fcc_stripe















! subroutine solve_along_path(hkr_model,Nlat,Nso,kpath,Nkpath,colors_name,points_name,file,pbc)
!   interface 
!      function hkr_model(kpoint,Nlat,Nso,pbc)
!        real(8),dimension(:)                    :: kpoint
!        integer                                 :: Nlat,Nso
!        logical                                 :: pbc
!        complex(8),dimension(Nlat*Nso,Nlat*Nso) :: hkr_model
!      end function hkr_model
!   end interface
!   integer                                      :: Nlat,Nso,Nlso
!   real(8),dimension(:,:)                       :: kpath
!   integer                                      :: Nkpath,Nktot
!   type(rgb_color),dimension(2*Nso)             :: colors_name
!   character(len=*),dimension(size(kpath,1))    :: points_name
!   character(len=*),optional                    :: file
!   logical                                      :: pbc
!   character(len=256)                           :: file_,xtics
!   integer                                      :: Npts,units(Nlat*Nso)
!   integer                                      :: ipts,ik,ic,unit,iorb,ilat,io,nrot,u1,u2
!   real(8)                                      :: coeff(Nlat*Nso)
!   type(rgb_color)                              :: corb(Nlat*Nso),c(Nlat*Nso)
!   real(8),dimension(size(kpath,2))             :: kstart,kstop,kpoint,kdiff
!   complex(8),dimension(Nlat*Nso,Nlat*Nso)      :: h
!   real(8),dimension(Nlat*Nso)                  :: Eval
!   character(len=64)                            :: fmt
!   !
!   Nlso  = Nlat*Nso
!   write(fmt,"(A,I0,A)")"(I12,",Nlso,"(F18.12,I18))"
!   file_ = "Eigenbands.tb";if(present(file))file_=file
!   Npts  = size(kpath,1)
!   Nktot = (Npts-1)*Nkpath
!   !
!   do io=1,Nso
!      corb(io) = colors_name(io)
!   enddo
!   do io=Nso+1,(Nlat-1)*Nso
!      corb(io) = gray88
!   enddo
!   do io=1,Nso
!      corb( (Nlat-1)*Nso + io) = colors_name(Nso + io)
!   enddo
!   !
!   units=free_units(Nlso)
!   do io=1,Nlso
!      open(units(io),file="site_"//reg(txtfy(io,4))//"_"//reg(file_))
!   enddo
!   open(free_unit(unit),file=reg(file_))
!   ic=0
!   do ipts=1,Npts-1
!      kstart = kpath(ipts,:)
!      kstop  = kpath(ipts+1,:)
!      kdiff  = (kstop-kstart)/Nkpath
!      do ik=1,Nkpath
!         ic=ic+1
!         kpoint = kstart + (ik-1)*kdiff
!         h = hkr_model(kpoint,Nlat,Nso,pbc)
!         call eigh(h,Eval)
!         do io=1,Nlso
!            coeff(:)=h(:,io)*conjg(h(:,io))
!            c(io) = coeff.dot.corb
!            write(units(io),"(I12,F18.12,I18)")ic,Eval(io),rgb(c(io))
!         enddo
!         write(unit,fmt)ic,(Eval(io),rgb(c(io)),io=1,Nlso)
!      enddo
!   enddo
!   close(unit)
!   do io=1,Nlso
!      close(units(io))
!   enddo
!   xtics="'"//reg(points_name(1))//"'1,"
!   do ipts=2,Npts-1
!      xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//reg(txtfy((ipts-1)*Nkpath+1))//","
!   enddo
!   xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//reg(txtfy((Npts-1)*Nkpath))//""
!   open(unit,file=reg(file_)//".gp")
!   write(unit,*)"set term wxt"
!   write(unit,*)"unset key"
!   write(unit,*)"set xtics ("//reg(xtics)//")"
!   write(unit,*)"set grid ytics xtics"
!   io=1
!   write(unit,*)"plot 'site_"//reg(txtfy(io,4))//"_"//reg(file_)//"' u 1:2:3 w l lw 2 lc rgb variable"
!   do io=2,Nlso
!      write(unit,*)"rep 'site_"//reg(txtfy(io,4))//"_"//reg(file_)//"' u 1:2:3 w l lw 2 lc rgb variable"
!   enddo
!   write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
!   write(unit,*)"#set out '"//reg(file_)//".png'"
!   write(unit,*)"#rep"
!   write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
!   write(unit,*)"#set out '"//reg(file_)//".svg'"
!   write(unit,*)"#rep"
!   write(unit,*)"set term postscript eps enhanced color 'Times'"
!   write(unit,*)"set output '|ps2pdf - "//reg(file_)//".pdf'"
!   write(unit,*)"rep"
!   close(unit)
!   call system("chmod +x "//reg(file_)//".gp")
! end subroutine solve_along_path
