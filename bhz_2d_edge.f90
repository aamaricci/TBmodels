program bhz_fcc_stripe
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none

  integer,parameter                               :: Norb=2,Nspin=2,Nso=Nspin*Norb
  integer                                         :: Lmats,Lreal
  integer                                         :: Nk,Nktot,Nkpath,Nkx,Ly,Npts,Nslat
  integer                                         :: Nx,Ny
  integer                                         :: i,iy,j,k,ik,ispin,iorb,jorb,ilat,jlat
  integer                                         :: iso,jso,io,jo,ii,jj
  real(8),dimension(2)                            :: bk1
  real(8),dimension(:),allocatable                :: kxgrid
  real(8),dimension(:,:),allocatable              :: kpath,dens,rho,zed
  complex(8),dimension(:,:,:),allocatable         :: Hk,Hkr,Gk
  real(8),dimension(:),allocatable                :: Wtk,wreal,kgrid
  complex(8),dimension(Nso,Nso)                   :: Gamma1,Gamma2,Gamma5
  type(rgb_color),dimension(:,:),allocatable      :: colors

  real(8)                                         :: mh,lambda,e0
  real(8)                                         :: xmu,beta,eps,wini,wfin
  complex(8),dimension(:,:,:,:,:,:),allocatable   :: Gmats,Greal
  complex(8),dimension(:,:,:,:,:,:,:),allocatable :: Gkreal
  real(8),dimension(:,:,:),allocatable            :: Akreal

  character(len=20)                               :: file,nkstring
  logical                                         :: iexist,nogf

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

  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  if(mod(Ly,2)/=0)stop "Error. Please choose use Ly%2=0"


  !SETUP THE GAMMA MATRICES:
  gamma1=kron_pauli( pauli_tau_z, pauli_sigma_x)
  gamma2=kron_pauli( pauli_tau_0,-pauli_sigma_y)
  gamma5=kron_pauli( pauli_tau_0, pauli_sigma_z)


  !START SOLVING THE INHOMO PROBLEM: EDGE STATES IN THE STRIPE GEOMETRY.
  !Generate the Hamiltonian for the full BZ (a set of 1D lines)
  Nslat=Ly*Nso
  allocate(Hkr(Nslat,Nslat,Nkx))
  allocate(Wtk(Nkx))

  bk1 = 2*pi*[1,0]

  print*,"Build H(k,R) model"
  call TB_set_bk(bk1)

  print*,"Solve H(k,R) along -pi:pi"
  Npts=3
  Nktot = (Npts-1)*Nkpath
  allocate(Kpath(Npts,1))  
  kpath(1,:)=[-1]*pi
  kpath(2,:)=[ 0]*pi
  kpath(3,:)=[ 1]*pi
  allocate(colors(Ly,Nso))
  colors = gray88
  colors(1,:) = [red1,gray88,blue1,gray88]
  colors(Ly,:) = [blue1,gray88,red1,gray88]
  call TB_solve_model(bhz_edge_model,Ly,Nso,kpath,Nkpath,&
       colors_name=colors,&
       points_name=[character(len=10) :: "-pi","0","pi"],&
       file="Eigenbands.nint",pbc=.false.)


  call TB_build_model(Hkr,bhz_edge_model,Ly,Nso,Nkvec=[Nkx,1,1],pbc=.false.)
  Wtk = 1d0/Nkx

  if(nogf) stop "Called with NOGF=TRUE! Got bands and exit."

  !Build the local GF:
  allocate(Gmats(Ly,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal(Ly,Nspin,Nspin,Norb,Norb,Lreal))
  Gmats=zero
  Greal=zero
  
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  call dmft_gloc_matsubara(Hkr,Wtk,Gmats,zeros(Ly,Nspin,Nspin,Norb,Norb,Lmats))
  call dmft_print_gf_matsubara(Gmats,"Gloc",4)
  call dmft_gloc_realaxis(Hkr,Wtk,Greal,zeros(Ly,Nspin,Nspin,Norb,Norb,Lreal))
  call dmft_print_gf_realaxis(Greal,"Gloc",4)



  !Build local observables:
  allocate(dens(Ly,Nso),rho(Ly,Nso))
  open(10,file="density.nint")
  open(11,file="rho.nint")
  do iy=1,Ly
     do ispin=1,Nspin
        do iorb=1,Norb
           i = iorb + (ispin-1)*Norb
           dens(iy,i) = fft_get_density(Gmats(iy,ispin,ispin,iorb,iorb,:),beta)
           rho(iy,i)  = -dimag(Gmats(iy,ispin,ispin,iorb,iorb,1))/pi
        enddo
     enddo
     write(10,"(I4,1000F20.12)")iy,(dens(iy,iorb),iorb=1,Nso)
     write(11,"(I4,1000F20.12)")iy,(rho(iy,iorb),iorb=1,Nso)
  enddo
  close(10)
  close(11)

  deallocate(Hkr)
  allocate(Hkr(Nslat,Nslat,Nktot))
  call TB_build_model(Hkr,bhz_edge_model,Ly,Nso,kpath,Nkpath,pbc=.false.)
  
  allocate(Gkreal(Nktot,Ly,Nspin,Nspin,Norb,Norb,Lreal))
  call start_timer
  do ik=1,Nktot
     call dmft_gk_realaxis(Hkr(:,:,ik),1d0,Gkreal(ik,:,:,:,:,:,:),zeros(Ly,Nspin,Nspin,Norb,Norb,Lreal))
     call eta(ik,Nktot)
  enddo
  call stop_timer
  !
  allocate(Akreal(Nktot,Ly,Lreal))
  Akreal = zero
  do ispin=1,Nspin
     do iorb=1,Norb
        Akreal = Akreal -dimag(Gkreal(:,:,ispin,ispin,iorb,iorb,:))/pi/Nspin/Norb
     enddo
  enddo
  !
  allocate(wreal(Lreal))
  allocate(kgrid(Nktot))
  wreal = linspace(wini,wfin,Lreal)
  kgrid = linspace(-pi,pi,Nktot,iend=.false.)
  call splot("GLoc_nso11.dat",wreal,sum(Akreal(:,1,:),1)/Nktot)
  call splot3d("zAkw_nso11.dat",kgrid,wreal,Akreal(:,1,:))
  call splot3d("zAkw_nso22.dat",kgrid,wreal,Akreal(:,2,:))
  call splot3d("zAkw_nso33.dat",kgrid,wreal,Akreal(:,3,:))
  call splot3d("zAkw_nso44.dat",kgrid,wreal,Akreal(:,4,:))
  call splot3d("zAkw_nso55.dat",kgrid,wreal,Akreal(:,5,:))

  
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
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat 
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
