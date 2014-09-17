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

program tdd_pam
  USE CONSTANTS
  USE PARSE_INPUT
  USE TOOLS
  USE ARRAYS
  USE IOTOOLS
  USE FUNCTIONS
  USE TIMER
  USE MATRIX
  implicit none

  integer,parameter    :: L=1000,Norb=4
  integer              :: i,j,k,ik,iorb,jorb
  real(8)              :: mh,rh,lambda,delta
  real(8)              :: xmu,beta,eps
  integer              :: dcshift,count
  real(8)              :: epsik,fmesh,n(Norb),eig(Norb)
  integer              :: Nkx,Nk
  real(8)              :: ix,iy
  real(8)              :: kx,ky
  real(8),dimension(L) :: wm,wr
  complex(8)           :: w,Hk(Norb,Norb),Hloc(Norb,Norb)
  complex(8)           :: fg(L,Norb,Norb),fgr(L,Norb,Norb)
  character(len=20)    :: file,nkstring
  logical              :: iexist,ibool,dcflag
  real,dimension(3)    :: color
  real,dimension(Norb,3) :: color_bands
  real,dimension(3),parameter :: &
       red=[1.,0.,0.],&
       green=[0.,1.,0.],&
       blue=[0.,0.,1.],&
       black=[0.,0.,0.],&
       bc=[0.5,1.,0.]


  call parse_input_variable(nkx,"NKX","inputBHZ.in",default=25)
  call parse_input_variable(mh,"MH","inputBHZ.in",default=3.d0)
  call parse_input_variable(rh,"RH","inputBHZ.in",default=0.d0)
  call parse_input_variable(lambda,"LAMBDA","inputBHZ.in",default=0.3d0)
  call parse_input_variable(delta,"DELTA","inputBHZ.in",default=0.d0)
  call parse_input_variable(xmu,"XMU","inputBHZ.in",default=0.d0)
  call parse_input_variable(eps,"EPS","inputBHZ.in",default=4.d-2)
  call parse_input_variable(beta,"BETA","inputBHZ.in",default=1000.d0)
  call parse_input_variable(file,"FILE","inputBHZ.in",default="hkfile_bhz.in")


  wm = pi/beta*real(2*arange(1,L)-1,8)

  wr = linspace(-10.d0,10.d0,L,mesh=fmesh)

  Nk=Nkx*Nkx
  write(*,*) "Using Nk="//txtfy(Nk)
  open(50,file=trim(file))
  write(nkstring,*)Nk
  write(50,*)trim(adjustl(trim(Nkstring)))," 2 1 1 1"

  fgr=zero ; fg =zero ;count=0; Hloc=zero
  call start_progress
  do ix=1,Nkx
     kx = -pi + 2.d0*pi*dble(ix-1)/dble(Nkx)
     do iy=1,Nkx
        ky = -pi + 2.d0*pi*dble(iy-1)/dble(Nkx)
        Hk(:,:) = hk_bhz(kx,ky)
        write(50,"(3(F10.7,1x))")kx,ky,pi
        do iorb=1,Norb
           write(50,"(10(2F10.7,1x))")(Hk(iorb,jorb),jorb=1,Norb)
        enddo
        Hloc=Hloc+Hk/dble(Nk)
        do i=1,L
           w = dcmplx(wr(i),eps)+xmu
           fgr(i,:,:)=fgr(i,:,:) + inverse_g0k(w,Hk)
           w = xi*wm(i)+xmu
           fg(i,:,:) =fg(i,:,:)  + inverse_g0k(w,Hk)
        enddo
        count=count+1
        call progress_bar(count,Nk)
     enddo
  enddo
  call stop_progress
  fgr= fgr/real(Nk,8)
  fg = fg/real(Nk,8)

  open(10,file="bhz_hloc.dat")
  do iorb=1,Norb
     write(10,"(90F21.12)")(dreal(Hloc(iorb,jorb)),jorb=1,Norb)
  enddo
  write(10,*)""
  do iorb=1,Norb
     write(10,"(90F21.12)")(dimag(Hloc(iorb,jorb)),jorb=1,Norb)
  enddo
  close(10)

  ik = 0
  color_bands(1,:)=black
  color_bands(2,:)=[1.0,1.0,0.0]
  color_bands(3,:)=red
  color_bands(4,:)=blue

  open(11,file="Eigenband_1.dat")
  open(12,file="Eigenband_2.dat")
  open(13,file="Eigenband_3.dat")
  open(14,file="Eigenband_4.dat")
  !From \Gamma=(0,0) to X=(pi,0): 100 steps
  do ix=1,Nk
     ik=ik+1
     kx = 0.d0 + pi*real(ix-1,8)/dble(Nk)
     ky = 0.d0
     Hk(:,:) = hk_bhz(kx,ky)
     eig = Eigk(Hk)
     color(1) = dot_product( dreal(Hk(:,1))*color_bands(:,1) , dreal(Hk(:,1))*color_bands(:,1) )
     color(2) = dot_product( dreal(Hk(:,1))*color_bands(:,2) , dreal(Hk(:,1))*color_bands(:,2) )
     color(3) = dot_product( dreal(Hk(:,1))*color_bands(:,3) , dreal(Hk(:,1))*color_bands(:,3) )
     write(10,"(I,20F25.12)")ik,eig(1),(color(i),i=1,3)

     write(11,"(I,3F25.12)")ik,eig(1),sum((Hk(:,1)))
     write(12,"(I,3F25.12)")ik,-eig(1),-sum((Hk(:,1)))
     write(13,"(I,3F25.12)")ik,eig(3),sum((Hk(:,3)))
     write(14,"(I,3F25.12)")ik,-eig(4),-sum((Hk(:,3)))
  enddo
  !From X=(pi,0) to M=(pi,pi): 100 steps
  do iy=1,Nk
     ik=ik+1
     kx = pi
     ky = 0.d0 + pi*real(iy-1,8)/dble(Nk)
     Hk(:,:) = hk_bhz(kx,ky)
     eig = Eigk(Hk)
     color(1) = dot_product( dreal(Hk(:,1))*color_bands(:,1) , dreal(Hk(:,1))*color_bands(:,1) )
     color(2) = dot_product( dreal(Hk(:,1))*color_bands(:,2) , dreal(Hk(:,1))*color_bands(:,2) )
     color(3) = dot_product( dreal(Hk(:,1))*color_bands(:,3) , dreal(Hk(:,1))*color_bands(:,3) )
     write(10,"(I,20F25.12)")ik,eig(1),(color(i),i=1,3)
     write(11,"(I,3F25.12)")ik,eig(1),sum((Hk(:,1)))
     write(12,"(I,3F25.12)")ik,-eig(1),-sum((Hk(:,1)))
     write(13,"(I,3F25.12)")ik,eig(3),sum((Hk(:,3)))
     write(14,"(I,3F25.12)")ik,-eig(4),-sum((Hk(:,3)))
  enddo
  !From M=(pi,pi) to \Gamma=(0,0): 100 steps
  do ix=1,Nk
     ik=ik+1
     iy=ix
     kx = pi - pi*real(ix-1,8)/dble(Nk)
     ky = pi - pi*real(iy-1,8)/dble(Nk)
     Hk(:,:) = hk_bhz(kx,ky)
     eig = Eigk(Hk)
     color(1) = dot_product( dreal(Hk(:,1))*color_bands(:,1) , dreal(Hk(:,1))*color_bands(:,1) )
     color(2) = dot_product( dreal(Hk(:,1))*color_bands(:,2) , dreal(Hk(:,1))*color_bands(:,2) )
     color(3) = dot_product( dreal(Hk(:,1))*color_bands(:,3) , dreal(Hk(:,1))*color_bands(:,3) )
     write(10,"(I,20F25.12)")ik,eig(1),(color(i),i=1,3)
     write(12,"(I,3F25.12)")ik,-eig(1),-sum((Hk(:,1)))
     write(13,"(I,3F25.12)")ik,eig(3),sum((Hk(:,3)))
     write(14,"(I,3F25.12)")ik,-eig(4),-sum((Hk(:,3)))
  enddo
  do i=1,4
     close(10+i)
  enddo

  call splot("DOS.dat",wr,-dimag(fgr(:,1,1))/pi,-dimag(fgr(:,2,2))/pi,-dimag(fgr(:,3,3))/pi,-dimag(fgr(:,4,4))/pi)
  call splot("G_iw.dat",wm,fg(:,1,1),fg(:,2,2),fg(:,3,3),fg(:,4,4))

  do iorb=1,Norb
     n(iorb) = -2.d0*sum(dimag(fgr(:,iorb,iorb))*fermi(wr(:),beta))*fmesh/pi
  enddo
  open(10,file="observables.dat")
  write(10,"(14F20.12)")mh,lambda,rh,delta,xmu,(n(iorb),iorb=1,Norb),sum(n)
  close(10)
  write(*,"(A,4F14.9)")"Occupations                   =",n


  print*,"Remember to open the file:"//trim(file)
  print*,"Parameters for "//trim(file)//" are in +paramaters4_"//trim(file)
  ! open(10,file="parameters_bhz_"//trim(file))
  ! write(10,nml=hkvars)
  ! close(10)
  close(50)




contains



  function hk_bhz(kx,ky) result(hk)
    real(8)                   :: kx,ky
    complex(8),dimension(4,4) :: hk
    Hk          = zero
    Hk(1:2,1:2) = hk_bhz2x2(kx,ky)
    Hk(3:4,3:4) = conjg(hk_bhz2x2(-kx,-ky))
    Hk(1,4) = -delta ; Hk(4,1)=-delta
    Hk(2,3) =  delta ; Hk(3,2)= delta
    Hk(1,3) = xi*rh*(sin(kx)-xi*sin(ky))
    Hk(3,1) =-xi*rh*(sin(kx)+xi*sin(ky))
  end function hk_bhz

  function hk_bhz2x2(kx,ky) result(hk)
    real(8)                   :: kx,ky
    complex(8),dimension(2,2) :: hk
    epsik   = cos(kx)+cos(ky)
    hk(1,1) = mh - epsik
    hk(2,2) =-hk(1,1)!-mh + epsik
    hk(1,2) = lambda*(sin(kx)-xi*sin(ky))
    hk(2,1) = lambda*(sin(kx)+xi*sin(ky))
  end function hk_bhz2x2


  function inverse_g0k(iw,hk) result(g0k)
    complex(8)                  :: iw
    complex(8),dimension(4,4)   :: hk
    complex(8),dimension(4,4)   :: g0k
    g0k=zero
    if(rh==0.AND.delta==0)then
       g0k(1:2,1:2) = inverse_g0k2x2(iw,hk(1:2,1:2))
       g0k(3:4,3:4) = inverse_g0k2x2(iw,hk(3:4,3:4))
    else
       g0k = -hk
       forall(i=1:4)g0k(i,i) = iw + xmu + g0k(i,i)
       call matrix_inverse(g0k)
    endif
  end function inverse_g0k
  !
  function inverse_g0k2x2(iw,hk) result(g0k)
    integer                     :: i
    complex(8),dimension(2,2)   :: hk
    complex(8)                  :: iw
    complex(8),dimension(2,2)   :: g0k
    complex(8)                  :: delta,ppi,vmix
    g0k=zero
    delta = iw - hk(1,1)
    ppi   = iw - hk(2,2)
    vmix  = -hk(1,2)
    g0k(1,1) = one/(delta - abs(vmix)**2/ppi)
    g0k(2,2) = one/(ppi - abs(vmix)**2/delta)
    g0k(1,2) = -vmix/(ppi*delta - abs(vmix)**2)
    g0k(2,1) = conjg(g0k(1,2))
  end function inverse_g0k2x2

  function Eigk(hk) result(eig)
    complex(8),dimension(4,4) :: hk
    real(8),dimension(4)      :: eig
    !call matrix_diagonalize(hk,eig)
    ! call matrix_diagonalize(hk(1:2,1:2),eig(1:2))
    ! call matrix_diagonalize(hk(3:4,3:4),eig(3:4))
    eig(1:2) = eigk2x2(hk(1:2,1:2))
    eig(3:4) = eigk2x2(hk(3:4,3:4))
  end function Eigk

  function Eigk2x2(hk) result(eig)
    complex(8),dimension(2,2) :: hk
    real(8),dimension(2)      :: eig
    complex(8)                :: delta0,vk2,sqrt0,u,v
    delta0 = (hk(1,1)-hk(2,2))
    vk2    = 4.d0*hk(1,2)*hk(2,1)
    sqrt0  = sqrt( abs(delta0)**2 + vk2 )
    eig(1)=hk(1,1)+hk(2,2) + sqrt0
    eig(2)=hk(1,1)+hk(2,2) - sqrt0
    eig = eig/2.d0
    if(sqrt0==0.d0)then
       u=sqrt(0.5d0)
       v=u
    else
       u = sqrt(0.5d0*(1.d0 + delta0/sqrt0))
       v = sqrt(0.5d0*(1.d0 - delta0/sqrt0))
    endif
    hk(1,1) = u
    hk(2,1) = -v
    hk(1,2) = v
    hk(2,2) = u
  end function Eigk2x2


  function matrix_uk(hk) result(uk)
    complex(8),dimension(2,2) :: hk
    complex(8),dimension(2,2) :: uk
    complex(8) :: delta0,ek,u,v
    delta0 = -(hk(1,1)-hk(2,2))
    ek     = sqrt( delta0**2 + 4.d0*hk(1,2)*hk(2,1) )
    u = sqrt(0.5d0*(1.d0 - delta0/ek))
    v = sqrt(0.5d0*(1.d0 + delta0/ek))
    uk(1,1) = u
    uk(1,2) = v
    uk(2,1) =-v
    uk(2,2) = u
  end function matrix_uk


end program tdd_pam


