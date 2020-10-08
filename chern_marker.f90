!------------------------------------------------------------------
! builds full_p_projector and full_q_projector
! we pass to real space representation and calculate the sum over ky
! of p_ky and q_ky, which are nx*nky x nx*nky matrices
! uses p_projectors and q_projectors at all k points (from the module)
!------------------------------------------------------------------

subroutine buildfullprojectors (p_projector,q_projector, nky,nx)
  implicit none
  integer, intent(in)           :: nky,nx
  complex(kind=idp), intent(in) :: p_projector(nky,nx,nx), q_projector(nky,nx,nx)
  integer                       :: ik, jy1, jy2, jx1, jx2, jr1, jr2
  complex(kind=idp)             :: bloch_phase
  full_p_projector = zzero
  full_q_projector = zzero
  loop_ky: do ik = 1,nky
     ky = (ik-1)*delta_k
     do jy1 = 1,nky
        do jy2 = 1,nky
           do jx1 = 1,nx
              do jx2 = 1,nx
                 jr1 = nx*(jy1-1) + jx1
                 jr2 = nx*(jy2-1) + jx2
                 bloch_phase = zexp( im * ky * (jy1 + coordy(jx1) - jy2 - coordy(jx2)) )
                 full_p_projector(jr1,jr2) = full_p_projector(jr1,jr2) +
                 bloch_phase * p_projector(ik,jx1,jx2)
                 full_q_projector(jr1,jr2) = full_q_projector(jr1,jr2) + bloch_phase*q_projector(ik,jx1,jx2)
              enddo
           enddo
        enddo
     enddo
  enddo loop_ky
  full_p_projector = full_p_projector * one_over_nky
  full_q_projector = full_q_projector * one_over_nky
  return
end subroutine buildfullprojectors

!------------------------------------------------------------------
! calcola il numero di chern localmente nello spazio reale,
! alla bianco & resta, per strip obc su x e pbc su y.
! prende in input i proiettori p e q sulle bande occupate complessive,
! che sono matrici nx*nky x nx*nky
! there are 5 versions (experimentally, they produce the same numbers!):
! get_local_chern (like get_local_chern_xy but less memory-consuming)
! get_local_chern_nomatmul (does not use matmul: 1.5 times more time-consuming)
! get_local_chern_xy (the unsymmetric one)
! get_local_chern_symm
! get_local_chern_comm
!------------------------------------------------------------------

!------------------------------------------------------------------
! get_local_chern is the unsymmetric version requiring less memory
! c(r) = -4 pi (2/v_uc) im tr (x_p y_q)
!------------------------------------------------------------------

subroutine get_local_chern (nx, nky, full_p_projector,
  full_q_projector, chern_marker, totalchern)
  implicit none
  integer,           intent(in)    :: nx, nky
  real(kind=idp),    intent(out)   :: chern_marker(nx,nky), totalchern
  complex(kind=idp), intent(inout) :: full_p_projector(nx*nky,nx*nky),full_q_projector(nx*nky,nx*nky)
  integer                          :: jx1, jx2, jx, jy, jy1, jy2, jr1, jr2, jr, js

  !...> now we construct the matrix t (puts in --> full_q_projector)
  do jy1 = 1,nky
     do jy2 = 1,nky
        do jx1 = 1,nx
           do jx2 = 1,nx
              jr1 = nx*(jy1-1)+jx1
              jr2 = nx*(jy2-1)+jx2
              full_q_projector(jr1,jr2) = coordx(jx1) * (coordy(jx2) +
              jy2-1) * full_q_projector(jr1,jr2)
           enddo
        enddo
     enddo
  enddo

  !....>  only diagonal elements are needed, but still matmul is much
  !more efficient!
  full_q_projector =  matmul (full_p_projector, matmul(full_q_projector, full_p_projector) )
  do jy = 1,nky
     do jx = 1,nx
        chern_marker(jx,jy) = -dimag(full_q_projector(nx*(jy-1)+jx,
        nx*(jy-1)+jx)) * four * twopi / unit_cell_area
     enddo
  enddo

  totalchern = zero
  do jr =1,nx*nky
     totalchern = totalchern + dimag(full_q_projector(jr,jr))
  enddo
  ! there is an extra 1/n_sites=1/(nx*nky) in the totalchern
  totalchern = -totalchern * four * twopi * one_over_nky / (nx*unit_cell_area)

  return
end subroutine get_local_chern
