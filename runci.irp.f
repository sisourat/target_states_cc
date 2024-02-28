!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gen_csfs(nmo_occ,nmo_act,nmo_virt)
 use mycsf_lib
 implicit none

 integer, intent(in) :: nmo_occ,nmo_act,nmo_virt
 integer :: nalpha, nbeta
 integer :: i, j, k

 integer, allocatable :: occ_a(:)
 integer, allocatable :: occ_b(:)

! only true for CIS of closed-shell systems
 myncsfmax = nmo_act*nmo_virt+1
 myndetmax = 2*(ncsfmax-1)+1 
 nalpha = nmo_occ
 nbeta = nmo_occ

 allocate(mycsf_basis(N_int,2, myndetmax, myncsfmax))
 allocate(mycoef_det(2,myncsfmax))

 allocate(occ_a(nmo_occ),occ_b(nmo_occ))

 mycsf_basis(:,:,:,:) = 0_bit_kind
 
! HF det
 mycoef_det(1,1) = 1.0d0
 do i = 1, nmo_occ
   occ_a(i) = i
   occ_b(i) = i
   call create_det(nalpha,nbeta,occ_a(1),occ_b(1),mycsf_basis(1,1,j,i))
 enddo
 
! do i = 1, ncsfmax ! first loop on the csf of the space ispace 
!   coef_det(1,i) =  0.70710678118654746
!   coef_det(2,i) =  0.70710678118654746
!   occ_a(
! enddo

 deallocate(occ_a,occ_b)

end subroutine gen_csfs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine free_csfs
 use mycsf_lib
 deallocate(mycsf_basis,mycoef_det)
end subroutine free_csfs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

