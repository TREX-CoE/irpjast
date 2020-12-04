BEGIN_PROVIDER [ double precision, jastrow_full ]
 implicit none
 BEGIN_DOC
 ! Complete jastrow factor 
 END_DOC
 ! integer :: i, j, k, l

 ! do l = 1, nnuc
 !    do k = 1, ncord
 !       do j = 0, ncord
 !          do i = 0, ncord
 !             write(*, *) cord_vect_0(i, j, k, l)
 !          end do
 !       end do
 !    end do
 ! end do


 print *, factor_ee
 print *, factor_en
 print *, factor_een

 jastrow_full = dexp(factor_ee + factor_en + factor_een)

END_PROVIDER

