BEGIN_PROVIDER [ double precision, jastrow_full ]
 implicit none
 BEGIN_DOC
 ! Complete jastrow factor 
 END_DOC
 integer :: i, j

 print *, "J_ee = ", factor_ee
 print *, "J_en = ", factor_en
 print *, "J_een = ", factor_een
 print *, "J = J_ee + J_en + J_een = ", factor_ee + factor_en + factor_een
 print *, "\nabla_ix J_een", factor_een_deriv_e(1, :)
 print *, "\nabla_iy J_een", factor_een_deriv_e(2, :)
 print *, "\nabla_iz J_een", factor_een_deriv_e(3, :)
 print *, "\nabla_i^2 J_een", factor_een_deriv_e(4, :)

 jastrow_full = dexp(factor_ee + factor_en + factor_een)

END_PROVIDER
