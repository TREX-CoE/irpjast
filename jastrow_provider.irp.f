BEGIN_PROVIDER [ double precision, jastrow_full ]
 implicit none
 BEGIN_DOC
 ! Complete jastrow factor 
 END_DOC
 integer :: i, j

 print *, "J_ee = ", factor_ee
 print *, "J_en = ", factor_en
 print *, "J_een = ", factor_een
 print *, "J_enn_naive = ", factor_een_naive
 print *, "J_enn_prog = ", factor_een_prog
 print *, "J = J_ee + J_en + J_een = ", factor_ee + factor_en + factor_een

 jastrow_full = dexp(factor_ee + factor_en + factor_een)

END_PROVIDER
