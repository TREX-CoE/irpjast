BEGIN_PROVIDER [ double precision, jastrow_full ]
 implicit none
 BEGIN_DOC
 ! Complete jastrow factor 
 END_DOC

 print *, factor_ee
 print *, factor_en
 print *, factor_een

 jastrow_full = dexp(factor_ee + factor_en + factor_een)

END_PROVIDER

