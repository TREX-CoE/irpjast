BEGIN_PROVIDER [ double precision, jastrow_full ]
 implicit none
 BEGIN_DOC
 ! Complete jastrow factor 
 END_DOC

 if (ncord == 0) then
    jastrow_full = factor_ee + factor_en
 else
    jastrow_full = factor_ee + factor_en + factor_een
 endif

 !print *, "J_ee = ", factor_ee
 !print *, "J_en = ", factor_en
 !print *, "J_een = ", factor_een
 !print *, "J = J_ee + J_en + J_een = ", factor_ee + factor_en + factor_een

END_PROVIDER

BEGIN_PROVIDER [ double precision, jastrow_derivs, (4, nelec) ]
 implicit none
 BEGIN_DOC
 ! Gradient and Laplacian
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC

 if (ncord == 0) then
    jastrow_derivs = factor_ee_deriv_e + factor_en_deriv_e
 else
    jastrow_derivs = factor_ee_deriv_e + factor_en_deriv_e + factor_een_deriv_e
 endif

 !print *, "\nabla J", jastrow_derivs(1:3, :)
 !print *, "\nabla^2 J = ", jastrow_derivs(4, :)

END_PROVIDER
