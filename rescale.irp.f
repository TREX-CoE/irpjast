BEGIN_PROVIDER [ double precision, kappa ]
 implicit none
 BEGIN_DOC
 ! Constant in rescaling
 END_DOC
 kappa = 0.6d0
END_PROVIDER

BEGIN_PROVIDER [ double precision, kappa_inv ]
 implicit none
 BEGIN_DOC
 ! inverse of kappa
 END_DOC
 kappa_inv = 1.0d0 / kappa
END_PROVIDER

BEGIN_PROVIDER [ double precision, rescale_ee, (nelec, nelec) ]
 implicit none
 BEGIN_DOC
 ! R = (1 - exp(-kappa r))/kappa for electron-electron for $J_{ee}$
 END_DOC
 integer :: i, j

 do j = 1, nelec
   do i = 1, nelec
     rescale_ee(i, j) = (1.0d0 - dexp(-kappa * elec_dist(i, j))) * kappa_inv
   enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, rescale_en, (nelec, nnuc) ]
 implicit none
 BEGIN_DOC
 ! R = (1 - exp(-kappa r))/kappa for electron-nucleus for $J_{en}$
 END_DOC
 integer :: i, j

 do j = 1, nnuc
   do i = 1, nelec
     rescale_en(i, j) = (1.d0 - dexp(-kappa * elnuc_dist(i, j))) * kappa_inv
   enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_e, (nelec, nelec, 0:ncord)]
 implicit none
 BEGIN_DOC
 ! R = exp(-kappa r) for electron-electron for $J_{een}$
 END_DOC
 integer :: i, j, l
 double precision :: kappa_l

 do l=0,ncord
  kappa_l = -dble(l) * kappa
  do j = 1, nelec
    do i = 1, nelec
      rescale_een_e(i, j, l) = kappa_l * elec_dist(i, j)
    enddo
  enddo
 enddo
 rescale_een_e = dexp(rescale_een_e)
END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_n, (nelec, nnuc, 0:ncord)]
 implicit none
 BEGIN_DOC
 ! R = exp(-kappa r) for electron-electron for $J_{een}$
 END_DOC
 integer :: i, j, l
 double precision :: kappa_l

 do l=0,ncord
   kappa_l = - dble(l) * kappa
   do j = 1, nnuc
     do i = 1, nelec
       rescale_een_n(i, j, l) = kappa_l * elnuc_dist(i, j)
     enddo
   enddo
 enddo
 rescale_een_n = dexp(rescale_een_n)
END_PROVIDER
