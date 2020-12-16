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

 do l = 0, ncord
    kappa_l = -dble(l) * kappa
    do j = 1, nelec
       do i = 1, nelec
          rescale_een_e(i, j, l) = kappa_l * elec_dist(i, j)
       enddo
    enddo
 enddo

 rescale_een_e = dexp(rescale_een_e)

 do l = 0, ncord
    do j = 1, nelec
       rescale_een_e(j, j, l) = 0.d0
    enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_n, (nelec, nnuc, 0:ncord)]
 implicit none
 BEGIN_DOC
 ! R = exp(-kappa r) for electron-electron for $J_{een}$
 END_DOC
 integer :: i, a, l
 double precision :: kappa_l

 do l = 0, ncord
    kappa_l = - dble(l) * kappa
    do a = 1, nnuc
       do i = 1, nelec
          rescale_een_n(i, a, l) = kappa_l * elnuc_dist(i, a)
       enddo
    enddo
 enddo

 rescale_een_n = dexp(rescale_een_n)

END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_n_deriv_e, (4, nelec, nnuc, 0:ncord)]
 implicit none
 BEGIN_DOC
 ! Derivative of the scaled distance J_{een} wrt R_{ia}
 END_DOC
 integer :: i, ii, j, l, a
 double precision :: kappa_l

 do l = 0, ncord
    kappa_l = - dble(l) * kappa
    do a = 1, nnuc
       do i = 1, nelec
          ! r'(x) \lor r''(x)
          do ii = 1, 4
             rescale_een_n_deriv_e(ii, i, a, l) = &
                  kappa_l * elnuc_dist_deriv_e(ii, i, a)
             !print *, "pp", ii, i, a, elnuc_dist_deriv_e(ii, i, a)
          enddo

          ! \left(r''(x)+r'(x)^2\right)
          rescale_een_n_deriv_e(4, i, a, l) = rescale_een_n_deriv_e(4, i, a, l) + &
          rescale_een_n_deriv_e(1, i, a, l) * rescale_een_n_deriv_e(1, i, a, l) + &
          rescale_een_n_deriv_e(2, i, a, l) * rescale_een_n_deriv_e(2, i, a, l) + &
          rescale_een_n_deriv_e(3, i, a, l) * rescale_een_n_deriv_e(3, i, a, l)

          ! \times e^{r(x)}
          do ii = 1, 4
             rescale_een_n_deriv_e(ii, i, a, l) = &
                 rescale_een_n_deriv_e(ii, i, a, l) * rescale_een_n(i, a, l)
          enddo
       enddo
    enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, elnuc_dist_deriv_e, (4, nelec, nnuc)]
 BEGIN_DOC
 ! Derivative of R_{ia} wrt x
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC
 implicit none
 integer :: i, ii, a
 double precision :: ria_inv

 do a = 1, nnuc 
    do i = 1, nelec
       ria_inv = 1.0d0 / elnuc_dist(i, a)
       do ii = 1, 3
          ! \frac{x-x0}{\sqrt{c+(x-x0)^2}}
          elnuc_dist_deriv_e(ii, i, a) = (elec_coord(i, ii) - nuc_coord(a, ii)) * ria_inv
       end do
       elnuc_dist_deriv_e(4, i, a) = 2.0d0 * ria_inv
    end do
 end do
END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_e_deriv_e, (4, nelec, nelec, 0:ncord)]
 BEGIN_DOC
 ! Derivative of the scaled distance J_{een} wrt R_{ia}
 END_DOC
 implicit none

 integer :: i, ii, j, l
 double precision :: kappa_l

 do l = 0, ncord
    kappa_l = - dble(l) * kappa
    do j = 1, nelec
       do i = 1, nelec
          ! r'(x) \lor r''(x)
          do ii = 1, 4
             rescale_een_e_deriv_e(ii, i, j, l) = &
                 kappa_l * elec_dist_deriv_e(ii, i, j)
          enddo

          ! \left(r''(x)+r'(x)^2\right)
          rescale_een_e_deriv_e(4, i, j, l) = rescale_een_e_deriv_e(4, i, j, l) + &
            rescale_een_e_deriv_e(1, i, j, l) * rescale_een_e_deriv_e(1, i, j, l) + &
            rescale_een_e_deriv_e(2, i, j, l) * rescale_een_e_deriv_e(2, i, j, l) + &
            rescale_een_e_deriv_e(3, i, j, l) * rescale_een_e_deriv_e(3, i, j, l)

          ! \times e^{r(x)}
          do ii = 1, 4
             rescale_een_e_deriv_e(ii, i, j, l) = &
                 rescale_een_e_deriv_e(ii, i, j, l) * rescale_een_e(i, j, l)
          enddo
       enddo
    enddo
 enddo
END_PROVIDER
 
BEGIN_PROVIDER [double precision, elec_dist_deriv_e, (4, nelec, nelec)]
 BEGIN_DOC
 ! Derivative of R_{ij} wrt x
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z = 2/rij
 END_DOC
 implicit none
 integer :: i, ii, j
 double precision :: rij_inv

 do j = 1, nelec
    do i = 1, nelec
       rij_inv = 1.0d0 / elec_dist(i, j)
       do ii = 1, 3
          ! \frac{x-x0}{\sqrt{c+(x-x0)^2}}
          elec_dist_deriv_e(ii, i, j) = (elec_coord(i, ii) - elec_coord(j, ii)) * rij_inv
       end do
       elec_dist_deriv_e(4, i, j) = 2.0d0 * rij_inv
    end do
    elec_dist_deriv_e(:, j, j) = 0.0d0
 end do

END_PROVIDER
