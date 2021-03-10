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

BEGIN_PROVIDER [ double precision, rescale_ee, (nelec_8, nelec) ]
 implicit none
 BEGIN_DOC
 ! R = (1 - exp(-kappa r))/kappa for electron-electron for $J_{ee}$
 END_DOC
 integer :: i, j
 double precision :: x

 do j = 1, nelec
    do i = 1, j-1
       x = (1.0d0 - dexp(-kappa * elec_dist(i, j))) * kappa_inv
       rescale_ee(i, j) = x
       rescale_ee(j, i) = x
    enddo
    rescale_ee(j, j) = 0.d0
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, rescale_ee_deriv_e, (4, nelec, nelec) ]
 implicit none
 BEGIN_DOC
 ! R = (1 - exp(-kappa r))/kappa derived wrt x
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC
 integer :: i, j, ii
 double precision :: f

 do j = 1, nelec
    do i = 1, nelec
       f = 1.d0 - kappa*rescale_ee(i,j) ! == dexp(-kappa * elec_dist(i, j))
       do ii = 1, 4
          rescale_ee_deriv_e(ii, i, j) = elec_dist_deriv_e(ii, i, j)
       end do
       rescale_ee_deriv_e(4, i, j) = rescale_ee_deriv_e(4, i, j) + &
       (-kappa * rescale_ee_deriv_e(1, i, j) * rescale_ee_deriv_e(1, i, j)) + &
       (-kappa * rescale_ee_deriv_e(2, i, j) * rescale_ee_deriv_e(2, i, j)) + &
       (-kappa * rescale_ee_deriv_e(3, i, j) * rescale_ee_deriv_e(3, i, j))
       do ii = 1, 4
          rescale_ee_deriv_e(ii, i, j) = rescale_ee_deriv_e(ii, i, j) &
               * f
       enddo
    enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, rescale_en, (nelec_8, nnuc) ]
 implicit none
 BEGIN_DOC
 ! R = (1 - exp(-kappa r))/kappa for electron-nucleus for $J_{en}$
 END_DOC
 integer :: i, a

 do a = 1, nnuc
    do i = 1, nelec
       rescale_en(i, a) = (1.0d0 - dexp(-kappa * elnuc_dist(i, a))) * kappa_inv
    enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [ double precision, rescale_en_deriv_e, (4, nelec, nnuc) ]
 implicit none
 BEGIN_DOC
 ! R = (1 - exp(-kappa r))/kappa derived wrt x
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC
 integer :: i, ii, a
 double precision :: f

 do a = 1, nnuc
    do i = 1, nelec
       f = 1.d0 - kappa*rescale_en(i,a) ! == dexp(-kappa * elnuc_dist(i, a))
       do ii = 1, 4
          rescale_en_deriv_e(ii, i, a) = elnuc_dist_deriv_e(ii, i, a)
       end do
       rescale_en_deriv_e(4, i, a) = rescale_en_deriv_e(4, i, a) + &
       (-kappa * rescale_en_deriv_e(1, i, a) * rescale_en_deriv_e(1, i, a)) + &
       (-kappa * rescale_en_deriv_e(2, i, a) * rescale_en_deriv_e(2, i, a)) + &
       (-kappa * rescale_en_deriv_e(3, i, a) * rescale_en_deriv_e(3, i, a))
       do ii = 1, 4
          rescale_en_deriv_e(ii, i, a) = rescale_en_deriv_e(ii, i, a) &
               * f
       enddo
    enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_e, (nelec_8, nelec, 0:ncord)]
 implicit none
 BEGIN_DOC
 ! R = exp(-kappa r) for electron-electron for $J_{een}$
 END_DOC
 integer :: i, j, k, l
 double precision :: x
 double precision, parameter :: f = dexp(1.d0)

 rescale_een_e(:, :, 0) = 1.d0

 do l = 1, ncord
    k=0
    do j = 1, nelec
       do i = 1, j-1
          k = k+1
          x = rescale_een_e_ij(k,l)
          rescale_een_e(i, j, l) = x
          rescale_een_e(j, i, l) = x
       enddo
    enddo
 enddo

 do l = 0, ncord
    do j = 1, nelec
       rescale_een_e(j, j, l) = 0.d0
    enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_e_ij, (nelec*(nelec-1)/2, 0:ncord)]
 implicit none
 BEGIN_DOC
 ! R = exp(-kappa r) for electron-electron for $J_{een}$
 END_DOC
 integer :: i, j, l,k
 double precision :: x
 double precision, parameter :: f = dexp(1.d0)

 rescale_een_e_ij(:, 0) = 1.d0

 k=0
 do j = 1, nelec
    do i = 1, j-1
       k = k+1
       rescale_een_e_ij(k, 1) = dexp(-kappa * elec_dist(i, j))
    enddo
 enddo

 do l = 2, ncord
    do k=1,(nelec*nelec-nelec)/2
      rescale_een_e_ij(k, l) = rescale_een_e_ij(k, l-1) * rescale_een_e_ij(k, 1)
    enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_n, (nelec_8, nnuc, 0:ncord)]
 implicit none
 BEGIN_DOC
 ! R = exp(-kappa r) for electron-electron for $J_{een}$
 END_DOC
 integer :: i, a, l
 double precision :: x
 double precision, parameter :: f = dexp(1.d0)

 rescale_een_n(:,:,0) = 1.d0

 do a = 1, nnuc
    do i = 1, nelec
       rescale_een_n(i, a, 1) = dexp(-kappa * elnuc_dist(i, a))
    enddo
 enddo

 do l = 2, ncord
    do a = 1, nnuc
       do i = 1, nelec
          rescale_een_n(i, a, l) = rescale_een_n(i, a, l-1) * rescale_een_n(i, a, 1)
       enddo
    enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_n_deriv_e, (nelec_8, 4, nnuc, 0:ncord)]
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
          do ii = 1, 4
          ! r'(x) \lor r''(x)
             rescale_een_n_deriv_e(i, ii, a, l) = &
                  kappa_l * elnuc_dist_deriv_e(ii, i, a)
             !print *, "pp", ii, i, a, elnuc_dist_deriv_e(ii, i, a)
          enddo

          ! \left(r''(x)+r'(x)^2\right)
          rescale_een_n_deriv_e(i, 4, a, l) = rescale_een_n_deriv_e(i, 4, a, l) + &
          rescale_een_n_deriv_e(i, 1, a, l) * rescale_een_n_deriv_e(i, 1, a, l) + &
          rescale_een_n_deriv_e(i, 2, a, l) * rescale_een_n_deriv_e(i, 2, a, l) + &
          rescale_een_n_deriv_e(i, 3, a, l) * rescale_een_n_deriv_e(i, 3, a, l)

          ! \times e^{r(x)}
          do ii = 1, 4
             rescale_een_n_deriv_e(i, ii, a, l) = &
                 rescale_een_n_deriv_e(i, ii, a, l) * rescale_een_n(i, a, l)
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

BEGIN_PROVIDER [double precision, rescale_een_e_deriv_e, (nelec_8, 4, nelec, 0:ncord)]
 BEGIN_DOC
 ! Derivative of the scaled distance J_{een} wrt R_{ia}
 END_DOC
 implicit none

 integer :: i, ii, j, l
 double precision :: kappa_l

 !TODO: Check if rescale_een_e_deriv_e(:,:,0) = 0.d0
 do l = 0, ncord
    kappa_l = - dble(l) * kappa
    do j = 1, nelec
       do i = 1, nelec
          ! r'(x) \lor r''(x)
          do ii = 1, 4
             rescale_een_e_deriv_e(i, ii, j, l) = &
                 kappa_l * elec_dist_deriv_e(ii, i, j)
          enddo

          ! \left(r''(x)+r'(x)^2\right)
            rescale_een_e_deriv_e(i, 4, j, l) = rescale_een_e_deriv_e(i, 4, j, l) + &
            rescale_een_e_deriv_e(i, 1, j, l) * rescale_een_e_deriv_e(i, 1, j, l) + &
            rescale_een_e_deriv_e(i, 2, j, l) * rescale_een_e_deriv_e(i, 2, j, l) + &
            rescale_een_e_deriv_e(i, 3, j, l) * rescale_een_e_deriv_e(i, 3, j, l)

          ! \times e^{r(x)}
          do ii = 1, 4
             rescale_een_e_deriv_e(i, ii, j, l) = &
                 rescale_een_e_deriv_e(i, ii, j, l) * rescale_een_e(i, j, l)
          enddo
       enddo
    enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, rescale_een_e_deriv_e_t, (nelec_8, 4, nelec, 0:ncord)]
  implicit none
  BEGIN_DOC
! Transposed rescale_een_e_deriv_e
  END_DOC
  integer :: i,j,k,l
  do l=0,ncord
    do j=1,nelec
      do i=1,nelec
        rescale_een_e_deriv_e_t(j,1:4,i,l) = rescale_een_e_deriv_e(i,1:4,j,l)
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
