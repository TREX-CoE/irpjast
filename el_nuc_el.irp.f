BEGIN_PROVIDER [double precision, factor_een]
 implicit none
 BEGIN_DOC
 ! Electron-electron nucleus contribution to Jastrow factor
 END_DOC
 integer :: i, j, alpha, p, k, l, lmax, cindex
 double precision :: x, y, z, t, c_inv, u, a, b, a2, b2, c, t0

 PROVIDE cord_vect
 factor_een = 0.0d0

 do alpha = 1, nnuc
    do j = 1, nelec
       b = rescale_een_n(j, alpha)
       do i = 1, nelec
          u = rescale_een_e(i, j)
          a = rescale_een_n(i, alpha)
          a2 = a * a
          b2 = b * b
          c = rescale_een_n(i, alpha) * rescale_een_n(j, alpha)
          c_inv = 1.0d0 / c
          cindex = 0
          do p = 2, ncord
             x = 1.0d0
             do k = 0, p - 1
                if ( k /= 0 ) then
                   lmax = p - k
                else
                   lmax = p - k - 2
                end if
                t = x
                do l = 1, rshift(p - k, 1)
                  t = t * c
                end do
                ! We have suppressed this from the following loop:
                ! if ( iand(p - k - l, 1) == 0 ) then
                !
                ! Start from l=0 when p-k is even
                ! Start from l=1 when p-k is odd
                if (iand(p - k, 1) == 0) then
                  y = 1.0d0
                  z = 1.0d0
                else
                  y = a
                  z = b
                endif
                do l = iand(p - k, 1), lmax, 2
                   ! This can be used in case of a flatten cord_vect
                   ! cidx = 1 + l + (ncord + 1) * k + (ncord + 1) * (ncord + 1) * (p - 1) + &
                   !      (ncord + 1) * (ncord + 1) * ncord * (alpha - 1)
                   cindex = cindex + 1
                   factor_een = factor_een + cord_vect(cindex, typenuc_arr(alpha)) * (y + z) * t
                   t = t * c_inv
                   y = y * a2
                   z = z * b2
                end do
                x = x * u
             end do
          end do
       end do
    end do
 end do

 factor_een = 0.5d0 * factor_een

END_PROVIDER

BEGIN_PROVIDER [double precision, factor_een_naive]
 implicit none
 BEGIN_DOC
 ! Electron-electron nucleus contribution to Jastrow factor in a naive way
 END_DOC
 integer :: i, j, alpha, p, k, l, lmax, cindex
 double precision :: ria, rja, rij

 PROVIDE cord_vect
 factor_een_naive = 0.0d0
 
 do alpha = 1, nnuc
    do j = 2, nelec
       rja = rescale_een_n(j, alpha)
       do i = 1, j - 1
          ria = rescale_een_n(i, alpha)
          rij = rescale_een_e(i, j)
          cindex = 0
          do p = 2, ncord
             do k = p - 1, 0, -1
                if ( k /= 0 ) then
                   lmax = p - k
                else
                   lmax = p - k - 2
                end if
                do l = lmax, 0, -1
                   if ( iand(p - k - l, 1) == 0 ) then
                      cindex = cindex + 1
                      factor_een_naive = factor_een_naive + &
                           cord_vect(cindex, typenuc_arr(alpha)) * &
                           rij ** k * (ria ** l + rja ** l) * (ria * rja) ** rshift(p - k - l, 1)
                      !factor_een_naive = factor_een_naive + &
                      !     cord_vect(cindex, typenuc_arr(alpha)) * &
                      !     rij(i, j, k) * (ria(i, alpha, l) + rja(j, alpha, l)) &
                      !     * (ria(i, alpha, l) * rja(j, alpha, l)) ** rshift(p - k - l, 1)
                   end if
                end do
             end do
          end do
       end do
    end do
 end do

END_PROVIDER

!BEGIN_PROVIDER [double precision, factor_een_prog]
! implicit none
! BEGIN_DOC
! ! Electron-electron nucleus contribution to Jastrow factor in a naive way
! END_DOC
! integer :: alpha, i, j, p, k, l, lmax, m, cindex
! double precision :: ria, rja, rij, rij_inv
! double precision :: c, c_inv, t, x, y, z ! Placeholders for optimization
!                                                                                                 
! PROVIDE cord_vect
! factor_een_prog = 0.0d0
!                                                                                                 
! do alpha = 1, nnuc
!    do j = 2, nelec
!       rja = rescale_een_n(j, alpha)
!       do i = 1, j - 1
!          ria = rescale_een_n(i, alpha)
!          rij = rescale_een_e(i, j)
!          rij_inv =  1.0d0 / (rij * rij)
!          c = ria * rja
!          c_inv = 1.0d0 / c
!          cindex = 0
!          do p = 2, ncord
!
!             x = 1.0d0
!             do l = 1, p
!               x = x * rij
!             end do
!
!             do k = p - 1, 0, -1
!                if ( k /= 0 ) then
!                   lmax = p - k
!                else
!                   lmax = p - k - 2
!                end if
!
!                t = 1.0d0
!                do l = 1, rshift(p - k, 1)
!                  t = t * c
!                end do
!
!                do l = lmax, iand(p - k, 1), -2
!                   cindex = cindex + 1
!                   factor_een_prog = factor_een_prog + &
!                        cord_vect(cindex, typenuc_arr(alpha)) * &
!                        x * (ria ** l + rja ** l) * t
!                   t = t * c_inv
!                   x = x * rij_inv
!                end do
!
!             end do
!          end do
!       end do
!    end do
! end do
!
!END_PROVIDER

!BEGIN_PROVIDER [double precision, rij, (nelec, nelec, ncord)]
!&BEGIN_PROVIDER [double precision, ria, (nelec, nnuc, ncord)]
!&BEGIN_PROVIDER [double precision, rja, (nelec, nnuc, ncord)]
! BEGIN_DOC
! ! Tables with powers
! END_DOC
! integer :: i, j, k, alpha
! double precision :: x, y, z
!
! rij(:, :, :) = 0.0d0
! ria(:, :, :) = 0.0d0
! rja(:, :, :) = 0.0d0
! 
! implicit none
! do alpha = 1, nnuc
!    do j = 2, nelec
!       z = 1.0d0
!       do k = 1, ncord
!          rja(j, alpha, k) = z
!          z = z * rescale_een_n(j, alpha)
!       end do
!       do i = 1, j - 1
!          y = 1.0d0
!          do k = 1, ncord
!             ria(i, alpha, k) = y
!             y = y * rescale_een_n(i, alpha)
!          end do
!          x = 1.0d0
!          do k = 1, ncord
!             rij(i, j, k) = x
!             x = x * rescale_een_e(i, j)
!          end do
!       end do
!    end do
! end do
!
!END_PROVIDER
 
