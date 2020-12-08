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
