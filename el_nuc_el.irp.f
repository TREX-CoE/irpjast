BEGIN_PROVIDER [double precision, factor_een]
 implicit none
 BEGIN_DOC
 ! Electron-electron nucleus contribution to Jastrow factor
 END_DOC
 integer :: i, j, alpha, p, k, l, lmax = 0
 double precision :: x, y, z, t, c_inv, u, a, b, a2, b2, c, t0
 factor_een = 0.0d0

 do alpha = 1, nnuc
    do j = 1, nelec
       b = rescale_een_n(j, alpha)
       do i = 1, nelec
        u = rescale_een_e(i,j)
        a = rescale_een_n(i, alpha)
        a2 = a*a
        b2 = b*b
        c = rescale_een_n(i, alpha) * rescale_een_n(j, alpha)
        c_inv = 1.d0/c
        do p = 2, ncord
           x = 1.d0
!           t0 = x* c** (rshift(p,1))
!           t = t0
           do k = 0, p - 1
              if ( k /= 0 ) then
                 lmax = p - k
              else
                 lmax = p - k - 2
              end if
              t = x* c** (rshift(p - k,1))
              ! We have suppressed this if from the following loop:
              ! if ( iand(p - k - l, 1) == 0 ) then
              !
              ! Start from l=0 when p-k is even
              ! Start from l=1 when p-k is odd
              if (iand(p-k,1) == 0) then
                y = 1.d0
                z = 1.d0
              else
                y = a
                z = b
              endif
!              print *, ''
              do l = iand(p-k,1), lmax, 2
!                 if (iand(p-k-l,1) == 0) then
                 factor_een = factor_een + cord_vect(l, k, p, alpha) * (y+z) * t
                 t = t * c_inv
                 y = y * a2
                 z = z * b2
              end do
              x = x * u
           end do
!           t0 = t0*c_inv
        end do
       end do
    end do
 end do

 factor_een = 0.5d0 * factor_een

END_PROVIDER
