BEGIN_PROVIDER [double precision, factor_een]
 implicit none
 BEGIN_DOC
 ! Electron-electron nucleus contribution to Jastrow factor
 END_DOC
 integer :: i, j, alpha, p, k, l, lmax = 0
 double precision :: x, y, z, t, t_inv
 factor_een = 0.0d0
 
 do alpha = 1, nnuc
    do j = 1, nelec
       do i = 1, nelec
        t_inv = 1.d0/(rescale_een_n(i, alpha) * rescale_een_n(j, alpha))
        do p = 2, ncord
           do k = p - 1, 0, -1
              if ( k == 0 ) then
                 lmax = p - k - 2
              else
                 lmax = p - k
              end if
              x = 1.d0
              y = 1.d0
              z = 1.d0
              t = (rescale_een_n(i, alpha) * rescale_een_n(j, alpha)) ** (rshift(p - k,1))
              do l = 0, lmax
                 if ( iand(p - k - l, 1) == 0 ) then
                    factor_een = factor_een + cord_vect(l, k, p, alpha) * x &
                         * (y + z) * t
                    t = t * t_inv
                 end if
                 x = x * rescale_een_e(i, j)
                 y = y * rescale_een_n(i, alpha) 
                 z = z * rescale_een_n(j, alpha) 
              end do
           end do
        end do
       end do
    end do
 end do

 factor_een = 0.5d0 * factor_een

END_PROVIDER
