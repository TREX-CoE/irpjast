BEGIN_PROVIDER [double precision, factor_een]
 implicit none
 BEGIN_DOC
 ! Electron-electron nucleus contribution to Jastrow factor
 END_DOC
 integer :: i, j, alpha, p, k, l, lmax = 0
 factor_een = 0.0d0
 
 do alpha = 1, nnuc
    do j = 1, nelec
       do i = 1, nelec
        do p = 2, ncord
           do k = p - 1, 0
              if ( k == 0 ) then
                 lmax = p - k - 2
              else
                 lmax = p - k
              end if
              do l = lmax, 0
                 if ( mod(p - k - l, 2) == 0 ) then
                    factor_een = factor_een + cord_vect(p, k, l) * rescale_een_e(i, j) ** k &
                         * (rescale_een_n(i, alpha) ** l + rescale_een_n(j, alpha) ** l) *  &
                         (rescale_een_n(i, alpha) * rescale_een_n(j, alpha)) ** ((p - k - l) * 0.5d0)
                 end if
              end do
           end do
        end do
       end do
    end do
 end do

END_PROVIDER
