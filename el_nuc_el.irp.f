BEGIN_PROVIDER [ double precision, factor_een ]
  implicit none
  BEGIN_DOC
  ! ElectronE-electron-nuclei contribution to Jastrow factor
  !
  ! 5436.20340250000
  END_DOC
  integer                        :: i, j, a, p, k, l, lmax, m, n
  double precision               :: cn, accu2, accu
  double precision               :: f(nnuc,0:ncord-2,0:ncord-2)
  double precision               :: tmp_c(nelec,nnuc,0:ncord,0:ncord-1)

!  factor_een = factor_een_blas
!  return

 factor_een = 0.0d0

 do p = 2, ncord
   do k = 0, p - 1
     if (k /= 0) then
       lmax = p - k
     else
       lmax = p - k - 2
     endif
     do l = 0, lmax
       if ( iand(p - k - l, 1) == 1) cycle

       m = (p - k - l) / 2

       do a = 1, nnuc
         accu2 = 0.d0
         cn = cord_vect_lkp(l, k, p, typenuc_arr(a))
         do j = 1, nelec
           accu = 0.d0
           do i = 1, nelec
             accu = accu +                                           &
                 rescale_een_e(i,j,k) *                              &
                 rescale_een_n(i,a,m) *                              &
                 rescale_een_n(j,a,m+l)
           enddo
           accu2 = accu2 + accu
         enddo
         factor_een = factor_een + accu2 * cn
       enddo

     enddo
   enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, factor_een_deriv_e, (4, nelec) ]
 implicit none
 BEGIN_DOC
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC
 integer :: i, ii, j, a, p, k, l, lmax, m
 double precision :: riam, rjam_cn, rial, rjal, rijk
 double precision, dimension(4) :: driam, drjam_cn, drial, drjal, drijk
 double precision :: cn, v1, v2, d1, d2, lap1, lap2

 factor_een_deriv_e = 0.0d0

 do p = 2, ncord
    do k = 0 , p - 1
       if (k /= 0) then
         lmax = p - k
       else
         lmax = p - k - 2
       endif

       do l = 0, lmax
          if ( iand(p - k - l, 1) == 1) cycle
          m = (p - k - l) / 2

          do a = 1, nnuc
             cn = cord_vect_lkp(l, k, p, typenuc_arr(a))

             do j = 1, nelec
                rjal = rescale_een_n(j, a, l)
                rjam_cn = rescale_een_n(j, a, m) * cn

                do ii = 1, 4
                   drjal(ii) = rescale_een_n_deriv_e(ii, j, a, l)
                   drjam_cn(ii) = rescale_een_n_deriv_e(ii, j, a, m) * cn
                enddo

                do i = 1, nelec
                   rial = rescale_een_n(i, a, l) + rjal
                   riam = rescale_een_n(i, a, m)
                   rijk = rescale_een_e(i, j, k)

                   do ii = 1, 4
                      drijk(ii) = rescale_een_e_deriv_e(ii, j, i, k)
                   enddo

                   v1 = rijk * rial    ! v(x)
                   v2 = rjam_cn * riam ! u(x)

                   lap1 = 0.0d0
                   lap2 = 0.0d0
                   do ii = 1, 3
                      d1 = drijk(ii) * rial + rijk * drjal(ii)
                      d2 = drjam_cn(ii) * riam
                      lap1 = lap1 + d1 * d2
                      lap2 = lap2 + drijk(ii) * drjal(ii)
                      factor_een_deriv_e(ii, j) = factor_een_deriv_e(ii, j) + v1 * d2 + d1 * v2
                   enddo

                   ! v(x) u''(x) + 2 * u'(x) v'(x) + u(x) v''(x)
                   ii = 4
                   d1 = drijk(ii) * rial + rijk * drjal(ii) + lap2 + lap2
                   d2 = drjam_cn(ii) * riam
                   factor_een_deriv_e(ii, j) = factor_een_deriv_e(ii, j) + v1 * d2 + d1 * v2 + lap1 + lap1

                enddo
             enddo
          enddo
       enddo
    enddo
 enddo

END_PROVIDER
