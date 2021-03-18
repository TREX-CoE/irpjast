BEGIN_PROVIDER [ double precision, factor_een ]
   implicit none
   BEGIN_DOC
   ! Electron -electron-nuclei contribution to Jastrow factor
   !
   ! 5436.20340250000
   END_DOC
   integer                        :: i, j, a, p, k, l, lmax, m, n
   double precision               :: cn, accu2, accu

   factor_een = 0.0d0

   do n = 1, dim_cord_vect

     l = lkpm_of_cindex(1,n)
     k = lkpm_of_cindex(2,n)
     p = lkpm_of_cindex(3,n)
     m = lkpm_of_cindex(4,n)

     do a = 1, nnuc
       accu2 = 0.d0
       cn = cord_vect_full(n, a)
       do j = 1, nelec
         accu = 0.d0
         do i = 1, nelec
           accu = accu +                                             &
               rescale_een_e(i,j,k) *                                &
               rescale_een_n(i,a,m)
         enddo
         accu2 = accu2 + accu*rescale_een_n(j,a,m+l)
       enddo
       factor_een = factor_een + accu2 * cn
     enddo

   enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, factor_een_deriv_e, (nelec,4) ]
  implicit none
  BEGIN_DOC
! Derivative of the Jeen
! 35533.115255
  END_DOC
  integer                        :: i, j, a, p, k, l, lmax, m, n
  double precision               :: cn, accu, accu2, daccu(1:4), daccu2(1:4)
  
!  factor_een_deriv_e(1:nelec,1:4) = factor_een_deriv_e_blas(1:4,1:nelec)
!  return
  
  factor_een_deriv_e(1:nelec,1:4) = 0.0d0
  
  do n = 1, dim_cord_vect
    
    l = lkpm_of_cindex(1,n)
    k = lkpm_of_cindex(2,n)
    p = lkpm_of_cindex(3,n)
    m = lkpm_of_cindex(4,n)
    
    do a = 1, nnuc
      cn = cord_vect_full(n, a)
      do j = 1, nelec
        accu=0.d0
        accu2 = 0.d0
        daccu (1:4) = 0.d0
        daccu2(1:4) = 0.d0
        do i = 1, nelec
          accu = accu +                                              &
              rescale_een_e(i,j,k) *                                 &
              rescale_een_n(i,a,m)
          accu2 = accu2 +                                            &
              rescale_een_e(i,j,k) *                                 &
              rescale_een_n(i,a,m+l)
          daccu(1:4) = daccu(1:4) +                                  &
              rescale_een_e_deriv_e_t(i,1:4,j,k) *                   &
              rescale_een_n(i,a,m)
          daccu2(1:4) = daccu2(1:4) +                                &
              rescale_een_e_deriv_e_t(i,1:4,j,k) *                   &
              rescale_een_n(i,a,m+l)
          
        enddo
        factor_een_deriv_e(j,1:4) = factor_een_deriv_e(j,1:4) +      &
            (accu * rescale_een_n_deriv_e(j,1:4,a,m+l) + daccu(1:4) * rescale_een_n(j,a,m+l) +&
            daccu2(1:4)* rescale_een_n(j,a,m) + accu2*rescale_een_n_deriv_e(j,1:4,a,m)) * cn
        
        factor_een_deriv_e(j,4) = factor_een_deriv_e(j,4) + 2.d0*(   &
            daccu (1) * rescale_een_n_deriv_e(j,1,a,m+l) +           &
            daccu (2) * rescale_een_n_deriv_e(j,2,a,m+l) +           &
            daccu (3) * rescale_een_n_deriv_e(j,3,a,m+l) +           &
            daccu2(1) * rescale_een_n_deriv_e(j,1,a,m  ) +           &
            daccu2(2) * rescale_een_n_deriv_e(j,2,a,m  ) +           &
            daccu2(3) * rescale_een_n_deriv_e(j,3,a,m  ) )*cn
      enddo
    enddo
  enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, factor_een_deriv_e_ref, (nelec,4) ]
 implicit none
 BEGIN_DOC
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC
 integer :: i, ii, j, a, p, k, l, lmax, m
 double precision :: riam, rjam_cn, rial, rjal, rijk
 double precision, dimension(4) :: driam, drjam_cn, drial, drjal, drijk
 double precision :: cn, v1, v2, d1, d2, lap1, lap2

 factor_een_deriv_e_ref = 0.0d0

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
                   drjal(ii) = rescale_een_n_deriv_e(j, ii, a, l)
                   drjam_cn(ii) = rescale_een_n_deriv_e(j, ii, a, m) * cn
                enddo

                do i = 1, nelec
                   rial = rescale_een_n(i, a, l) + rjal
                   riam = rescale_een_n(i, a, m)
                   rijk = rescale_een_e(i, j, k)

                   do ii = 1, 4
                      drijk(ii) = rescale_een_e_deriv_e(j, ii, i, k)
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
                      factor_een_deriv_e_ref(j,ii) = factor_een_deriv_e_ref(j,ii) + v1 * d2 + d1 * v2
                   enddo

                   ! v(x) u''(x) + 2 * u'(x) v'(x) + u(x) v''(x)
                   ii = 4
                   d1 = drijk(ii) * rial + rijk * drjal(ii) + lap2 + lap2
                   d2 = drjam_cn(ii) * riam
                   factor_een_deriv_e_ref(j,ii) = factor_een_deriv_e_ref(j,ii) + v1 * d2 + d1 * v2 + lap1 + lap1

                enddo
             enddo
          enddo
       enddo
    enddo
 enddo

END_PROVIDER
