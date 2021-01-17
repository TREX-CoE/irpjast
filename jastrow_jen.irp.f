BEGIN_PROVIDER [ double precision, elnuc_dist, (nelec, nnuc) ]
 implicit none
 BEGIN_DOC
 ! e-n distance
 END_DOC
 integer :: i, j
 double precision :: x, y, z
 do j = 1, nnuc
    do i = 1, nelec
       x = elec_coord(i, 1) - nuc_coord(j, 1)
       y = elec_coord(i, 2) - nuc_coord(j, 2)
       z = elec_coord(i, 3) - nuc_coord(j, 3)
       elnuc_dist(i, j) = dsqrt( x*x + y*y + z*z )
    enddo
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, factor_en]
 implicit none
 BEGIN_DOC
 ! Electron-nuclei contribution to Jastrow factor
 END_DOC
 integer :: i, a, p
 double precision :: pow_ser

 factor_en = 0.0d0

 do a = 1 , nnuc
    do i = 1, nelec
       pow_ser = 0.0d0

       do p = 2, naord
          pow_ser = pow_ser + aord_vect(p + 1, typenuc_arr(a)) * rescale_en_stored(p, i, a)
       end do

       factor_en = factor_en + aord_vect(1, typenuc_arr(a)) * rescale_en(i, a) &
            / (1.0d0 + aord_vect(2, typenuc_arr(a)) * rescale_en(i, a)) + pow_ser

    end do
 end do

END_PROVIDER

BEGIN_PROVIDER [double precision, factor_en_deriv_e, (4, nelec) ]
 implicit none
 BEGIN_DOC
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC
 integer :: i, ii, a, p
 double precision :: y, den, invden, lap1, lap2, lap3, third
 double precision, dimension(3) :: pow_ser_g
 double precision, dimension(4) :: dx

 factor_en_deriv_e = 0.0d0
 third = 1.0d0 / 3.0d0

 do a = 1 , nnuc
    do i = 1, nelec
       pow_ser_g = 0.0d0
       den = 1.0d0 + aord_vect(2, typenuc_arr(a)) * rescale_en(i, a)
       invden = 1.0d0 / den

       do ii = 1, 4
          dx(ii) = rescale_en_deriv_e(ii, i, a)
       enddo

       lap1 = 0.0d0
       lap2 = 0.0d0
       lap3 = 0.0d0
       do ii = 1, 3
          do p = 2, naord
             ! p a_{p+1} r[i,a]^(p-1)
             y = p * aord_vect(p + 1, typenuc_arr(a)) * rescale_en_stored(p - 1, i, a)
             pow_ser_g(ii) += y * dx(ii)
             ! (p-1) p a_{p+1} r[i,a]^(p-2) r'[i,a]^2
             lap1 += (p - 1) * p * aord_vect(p + 1, typenuc_arr(a)) * &
		     rescale_en_stored(p - 2, i, a) * dx(ii) * dx(ii)
             ! p a_{p+1} r[i,a]^(p-1) r''[i,a]
             lap2 += y
          end do

          ! (a1 (-2 a2 r'[i,a]^2+(1+a2 r[i,a]) r''[i,a]))/(1+a2 r[i,a])^3
          lap3 += -2.0d0 * aord_vect(2, typenuc_arr(a)) * dx(ii) * dx(ii)

          ! \frac{a1 * r'(i,a)}{(a2 * r(i,a)+1)^2}
          factor_en_deriv_e(ii, i) += aord_vect(1, typenuc_arr(a)) &
               * dx(ii) * invden * invden + pow_ser_g(ii)
       enddo

       ii = 4
       lap2 *= dx(ii) * third
       lap3 += den * dx(ii)
       lap3 *= aord_vect(1, typenuc_arr(a)) * invden * invden * invden
       factor_en_deriv_e(ii, i) += lap1 + lap2 + lap3

    end do
 end do

END_PROVIDER
