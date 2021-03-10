integer function size_8(n) 
 implicit none
 integer, intent(in) :: n
 integer :: n8

 n8 = ((n-1)/8+1) * 8
 if (popcnt(n8) == 1) then 
    ! Power of two, shift by 8
    n8 = n8 + 8
 endif
 size_8 = n8
end

BEGIN_PROVIDER [ integer, nelec_8 ]
 implicit none
 integer, external :: size_8
 nelec_8 = size_8(nelec)
END_PROVIDER




BEGIN_PROVIDER [ integer, nelec ]
 implicit none
 BEGIN_DOC
 ! Number of electrons
 END_DOC

 character*(32) :: buffer
 integer, external :: iargc
 if (iargc() == 0) then
   nelec = 10
 else
  call getarg(1,buffer)
  read(buffer,*)nelec
 endif

END_PROVIDER

BEGIN_PROVIDER [ integer, nelec_up ]
 implicit none
 BEGIN_DOC
 ! Number of alpha and beta electrons
 END_DOC
 nelec_up = nelec/2
END_PROVIDER


BEGIN_PROVIDER [ double precision, elec_coord, (nelec, 3) ]
 implicit none
 BEGIN_DOC
 ! Electron coordinates
 END_DOC
 character(len=*), parameter :: FILE_NAME = "elec_coord.txt"
 integer :: fu, rc, i, j

 open(action='read', file=FILE_NAME, iostat=rc, newunit=fu)

 do i = 1, nelec
    read(fu, *) elec_coord(i, :)
 end do

 close(fu)

END_PROVIDER

BEGIN_PROVIDER [ double precision, elec_dist, (nelec, nelec) ]
 implicit none
 BEGIN_DOC
 ! e-e distance
 END_DOC
 integer :: i, j
 double precision :: x, y, z

 do j = 1, nelec
    do i = 1, nelec
       x = elec_coord(i, 1) - elec_coord(j, 1)
       y = elec_coord(i, 2) - elec_coord(j, 2)
       z = elec_coord(i, 3) - elec_coord(j, 3)
       elec_dist(i, j) = dsqrt( x*x + y*y + z*z )
    enddo
!    elec_dist(j, j) = 1.d-10
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, asymp_jasb, (2)]
 BEGIN_DOC
 ! Asymptotic component subtracted from J_ee
 END_DOC
 implicit none
 integer :: i, p
 double precision :: asym_one, x

 asym_one = bord_vect(1) * kappa_inv / (1.0d0 + bord_vect(2) * kappa_inv)
 asymp_jasb(:) = (/asym_one, 0.5d0 * asym_one/)
 
 do i = 1, 2
    x = kappa_inv
    do p = 2, nbord
       x = x * kappa_inv
       asymp_jasb(i) = asymp_jasb(i) + bord_vect(p + 1) * x
    end do
 end do

END_PROVIDER

BEGIN_PROVIDER [double precision, factor_ee]
 implicit none
 BEGIN_DOC
 ! Electron-electron contribution to Jastrow factor
 END_DOC
 integer :: i, j, p, ipar
 double precision :: pow_ser, x, spin_fact

 factor_ee = 0.0d0

 do j = 1, nelec
    do i = 1, j - 1
       x = rescale_ee(i, j) 
       pow_ser = 0.0d0
       spin_fact = 1.0d0
       ipar = 1

       do p = 2, nbord
          x = x * rescale_ee(i, j)
          pow_ser = pow_ser + bord_vect(p + 1) * x
       end do

       if (j.le.nelec_up .or. i.gt.nelec_up) then
           spin_fact = 0.5d0
	   ipar = 2
       end if

       factor_ee = factor_ee + spin_fact * bord_vect(1) * rescale_ee(i, j) &
            / (1.0d0 + bord_vect(2) * rescale_ee(i, j)) - asymp_jasb(ipar) + pow_ser

    end do
 end do

END_PROVIDER

BEGIN_PROVIDER [double precision, factor_ee_deriv_e, (4, nelec) ]
 implicit none
 BEGIN_DOC
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC
 integer :: i, ii, j, p
 double precision :: x, x_inv, y, den, invden, lap1, lap2, lap3, third, spin_fact
 double precision, dimension(3) :: pow_ser_g
 double precision, dimension(4) :: dx

 factor_ee_deriv_e = 0.0d0
 third = 1.0d0 / 3.0d0

 do j = 1 , nelec
    do i = 1, nelec
       pow_ser_g = 0.0d0
       spin_fact = 1.0d0
       den = 1.0d0 + bord_vect(2) * rescale_ee(i, j)
       invden = 1.0d0 / den
       x_inv = 1.0d0 / (rescale_ee(i, j) + 1.0d-18)

       do ii = 1, 4
          dx(ii) = rescale_ee_deriv_e(ii, j, i)
       enddo

       if ((i.le.nelec_up .and. j.le.nelec_up) .or. &
            (i.gt.nelec_up .and. j.gt.nelec_up)) then
           spin_fact = 0.5d0
       end if

       lap1 = 0.0d0
       lap2 = 0.0d0
       lap3 = 0.0d0
       do ii = 1, 3
          x = rescale_ee(i, j)
          do p = 2, nbord
             ! p a_{p+1} r[i,j]^(p-1)
             y = p * bord_vect(p + 1) * x
             pow_ser_g(ii) += y * dx(ii)
             ! (p-1) p a_{p+1} r[i,j]^(p-2) r'[i,j]^2
             lap1 += (p - 1) * y * x_inv * dx(ii) * dx(ii)
             ! p a_{p+1} r[i,j]^(p-1) r''[i,j]
             lap2 += y
             x = x * rescale_ee(i, j)
          end do

          ! (a1 (-2 a2 r'[i,j]^2+(1+a2 r[i,j]) r''[i,j]))/(1+a2 r[i,j])^3
          lap3 += -2.0d0 * bord_vect(2) * dx(ii) * dx(ii)

          ! \frac{a1 * r'(i,j)}{(a2 * r(i,j)+1)^2}
          factor_ee_deriv_e(ii, j) += spin_fact * bord_vect(1) &
               * dx(ii) * invden * invden + pow_ser_g(ii)
       enddo

       ii = 4
       lap2 *= dx(ii) * third
       lap3 += den * dx(ii)
       lap3 *= spin_fact * bord_vect(1) * invden * invden * invden
       factor_ee_deriv_e(ii, j) += lap1 + lap2 + lap3

    end do
 end do

END_PROVIDER
