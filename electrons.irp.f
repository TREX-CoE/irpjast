BEGIN_PROVIDER [ integer, nelec ]
 implicit none
 BEGIN_DOC
 ! Number of electrons
 END_DOC
 nelec = 10
END_PROVIDER

BEGIN_PROVIDER [ integer, nelec_up ]
 implicit none
 BEGIN_DOC
 ! Number of alpha and beta electrons
 END_DOC
 nelec_up = 5
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
