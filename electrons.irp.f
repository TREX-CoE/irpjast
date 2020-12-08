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
 enddo
END_PROVIDER

BEGIN_PROVIDER [double precision, factor_ee]
 implicit none
 BEGIN_DOC
 ! Electron-electron contribution to Jastrow factor
 END_DOC
 integer :: i, j, p
 double precision :: pow_ser, x, b_one

 factor_ee = 0.0d0

 do j = 1, nelec
    do i = 1, nelec
       x = rescale_ee(i, j) 
       pow_ser = 0.0d0
       do p = 2, nbord
          x = x * rescale_ee(i, j)
          pow_ser = pow_ser + bord_vect(p + 1) * x
       end do

       if (i <= nelec_up .or. j >= nelec_up) then
           b_one = bord_vect(1) * 0.5d0
       else
	   b_one = bord_vect(1)
       end if

       factor_ee = factor_ee + b_one * rescale_ee(i, j) &
            / (1.0d0 + bord_vect(2) * rescale_ee(i, j)) + pow_ser
    end do
 end do

 factor_ee = 0.5d0 * factor_ee

END_PROVIDER
