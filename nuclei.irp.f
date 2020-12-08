BEGIN_PROVIDER [ integer, nnuc ]
 implicit none
 BEGIN_DOC
 ! Number of nuclei
 END_DOC
 nnuc = 2
END_PROVIDER


BEGIN_PROVIDER [ integer, typenuc ]
&BEGIN_PROVIDER [integer, typenuc_arr, (nnuc)]
 implicit none
 BEGIN_DOC
 ! Number of nuclei
 END_DOC
 typenuc = 1
 typenuc_arr = (/1, 1/)
END_PROVIDER


BEGIN_PROVIDER [ double precision, nuc_coord, (nnuc, 3) ]
 implicit none
 BEGIN_DOC
 ! Nuclei coordinates
 END_DOC
 character(len=*), parameter :: FILE_NAME = "geometry.txt"
 integer :: fu, rc, i
 
 open(action='read', file=FILE_NAME, iostat=rc, newunit=fu)

 do i = 1, nnuc
    read(fu, *) nuc_coord(i, :)
 end do
 
 close(fu)

END_PROVIDER

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
 integer :: i, j, p, q
 double precision :: pow_ser, x

 factor_en = 0.0d0

 do j = 1 , nnuc
    do i = 1, nelec
       x = rescale_en(i, j) 
       pow_ser = 0.0d0
       do p = 2, naord
          x = x * rescale_en(i, j) 
          pow_ser = pow_ser + aord_vect(p + 1, typenuc_arr(j)) * x
       end do
       factor_en = factor_en + aord_vect(1, typenuc_arr(j)) * rescale_en(i, j) &
            / (1 + aord_vect(2, typenuc_arr(j)) * rescale_en(i, j)) + pow_ser
    end do
 end do

END_PROVIDER
