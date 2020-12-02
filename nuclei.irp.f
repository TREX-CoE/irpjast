BEGIN_PROVIDER [ integer, nnuc ]
 implicit none
 BEGIN_DOC
 ! Number of nuclei
 END_DOC
 nnuc = 10
END_PROVIDER


BEGIN_PROVIDER [ double precision, nuc_coord, (nnuc, 3) ]
 implicit none
 BEGIN_DOC
 ! Nuclei coordinates
 END_DOC
 integer :: i, j
 PROVIDE seed
 do j = 1 , 3
   do i = 1, nnuc
     call random_number(nuc_coord(i, j))
   enddo
 enddo
 FREE seed

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
 integer :: i, j, p
 double precision :: pow_ser = 0.0d0
 factor_en = 0.0d0

 do j = 1 , nnuc
    do i = 1, nnuc
       do p = 2, naord
          pow_ser = pow_ser + aord_vect(p) * rescale_en(i, j) ** p
       end do
       factor_en = factor_en + aord_vect(1) * rescale_en(i, j) &
            / (1 + aord_vect(2) * rescale_en(i, j)) + pow_ser
    end do
 end do

 factor_en = 0.5d0 * factor_en
END_PROVIDER
