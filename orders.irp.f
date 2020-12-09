BEGIN_PROVIDER [integer, naord]
 implicit none
 BEGIN_DOC
 ! Expansion order for f_en
 END_DOC
 naord = 5
END_PROVIDER

BEGIN_PROVIDER [integer, nbord]
 implicit none
 BEGIN_DOC
 ! Expansion order for f_ee
 END_DOC
 nbord = 5
END_PROVIDER

BEGIN_PROVIDER [integer, ncord]
 implicit none
 BEGIN_DOC
 ! Expansion order for f_een
 END_DOC
 ncord = 5
END_PROVIDER
 
BEGIN_PROVIDER [integer, dim_cord_vect]
 implicit none
 BEGIN_DOC
 ! Recomputes the length of the unique C coefficients
 END_DOC
 integer :: k, p, l, lmax

 dim_cord_vect = 0

 do p = 2, ncord
    do k = 0, p - 1
       if ( k /= 0 ) then
          lmax = p - k
       else
          lmax = p - k - 2
       end if
       do l = iand(p - k, 1), lmax, 2
          dim_cord_vect = dim_cord_vect + 1
       end do
    end do
 end do

END_PROVIDER
 

BEGIN_PROVIDER [double precision, aord_vect, (naord + 1, typenuc)]
&BEGIN_PROVIDER [double precision, bord_vect, (nbord + 1)]
&BEGIN_PROVIDER [double precision, cord_vect, (dim_cord_vect, typenuc)]
 implicit none
 BEGIN_DOC
 ! Read Jastow coefficients from file
 END_DOC
 PROVIDE naord
 PROVIDE nbord
 PROVIDE ncord
 character(len=*), parameter :: FILE_NAME = "jast_coeffs.txt"
 integer :: i, fu, rc
 
 open(action='read', file=FILE_NAME, iostat=rc, newunit=fu)
 
 read(fu, *) aord_vect
 read(fu, *) bord_vect
 read(fu, *) cord_vect
 
 close(fu)

END_PROVIDER
