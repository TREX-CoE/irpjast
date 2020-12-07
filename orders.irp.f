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

BEGIN_PROVIDER [double precision, aord_vect, (naord, typenuc)]
&BEGIN_PROVIDER [double precision, bord_vect, (nbord)]
&BEGIN_PROVIDER [double precision, cord_vect, (0:ncord , 0:ncord  , ncord , typenuc)]
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
