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

BEGIN_PROVIDER [double precision, aord_vect, (naord)]
 implicit none
 BEGIN_DOC
 ! Vector of the `a' coefficients
 END_DOC
 integer :: i
 call random_number(aord_vect)
 aord_vect = aord_vect*.1e-2
END_PROVIDER

BEGIN_PROVIDER [double precision, bord_vect, (nbord)]
 implicit none
 BEGIN_DOC
 ! Vector of the `b' coefficients
 END_DOC
 integer :: i
 call random_number(bord_vect)
 bord_vect = bord_vect*.1e-6
END_PROVIDER

BEGIN_PROVIDER [double precision, cord_vect, (0:ncord,0:ncord,ncord,nnuc)]
 implicit none
 BEGIN_DOC
 ! Vector of the `c' coefficients
 END_DOC
 call random_number(cord_vect)
 cord_vect = cord_vect*.1e-4
END_PROVIDER
