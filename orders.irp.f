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
 do i = 1, naord 
    call random_number(aord_vect)
 end do
END_PROVIDER

BEGIN_PROVIDER [double precision, bord_vect, (nbord)]
 implicit none
 BEGIN_DOC
 ! Vector of the `b' coefficients
 END_DOC
 do i = 1, nbord 
    call random_number(bord_vect)
 end do
END_PROVIDER

BEGIN_PROVIDER [double precision, cord_vect, (ncord)]
 implicit none
 BEGIN_DOC
 ! Vector of the `b' coefficients
 END_DOC
 do i = 1, ncord 
    call random_number(cord_vect)
 end do
END_PROVIDER
