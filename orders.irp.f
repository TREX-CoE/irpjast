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
&BEGIN_PROVIDER [double precision, bord_vect, (nbord)]
&BEGIN_PROVIDER [double precision, cord_vect, (ncord * ncord * ncord * nnuc)]
implicit none
PROVIDE naord
PROVIDE nbord
PROVIDE ncord
BEGIN_DOC
! Read Jastow coefficients from file (NEEDS OPTIMIZATION!)
END_DOC
character(len=*), parameter :: FILE_NAME = "orders_inp"
integer :: i, fu, rc
double precision, dimension(naord + nbord + ncord * ncord * ncord * nnuc) :: allord_vect

open(action='read', file=FILE_NAME, iostat=rc, newunit=fu)

do i = 1, naord + nbord + ncord * ncord * ncord * nnuc
   read(fu, *) allord_vect(i)
end do

aord_vect = allord_vect(1:naord)
bord_vect = allord_vect(naord + 1: naord + nbord)
cord_vect = allord_vect(naord + nbord + 1:)

close(fu)

END_PROVIDER

! BEGIN_PROVIDER [double precision, aord_vect, (naord)]
!  implicit none
!  BEGIN_DOC
!  ! Vector of the `a' coefficients
!  END_DOC
!  integer :: i
!  PROVIDE seed
!  call random_number(aord_vect)
!  aord_vect = aord_vect*.1d-2
!  FREE seed
! END_PROVIDER
! 
! BEGIN_PROVIDER [double precision, bord_vect, (nbord)]
!  implicit none
!  BEGIN_DOC
!  ! Vector of the `b' coefficients
!  END_DOC
!  integer :: i
!  PROVIDE seed
!  call random_number(bord_vect)
!  bord_vect = bord_vect*.1d-6
!  FREE seed
! END_PROVIDER
! 
BEGIN_PROVIDER [double precision, cord_vect_0, (0:ncord,0:ncord,ncord,nnuc)]
 implicit none
 BEGIN_DOC
 ! Vector of the `c' coefficients
 END_DOC
 PROVIDE seed
 call random_number(cord_vect_0)
 cord_vect_0 = cord_vect_0 * .1d-4
 FREE seed
END_PROVIDER
