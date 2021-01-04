subroutine irpinp(nelec, elec_coord, nnuc, typenuc, nuc_coord, &
       	naord, nbord, ncord, aord_vect, bord_vect, cord_vect)
 ! This is a ghost subroutine to interop with CHAMP variables
 implicit none
 integer, intent(in) :: nelec
 integer, intent(in) :: nnuc
 integer, intent(in) :: typenuc
 integer, intent(in) :: naord
 integer, intent(in) :: nbord
 integer, intent(in) :: ncord
 double precision, dimension(nelec, 3), intent(inout) :: elec_coord
 double precision, dimension(nnuc, 3), intent(inout) :: nuc_coord
 double precision, dimension(naord + 1, typenuc), intent(inout) :: aord_vect
 double precision, dimension(nbord + 1), intent(inout) :: bord_vect
 double precision, dimension(ncord, typenuc), intent(inout) :: cord_vect
end subroutine irpinp
