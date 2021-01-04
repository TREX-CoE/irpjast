BEGIN_PROVIDER [ double precision, elec_coord, (nelec, 3) ]
&BEGIN_PROVIDER [ double precision, nuc_coord, (nnuc, 3) ]
&BEGIN_PROVIDER [double precision, aord_vect, (naord + 1, typenuc)]
&BEGIN_PROVIDER [double precision, bord_vect, (nbord + 1)]
&BEGIN_PROVIDER [double precision, cord_vect, (dim_cord_vect, typenuc)]
 implicit none
 BEGIN_DOC
 ! Reads all input requested for Jatrow computation from an external procedure.
 ! Can be be used to interface with CHAMP
 END_DOC

 call irpinp(nelec, elec_coord, nnuc, typenuc, nuc_coord, &
	 naord, nbord, dim_cord_vect, aord_vect, bord_vect, cord_vect)

 TOUCH elec_coord
 TOUCH nuc_coord
 TOUCH aord_vect
 TOUCH bord_vect
 TOUCH cord_vect
 
END_PROVIDER
