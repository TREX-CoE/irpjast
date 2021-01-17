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
 call jast_elec_champ(nelec, elec_coord)
END_PROVIDER

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
 ! Type of the nuclei
 END_DOC
 typenuc = 1
 typenuc_arr = (/1, 1/)
END_PROVIDER

BEGIN_PROVIDER [ double precision, nuc_coord, (nnuc, 3) ]
 implicit none
 BEGIN_DOC
 ! Nuclei coordinates
 END_DOC
 call jast_nuc_champ(nnuc, nuc_coord)
END_PROVIDER

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
 
BEGIN_PROVIDER [integer, ndim_cord_vect]
 implicit none
 BEGIN_DOC
 ! Recomputes the length of the unique C coefficients
 END_DOC
 integer :: k, p, l, lmax

 ndim_cord_vect = 0

 do p = 2, ncord
    do k = 0, p - 1
       if ( k /= 0 ) then
          lmax = p - k
       else
          lmax = p - k - 2
       end if
       do l = iand(p - k, 1), lmax, 2
          ndim_cord_vect += 1
       end do
    end do
 end do

END_PROVIDER
 
BEGIN_PROVIDER [double precision, aord_vect, (naord + 1, typenuc)]
&BEGIN_PROVIDER [double precision, bord_vect, (nbord + 1)]
&BEGIN_PROVIDER [double precision, cord_vect, (ndim_cord_vect, typenuc)]
 implicit none
 BEGIN_DOC
 ! Read Jastow coefficients from file
 END_DOC
 call jast_pars_champ(naord, typenuc, aord_vect, nbord, bord_vect, ndim_cord_vect, cord_vect)
END_PROVIDER

BEGIN_PROVIDER [ double precision, cord_vect_lkp, (0:ncord-1, 0:ncord-1, 2:ncord, typenuc) ]
 implicit none
 BEGIN_DOC
 ! Creates c-tensor with right order of the indexes p, k, l
 END_DOC
 integer :: alpha, l, k, p, lmax, cindex

 cord_vect_lkp = 0.0d0
 cindex = 0
 do alpha = 1, typenuc
   do p = 2, ncord
     do k = p - 1, 0, -1
       if ( k /= 0 ) then
         lmax = p - k
       else
         lmax = p - k - 2
       end if
       do l = lmax, 0, -1
         if (iand(p - k - l, 1) == 1) cycle
         cindex = cindex + 1
         cord_vect_lkp(l, k, p, alpha) = cord_vect(cindex, alpha)
       end do
     end do
   end do
 end do

END_PROVIDER
