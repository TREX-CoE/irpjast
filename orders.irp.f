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
   do k = p - 1, 0, -1
     if ( k /= 0 ) then
       lmax = p - k
     else
       lmax = p - k - 2
     end if
     do l = lmax, 0, -1
       if (iand(p - k - l, 1) == 1) cycle
       dim_cord_vect = dim_cord_vect + 1
       end do
    end do
 end do

END_PROVIDER


BEGIN_PROVIDER [ integer, lkpm_of_cindex, (4,dim_cord_vect) ]
  implicit none
  BEGIN_DOC
! Transform l,k,p into a consecutive index
  END_DOC
  integer                        :: p,k,l,lmax,m
  integer                        :: kk
  kk=0
  do p = 2, ncord
    do k = p - 1, 0, -1
      if ( k /= 0 ) then
        lmax = p - k
      else
        lmax = p - k - 2
      end if
      do l = lmax, 0, -1
        if (iand(p - k - l, 1) == 1) cycle
        m = (p - k - l) / 2
        kk = kk+1
        lkpm_of_cindex(1,kk) = l
        lkpm_of_cindex(2,kk) = k
        lkpm_of_cindex(3,kk) = p
        lkpm_of_cindex(4,kk) = m
      enddo
    enddo
  enddo

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

BEGIN_PROVIDER [ double precision, cord_vect_full, (dim_cord_vect, nnuc) ]
 implicit none
 BEGIN_DOC
 ! cord_vect for all atoms
 END_DOC
 integer :: a
 do a=1,nnuc
   cord_vect_full(1:dim_cord_vect,a) = cord_vect(1:dim_cord_vect,typenuc_arr(a))
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, cord_vect_lkp, (0:ncord-1, 0:ncord-1, 2:ncord, typenuc) ]
 implicit none
 BEGIN_DOC
 ! Creates c-tensor with right order of the indexes p, k, l
 END_DOC
 integer :: alpha, l, k, p, lmax, cindex

 cord_vect_lkp = 0.0d0
 do alpha = 1, typenuc
   cindex = 0
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
