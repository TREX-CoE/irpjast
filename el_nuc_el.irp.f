BEGIN_PROVIDER [ double precision, factor_een ]
 implicit none
 BEGIN_DOC
 !
 END_DOC
 integer :: i,j,a,p,k,l,lmax,m
 double precision :: riam, rjam_cn, rial, rjal, rijk
 double precision :: cn

 factor_een = 0.0d0

 do p=2,ncord
   do k=0,p-1
     if (k /= 0) then
       lmax = p-k
     else
       lmax = p-k-2
     endif

     do l=0,lmax
       if ( iand(p-k-l,1) == 1) then
         cycle
       endif
       m = (p-k-l)/2

       do a=1, nnuc
         cn = cord_vect_lkp(l,k,p,typenuc_arr(a))
         do j=1, nelec
           rjal    = rescale_een_n(j,a,l)
           rjam_cn = rescale_een_n(j,a,m) * cn
           do i=1, j-1
             rial = rescale_een_n(i,a,l)
             riam = rescale_een_n(i,a,m)
             rijk = rescale_een_e(i,j,k)
             factor_een = factor_een + &
                            rijk * (rial+rjal) * riam * rjam_cn
           enddo
         enddo
       enddo

     enddo
   enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, factor_een_deriv_e, (4,nelec) ]
 implicit none
 BEGIN_DOC
 ! dimensions 1-3 : dx,dy,dz
 !
 ! Rdimension 4 : d2x + d2y + d2z
 END_DOC
 integer :: i,j,a,p,k,l,lmax,m
 double precision :: riam, rjam_cn, rial, rjal, rijk
 double precision, dimension(4) :: driam, drjam_cn, drial, drjal, drijk
 double precision :: cn, v1, v2, l1, l2, d1, d2

 factor_een_deriv_e = 0.0d0

 do p=2,ncord
   do k=0,p-1
     if (k /= 0) then
       lmax = p-k
     else
       lmax = p-k-2
     endif

     do l=0,lmax
       if ( iand(p-k-l,1) == 1) then
         cycle
       endif
       m = (p-k-l)/2

       do a=1, nnuc
         cn = cord_vect_lkp(l,k,p,typenuc_arr(a))
         do j=1, nelec
           rjal    = rescale_een_n(j,a,l)
           rjam_cn = rescale_een_n(j,a,m) * cn
           do ii=1,4
             drjal(ii)    = rescale_een_n_deriv_e(ii,j,a,l)
             drjam_cn(ii) = rescale_een_n_deriv_e(ii,j,a,m) * cn
           enddo
           factor_een_deriv_e(:,j) = 0.d0
           do i=1, nelec
             rial = rescale_een_n(i,a,l)
             riam = rescale_een_n(i,a,m)
             rijk = rescale_een_e(i,j,k)
             do ii=1,4
               drijk(ii) = rescale_een_e_deriv_e(ii,i,j,k)
             enddo
             l1 = 0.d0
             l2 = 0.d0
             x(1:3) = 0.d0
             x(4) = 2.d0
             do ii=1,4
               v1 = rijk * (rial+rjal)
               v2 = rjam_cn * riam

               d1 = drijk(ii) * (rial+rjal) + rijk * (rial+drjal(ii))
               d2 = drjam_cn(ii) * riam

               l1 = l1 + drijk(ii) * (rial+drjal(ii))
               l2 = l2 + drjam_cn(ii) * riam

               factor_een_deriv_e(ii,j) = factor_een_deriv_e(ii,j) + &
                  v1 * d2 + d1 * v2 + x(ii) * (l1 + l2)
             enddo
             factor_een_deriv_e(ii,j) = factor_een_deriv_e(ii,j) + &
                v1 * d2 + d1 * v2
           enddo
         enddo
       enddo

     enddo
   enddo
 enddo
 factor_een_deriv_e = 0.5d0 * factor_een_deriv_e

END_PROVIDER

