 BEGIN_PROVIDER [ double precision,  tmp_c, (nelec_16,nnuc_16,0:ncord,0:ncord-1) ]
&BEGIN_PROVIDER [ double precision, dtmp_c, (nelec_16,4,nnuc_16,0:ncord,0:ncord-1) ]
 implicit none
 BEGIN_DOC
 ! Calculate the intermediate buffers
 ! tmp_c:
 ! r_{ij}^k . R_{ja}^l -> tmp_c_{ia}^{kl}
 !
 ! dtmp_c:
 ! dr_{ij}^k . R_{ja}^l -> dtmp_c_{ia}^{kl}
 END_DOC
 integer :: k

 ! r_{ij}^k . R_{ja}^l -> tmp_c_{ia}^{kl}
 do k=0,ncord-1
   call dgemm('N','N', nelec_16, nnuc_16*(ncord+1), nelec_16, 1.d0,           &
       rescale_een_e(1,1,k), size(rescale_een_e,1),                  &
       rescale_een_n(1,1,0), size(rescale_een_n,1), 0.d0,            &
       tmp_c(1,1,0,k), size(tmp_c,1))
 enddo

 ! dr_{ij}^k . R_{ja}^l -> dtmp_c_{ia}^{kl}
 do k=0,ncord-1
   call dgemm('N','N', 4*nelec_16, nnuc_16*(ncord+1), nelec, 1.d0,         &
       rescale_een_e_deriv_e(1,1,1,k), 4*size(rescale_een_e_deriv_e,1),&
       rescale_een_n(1,1,0), size(rescale_een_n,1), 0.d0,            &
       dtmp_c(1,1,1,0,k), 4*size(dtmp_c,1))
 enddo


END_PROVIDER


 BEGIN_PROVIDER [ double precision, factor_een_blas ]
&BEGIN_PROVIDER [ double precision, factor_een_deriv_e_blas, (nelec_16,4) ]
 implicit none
 BEGIN_DOC
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC

 integer                        :: i, j, a, p, k, l, lmax, m, n, ii
 double precision               :: accu, cn
! double precision,dimension(:),allocatable :: cn

 factor_een_blas = 0.0d0
 factor_een_deriv_e_blas(1:nelec_16,1:4) = 0.0d0

 do n = 1, dim_cord_vect

   l = lkpm_of_cindex(1,n)
   k = lkpm_of_cindex(2,n)
   p = lkpm_of_cindex(3,n)
   m = lkpm_of_cindex(4,n)

   do a = 1, nnuc
     cn = cord_vect_full(n, a)
     if (cn == 0.d0) cycle

     accu = 0.d0
     do j=1,nelec
       accu = accu + rescale_een_n(j,a,m) * tmp_c(j,a,m+l,k)
     enddo
     factor_een_blas = factor_een_blas + accu * cn

     do ii=1,4
       do j=1,nelec
         factor_een_deriv_e_blas(j,ii) = factor_een_deriv_e_blas(j,ii) + (&
           tmp_c(j,a,m,k) * rescale_een_n_deriv_e(j,ii,a,m+l) +     &
           dtmp_c(j,ii,a,m,k)   * rescale_een_n(j,a,m+l) +            &
           dtmp_c(j,ii,a,m+l,k) * rescale_een_n(j,a,m) +            &
           tmp_c(j,a,m+l,k)*rescale_een_n_deriv_e(j,ii,a,m)         &
           ) * cn
       enddo
     enddo

     cn = cn+cn
     do j=1,nelec
       factor_een_deriv_e_blas(j,4) = factor_een_deriv_e_blas(j,4) + (&
           dtmp_c(j,1,a,m  ,k) * rescale_een_n_deriv_e(j,1,a,m+l) +  &
           dtmp_c(j,2,a,m  ,k) * rescale_een_n_deriv_e(j,2,a,m+l) +  &
           dtmp_c(j,3,a,m  ,k) * rescale_een_n_deriv_e(j,3,a,m+l) +  &
           dtmp_c(j,1,a,m+l,k) * rescale_een_n_deriv_e(j,1,a,m  ) +  &
           dtmp_c(j,2,a,m+l,k) * rescale_een_n_deriv_e(j,2,a,m  ) +  &
           dtmp_c(j,3,a,m+l,k) * rescale_een_n_deriv_e(j,3,a,m  )    &
           )*cn
     enddo

   enddo
 enddo

END_PROVIDER
