 BEGIN_PROVIDER [ double precision, factor_een_blas ]
&BEGIN_PROVIDER [ double precision, factor_een_deriv_e_blas, (4, nelec) ]
 implicit none
 BEGIN_DOC
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC

 integer                        :: i, j, a, p, k, l, lmax, m, n
 double precision               :: cn(ncord), accu
 double precision               :: f(nnuc,0:ncord-2,0:ncord-2)
 double precision               :: tmp_c(nelec,nnuc,0:ncord,0:ncord-1)
 double precision               :: dtmp_c(4,nelec,nnuc,0:ncord,0:ncord-1)

 factor_een_blas = 0.0d0
 factor_een_deriv_e_blas(1:4,1:nelec) = 0.0d0

 ! r_{ij}^k . R_{ja}^l -> tmp_c_{ia}^{kl}
 do k=0,ncord-1
   call dgemm('N','N', nelec, nnuc*(ncord+1), nelec, 1.d0,           &
       rescale_een_e(1,1,k), size(rescale_een_e,1),                  &
       rescale_een_n(1,1,0), size(rescale_een_n,1), 0.d0,            &
       tmp_c(1,1,0,k), size(tmp_c,1))
 enddo

 ! dr_{ij}^k . R_{ja}^l -> dtmp_c_{ia}^{kl}
 do k=0,ncord-1
   call dgemm('N','N', 4*nelec, nnuc*(ncord+1), nelec, 1.d0,         &
       rescale_een_e_deriv_e(1,1,1,k), 4*size(rescale_een_e_deriv_e,2),&
       rescale_een_n(1,1,0), size(rescale_een_n,1), 0.d0,            &
       dtmp_c(1,1,1,0,k), 4*size(dtmp_c,2))
 enddo


 do n = 1, dim_cord_vect

   l = lkpm_of_cindex(1,n)
   k = lkpm_of_cindex(2,n)
   p = lkpm_of_cindex(3,n)
   m = lkpm_of_cindex(4,n)

   do a = 1, nnuc
     cn(a) = cord_vect_full(n, a)
   enddo

   do a = 1, nnuc
     accu = 0.d0

     do j=1,nelec
       accu = accu + rescale_een_n(j,a,m) * tmp_c(j,a,m+l,k)

       factor_een_deriv_e_blas(1:4,j) = factor_een_deriv_e_blas(1:4,j) + (&
           tmp_c(j,a,m,k) * rescale_een_n_deriv_e(1:4,j,a,m+l) +     &
           dtmp_c(1:4,j,a,m,k) * rescale_een_n(j,a,m+l) +            &
           dtmp_c(1:4,j,a,m+l,k) * rescale_een_n(j,a,m) +            &
           tmp_c(j,a,m+l,k)*rescale_een_n_deriv_e(1:4,j,a,m)         &
           ) * cn(a)

       factor_een_deriv_e_blas(4,j) = factor_een_deriv_e_blas(4,j) + (&
           dtmp_c(1,j,a,m  ,k) * rescale_een_n_deriv_e(1,j,a,m+l) +  &
           dtmp_c(2,j,a,m  ,k) * rescale_een_n_deriv_e(2,j,a,m+l) +  &
           dtmp_c(3,j,a,m  ,k) * rescale_een_n_deriv_e(3,j,a,m+l) +  &
           dtmp_c(1,j,a,m+l,k) * rescale_een_n_deriv_e(1,j,a,m  ) +  &
           dtmp_c(2,j,a,m+l,k) * rescale_een_n_deriv_e(2,j,a,m  ) +  &
           dtmp_c(3,j,a,m+l,k) * rescale_een_n_deriv_e(3,j,a,m  )    &
           )*cn(a)*2.d0

     enddo
     factor_een_blas = factor_een_blas + accu * cn(a)

   enddo
 enddo

END_PROVIDER
