
 BEGIN_PROVIDER [ double precision, factor_een_blas_tiled_simd ]
&BEGIN_PROVIDER [ double precision, factor_een_deriv_e_blas_tiled_simd, (nelec_16,4) ]
 implicit none
 BEGIN_DOC
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC

 integer                        :: i, j, a, p, k, l, lmax, m, n, ii, jj
 integer                        :: idxi, idxj, idxa, aa
 double precision               :: accu, cn
! double precision,dimension(:),allocatable :: cn

 factor_een_blas_tiled_simd = 0.0d0
 factor_een_deriv_e_blas_tiled_simd(1:nelec_16,1:4) = 0.0d0

 do n = 1, dim_cord_vect

   l = lkpm_of_cindex(1,n)
   k = lkpm_of_cindex(2,n)
   p = lkpm_of_cindex(3,n)
   m = lkpm_of_cindex(4,n)

   do a = 0, ntiles_nnuc - 1
    do aa = 1, tile_size
      idxa = a*tile_size + aa

      cn = cord_vect_full(n, idxa)
      if (cn == 0.d0) cycle

      accu = 0.d0
      do j = 0, ntiles_nelec - 1
        do jj = 1, tile_size
          idxj = j*tile_size + jj
          accu = accu + rescale_een_n_tiled(jj,aa,m,j,a) * tmp_c_tiled_simd(jj,aa,m+l,j,a,k)
        enddo
      enddo
      factor_een_blas_tiled_simd = factor_een_blas_tiled_simd + accu * cn

      do j = 0, ntiles_nelec - 1
        do ii=1,4
        do jj = 1, tile_size
          idxj = j*tile_size + jj
          factor_een_deriv_e_blas_tiled_simd(idxj,ii) = factor_een_deriv_e_blas_tiled_simd(idxj,ii) + (&
            tmp_c_tiled_simd(jj,aa,m,j,a,k) * rescale_een_n_deriv_e(idxj,ii,idxa,m+l) +     &
            dtmp_c_tiled_simd(jj,aa,ii,m,j,a,k)   * rescale_een_n_tiled(jj,aa,m+l,j,a) +            &
            dtmp_c_tiled_simd(jj,aa,ii,m+l,j,a,k) * rescale_een_n_tiled(jj,aa,m  ,j,a) +            &
            tmp_c_tiled_simd(jj,aa,m+l,j,a,k)*rescale_een_n_deriv_e(idxj,ii,idxa,m)         &
            ) * cn
        enddo
        enddo
      enddo

      cn = cn+cn
      do j = 0, ntiles_nelec - 1
        do jj = 1, tile_size
          idxj = j*tile_size + jj
        factor_een_deriv_e_blas_tiled_simd(idxj,4) = factor_een_deriv_e_blas_tiled_simd(idxj,4) + (&
            dtmp_c_tiled_simd(jj,aa,1,m  ,j,a,k) * rescale_een_n_deriv_e(idxj,1,idxa,m+l) +  &
            dtmp_c_tiled_simd(jj,aa,2,m  ,j,a,k) * rescale_een_n_deriv_e(idxj,2,idxa,m+l) +  &
            dtmp_c_tiled_simd(jj,aa,3,m  ,j,a,k) * rescale_een_n_deriv_e(idxj,3,idxa,m+l) +  &
            dtmp_c_tiled_simd(jj,aa,1,m+l,j,a,k) * rescale_een_n_deriv_e(idxj,1,idxa,m  ) +  &
            dtmp_c_tiled_simd(jj,aa,2,m+l,j,a,k) * rescale_een_n_deriv_e(idxj,2,idxa,m  ) +  &
            dtmp_c_tiled_simd(jj,aa,3,m+l,j,a,k) * rescale_een_n_deriv_e(idxj,3,idxa,m  )    &
            )*cn
        enddo
      enddo

    enddo ! aa
   enddo ! a
 enddo ! n

END_PROVIDER
