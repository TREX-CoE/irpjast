BEGIN_PROVIDER [ integer, tile_size ]
implicit none
BEGIN_DOC
! Tile size for tiling tables
END_DOC
tile_size = 16
END_PROVIDER

 BEGIN_PROVIDER [double precision, rescale_een_e_tiled, (tile_size, tile_size, 0:ntiles_nelec, 0:ntiles_nelec, 0:ncord)]
&BEGIN_PROVIDER [double precision, rescale_een_n_tiled, (tile_size, tile_size, 0:ncord, 0:ntiles_nelec, 0:ntiles_nnuc)]
&BEGIN_PROVIDER [double precision, rescale_een_e_deriv_e_tiled, (tile_size, 4, tile_size, 0:ntiles_nelec, 0:ntiles_nelec, 0:ncord )]
 implicit none
 BEGIN_DOC
 ! R = exp(-kappa r) for electron-electron for $J_{een}$
 END_DOC
 integer :: i, j, k, l, a, ii, jj, aa
 integer :: idxi, idxj, idxa

 ! Fill up rescale_een_e_tiled
 do l = 0, ncord
    do j = 0, ntiles_nelec - 1
     do i = 0, ntiles_nelec - 1
       do jj = 1, tile_size
         idxj = j*tile_size + jj
         do ii = 1, tile_size
           idxi = i*tile_size + ii
           rescale_een_e_tiled(ii,jj,i,j,l) = rescale_een_e(idxi,idxj,l)
         enddo
       enddo
    enddo
    enddo
 enddo

 ! Fill up rescale_een_n_tiled
 do l = 0, ncord
    do a = 0, ntiles_nnuc - 1
       do i = 0, ntiles_nelec - 1
         do ii = 1, tile_size
           idxi = i*tile_size + ii
           do aa = 1, tile_size
             idxa = a*tile_size + aa
             rescale_een_n_tiled(ii, aa, l, i, a) =  rescale_een_n(idxi, idxa, l)
           enddo
         enddo
       enddo
    enddo
 enddo

 ! Fill up rescale_een_e_deriv_e_tiled
 do l = 0, ncord
    do j = 0, ntiles_nelec - 1
    do jj = 1, tile_size
       idxj = j*tile_size + jj
       do i = 0, ntiles_nelec - 1
       do ii = 1, tile_size
          idxi = i*tile_size + ii
          do k = 1, 4
             rescale_een_e_deriv_e_tiled(ii,k,jj,i,j,l) = rescale_een_e_deriv_e(idxi, k, idxj, l) 
          enddo
       enddo
       enddo
    enddo
    enddo
enddo

END_PROVIDER

 BEGIN_PROVIDER [ double precision,  tmp_c_tiled, (tile_size, tile_size,0:ncord, 0:ntiles_nelec, 0:ntiles_nnuc,0:ncord-1) ]
&BEGIN_PROVIDER [ double precision, dtmp_c_tiled, (tile_size, 4,tile_size,0:ncord, 0:ntiles_nelec, 0:ntiles_nnuc,0:ncord-1) ]
 !use tiling_interface
 implicit none
 BEGIN_DOC
 ! Calculate the intermediate buffers
 ! tmp_c:
 ! r_{ij}^k . R_{ja}^l -> tmp_c_{ia}^{kl}
 !
 ! dtmp_c:
 ! dr_{ij}^k . R_{ja}^l -> dtmp_c_{ia}^{kl}
 END_DOC
 integer :: k, i, j, a, l,m
 integer :: ii, jj, aa, kk, ll
 integer :: res

 ! r_{ij}^k . R_{ja}^l -> tmp_c_{ia}^{kl}
 do k=0,ncord-1
   do j = 0, ntiles_nelec - 1
     do i = 0, ntiles_nelec - 1
       do a = 0, ntiles_nnuc - 1
        do m = 0, ncord
              !DIR$ vector aligned
          do ii = 1, 16
              !DIR$ vector aligned
           do jj = 1, 16
              !DIR$ vector aligned
             do kk = 1, 16
                 tmp_c_tiled(ii,jj,m,j,a,k) = tmp_c_tiled(ii,jj,m,j,a,k) + &
                                          rescale_een_e_tiled(ii,kk,j,i,k)*&
                                          rescale_een_n_tiled(kk,jj,m,i,a)
            enddo
           enddo
         enddo
        enddo
   !call dgemm('N','N', tile_size, tile_size*(ncord+1), tile_size, 1.d0,           &
   !    rescale_een_e_tiled(1,1,j,i,k), size(rescale_een_e_tiled,1),                  &
   !    rescale_een_n_tiled(1,1,0,i,a), size(rescale_een_n_tiled,1), 1.d0,            &
   !    tmp_c_tiled(1,1,0,j,a,k), size(tmp_c_tiled,1))
   !call run_magma_dgemm_async_gpu_c(rescale_een_e_tiled(1,1,j,i,k),       &
   !                                rescale_een_n_tiled(1,1,0,i,a), &
   !                                tmp_c_tiled(1,1,0,j,a,k),       &
   !                                tile_size, tile_size*(ncord+1), &
   !                                tile_size,                      &
   !                                size(rescale_een_e_tiled,1),    &
   !                                size(rescale_een_n_tiled,1),    &
   !                                size(tmp_c_tiled,1))
       enddo
     enddo
   enddo
 enddo

 ! dr_{ij}^k . R_{ja}^l -> dtmp_c_{ia}^{kl}
 do k=0,ncord-1
   do i = 0, ntiles_nelec - 1
     do a = 0, ntiles_nnuc - 1
       do j = 0, ntiles_nelec - 1
        do m = 0, ncord
              !DIR$ vector aligned
          do ll = 1, 4
              !DIR$ vector aligned
           do ii = 1, 16
              !DIR$ vector aligned
            do jj = 1, 16
              !DIR$ vector aligned
             do kk = 1, 16
                 dtmp_c_tiled(ii,jj,ll,m,j,a,k) =       dtmp_c_tiled(ii,jj,ll,m,j,a,k) + &
                                          rescale_een_e_deriv_e_tiled(ii,kk,ll,j,i,k)*&
                                                  rescale_een_n_tiled(kk,jj,m,i,a)
   !call dgemm('N','N', 4*tile_size, tile_size*(ncord+1), tile_size, 1.d0,         &
   !    rescale_een_e_deriv_e_tiled(1,1,1,k,j,i), 4*size(rescale_een_e_deriv_e_tiled,1),&
   !    rescale_een_n_tiled(1,1,0,i,a), size(rescale_een_n_tiled,1), 1.d0,            &
   !    dtmp_c_tiled(1,1,1,0,j,a,k), 4*size(dtmp_c_tiled,1))
   !call run_magma_dgemm_async_gpu_c(rescale_een_e_deriv_e_tiled(1,1,1,k,j,i), &
   !                                rescale_een_n_tiled(1,1,0,i,a),            &
   !                                dtmp_c_tiled(1,1,1,0,j,a,k),               &
   !                                4*tile_size, tile_size*(ncord+1),          &
   !                                tile_size,                                 &
   !                                4*size(rescale_een_e_deriv_e_tiled,1),     &
   !                                size(rescale_een_n_tiled,1),               &
   !                                4*size(dtmp_c_tiled,1))
             enddo
            enddo
           enddo
         enddo
        enddo
       enddo
     enddo
   enddo
 enddo


END_PROVIDER

