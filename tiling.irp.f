BEGIN_PROVIDER [ integer, tile_size ]
implicit none
BEGIN_DOC
! Tile size for tiling tables
END_DOC
tile_size = 16
END_PROVIDER

 BEGIN_PROVIDER [double precision, rescale_een_e_tiled, (tile_size, tile_size, ntiles_nelec, ntiles_nelec, 0:ncord)]
&BEGIN_PROVIDER [double precision, rescale_een_n_tiled, (tile_size, tile_size, ntiles_nelec, ntiles_nnuc, 0:ncord)]
&BEGIN_PROVIDER [double precision, rescale_een_e_deriv_e_tiled, (tile_size, nelec_16, 4, tile_size, nelec_16, 0:ncord)]
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
         idxj = ntiles_nelec*tile_size + jj
         do ii = 1, tile_size
           idxi = ntiles_nelec*tile_size + ii
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
             rescale_een_n_tiled(ii, aa, i, a, l) =  rescale_een_n(i, a, l)
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
             rescale_een_e_deriv_e_tiled(ii,i,k,jj,j,l) = rescale_een_e_deriv_e(i, k, j, l) 
          enddo
       enddo
       enddo
    enddo
    enddo
enddo

END_PROVIDER

