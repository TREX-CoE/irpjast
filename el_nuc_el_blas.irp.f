BEGIN_PROVIDER [ double precision, factor_een_blas ]
  implicit none
  BEGIN_DOC
  ! ElectronE-electron-nuclei contribution to Jastrow factor
  !
  ! 4124.84239750000
  END_DOC
  integer                        :: i, j, a, p, k, l, lmax, m, n
  double precision               :: cn(nnuc), accu
  double precision               :: f(nnuc,0:ncord-2,0:ncord-2)
  double precision               :: tmp_c(nelec,nnuc,0:ncord,0:ncord-1)

 factor_een_blas = 0.0d0

 ! r_{ij}^k . R_{ja}^l -> tmp_c_{ia}^{kl}
 do k=0,ncord-1
   call dgemm('N','N', nelec, nnuc*(ncord+1), nelec, 1.d0,           &
       rescale_een_e(1,1,k), size(rescale_een_e,1),                  &
       rescale_een_n(1,1,0), size(rescale_een_n,1), 0.d0,            &
       tmp_c(1,1,0,k), size(tmp_c,1))
 enddo

 do p = 2, ncord
   do k = 0, p - 1
     m = p-k
     if (k > 0) then
       lmax = m
     else
       lmax = m - 2
     endif

     n = shiftr(m,1)
     do l = iand(m, 1), lmax, 2

       do a = 1, nnuc
         cn(a) = cord_vect_lkp(l, k, p, typenuc_arr(a))
       enddo

       do a = 1, nnuc
         accu = 0.d0
         do i=1,nelec
           accu = accu + rescale_een_n(i,a,n) * tmp_c(i,a,n+l,k)
         enddo
         factor_een_blas = factor_een_blas + accu * cn(a)
       enddo
       n = n-1

     enddo
   enddo
 enddo

END_PROVIDER

BEGIN_PROVIDER [ double precision, factor_een_deriv_e_blas, (4, nelec) ]
 implicit none
 BEGIN_DOC
 ! Dimensions 1-3 : dx, dy, dz
 ! Dimension 4 : d2x + d2y + d2z
 END_DOC

 integer                        :: i, j, a, p, k, l, lmax, m, n
 double precision               :: cn(ncord), x
 double precision               :: f(nnuc,0:ncord-2,0:ncord-2)
 double precision               :: tmp_c(nelec,nnuc,0:ncord,0:ncord-1)
 double precision               :: dtmp_c(4,nelec,nnuc,0:ncord,0:ncord-1)

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


 do p = 2, ncord
   do k = 0, p - 1
      m = p-k
      if (k > 0) then
        lmax = m
      else
        lmax = m - 2
      endif

     n = shiftr(m,1)
     do l = iand(m, 1), lmax, 2

       do a = 1, nnuc
         cn(a) = cord_vect_lkp(l, k, p, typenuc_arr(a))
       enddo

       do a = 1, nnuc
         do j=1,nelec

           factor_een_deriv_e_blas(1:4,j) = factor_een_deriv_e_blas(1:4,j) + (&
               tmp_c(j,a,n,k) * rescale_een_n_deriv_e(1:4,j,a,n+l) + &
               dtmp_c(1:4,j,a,n,k) * rescale_een_n(j,a,n+l) +        &
               dtmp_c(1:4,j,a,n+l,k) * rescale_een_n(j,a,n) +        &
               tmp_c(j,a,n+l,k)*rescale_een_n_deriv_e(1:4,j,a,n)     &
               ) * cn(a)

           factor_een_deriv_e_blas(4,j) = factor_een_deriv_e_blas(4,j) + 2.d0*(&
               dtmp_c(1,j,a,n  ,k) * rescale_een_n_deriv_e(1,j,a,n+l) +&
               dtmp_c(2,j,a,n  ,k) * rescale_een_n_deriv_e(2,j,a,n+l) +&
               dtmp_c(3,j,a,n  ,k) * rescale_een_n_deriv_e(3,j,a,n+l) +&
               dtmp_c(1,j,a,n+l,k) * rescale_een_n_deriv_e(1,j,a,n  ) +&
               dtmp_c(2,j,a,n+l,k) * rescale_een_n_deriv_e(2,j,a,n  ) +&
               dtmp_c(3,j,a,n+l,k) * rescale_een_n_deriv_e(3,j,a,n  )&
               )*cn(a)

         enddo
       enddo
       n = n-1

     enddo
   enddo
 enddo
END_PROVIDER
