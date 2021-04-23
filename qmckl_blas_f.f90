! Generated from qmckl_dgemm.org

module qmckl_blas
  use :: iso_c_binding

  interface
     subroutine qmckl_dgemm(transa, transb, m, n, k, &
          alpha, A, lda, B, ldb, beta, C, ldc) bind(C)
       use :: iso_c_binding
       implicit none
       character(kind=c_char  ), value :: transa, transb
       integer  (kind=c_int   ), value :: m, n, k, lda, ldb, ldc
       real     (kind=c_double), value :: alpha, beta
       real     (kind=c_double)        :: A(lda,*), B(ldb,*), C(ldc,*)
     end subroutine qmckl_dgemm
  end interface

end module qmckl_blas
