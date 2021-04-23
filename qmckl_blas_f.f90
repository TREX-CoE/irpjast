! Generated from qmckl_dgemm.org

module qmckl_blas
  use :: iso_c_binding

  interface
     subroutine qmckl_dgemm(transa, transb, m, n, k, &
          alpha, A, lda, B, ldb, beta, C, ldc, res) bind(C)
       use :: iso_c_binding
       implicit none
       character(kind=c_char  ), value :: transa, transb
       integer  (kind=c_int   ), value :: m, n, k, lda, ldb, ldc
       real     (kind=c_double), value :: alpha, beta
       real     (kind=c_double)        :: A(lda,*), B(ldb,*), C(ldc,*)
       integer  (kind=c_int64_t)       :: res
     end subroutine qmckl_dgemm
  end interface

  interface
     subroutine qmckl_tasks_run(gemms, ngemms) bind(C)
       use :: iso_c_binding
       implicit none
       integer (kind=c_int32_t), value :: ngemms
       integer (kind=c_int64_t)        :: gemms(ngemms)
     end subroutine qmckl_tasks_run
  end interface

end module qmckl_blas
