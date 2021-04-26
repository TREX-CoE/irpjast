! Generated from qmckl_dgemm.org

module qmckl_blas
  use :: iso_c_binding

  interface
     subroutine qmckl_dgemm(transa, transb, m, n, k, &
          alpha, A, lda, B, ldb, beta, C, ldc, tasks, ntasks) bind(C)
       use :: iso_c_binding
       implicit none
       character(kind=c_char  ), value :: transa, transb
       integer  (kind=c_int   ), value :: m, n, k, lda, ldb, ldc
       real     (kind=c_double), value :: alpha, beta
       real     (kind=c_double)        :: A(lda,*), B(ldb,*), C(ldc,*)
       integer  (kind=c_int64_t)       :: tasks(*)
       integer  (kind=c_int64_t)       :: ntasks
     end subroutine qmckl_dgemm
  end interface

  interface
     subroutine qmckl_tasks_run(tasks, ntasks) bind(C)
       use :: iso_c_binding
       implicit none
       integer (kind=c_int64_t), value :: ntasks
       integer (kind=c_int64_t)        :: tasks(ntasks)
     end subroutine qmckl_tasks_run
  end interface

end module qmckl_blas
