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

  interface 
    subroutine alloc(A, sze) bind(C)
      use :: iso_c_binding
      implicit none
      type(c_ptr) :: A
      integer(c_size_t), value :: sze
    end subroutine
  end interface

  interface 
    subroutine free(A) bind(C,name='starpu_free')
      use :: iso_c_binding
      implicit none
      type(c_ptr), value :: A
    end subroutine
  end interface

end module qmckl_blas

subroutine f_dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC) &
          bind(C, name='f_dgemm')
  use iso_c_binding
  implicit none
  character, intent(in), value :: TRANSA, TRANSB
  integer, intent(in), value   :: M,N,K,LDA,LDB,LDC
  double precision, intent(in), value :: ALPHA, BETA
  double precision, intent(in) :: A(LDA,*), B(LDB,*)
  double precision, intent(out) :: C(LDC,*)
  call dgemm(TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC)
end subroutine

