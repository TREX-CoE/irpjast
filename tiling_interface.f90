         module tiling_interface
         use :: ISO_C_BINDING
         interface
           subroutine run_magma_dgemm_async_gpu_c(A, B, C, n, m, k,lda,&
             ldb, ldc) bind(C)
           use :: ISO_C_BINDING
           implicit none
           integer, value :: n,m,k
           integer, value :: lda, ldb, ldc
           real(KIND=C_DOUBLE):: A(n,m)
           real(KIND=C_DOUBLE):: B(n,m)
           real(KIND=C_DOUBLE):: C(n,m)
           end subroutine run_magma_dgemm_async_gpu_c
         end interface
       end module

