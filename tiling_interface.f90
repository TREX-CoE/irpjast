#if defined(__AVX512VL__)
#define VL 8
#elif defined(__AVX__) || defined(__AVX2__)
#define VL 4
#elif defined(__SSE2__) || defined(__SSE3__) || defined(__SSE4_1__) || defined(__SSE4_2__)
#define VL 2
#else
#define VL 1
#endif
         module tiling_interface
         use :: ISO_C_BINDING
         interface
           subroutine run_magma_dgemm_async_gpu_c(A, B, C, m, n, k,lda,&
             ldb, ldc) bind(C)
           use :: ISO_C_BINDING
           implicit none
           integer(KIND=C_INT), value :: n,m,k
           integer(KIND=C_INT), value :: lda, ldb, ldc
           real(KIND=C_DOUBLE):: A(m,k)
           real(KIND=C_DOUBLE):: B(k,n)
           real(KIND=C_DOUBLE):: C(m,n)
           end subroutine run_magma_dgemm_async_gpu_c
         end interface
       end module

       module simd
             use :: iso_c_binding
             type, public :: simd_real8
             real(c_double) :: x(0:VL-1)
             end type simd_real8
             type, public :: simd_int4
             integer(c_int) :: x(0:VL-1)
             end type simd_int4
             type, public :: simd_mask8
             logical :: x(0:VL-1)
             end type simd_mask8
       end module simd

