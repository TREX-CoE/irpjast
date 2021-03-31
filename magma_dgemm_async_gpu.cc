#include <stdio.h>
#include <stdlib.h>
#include "/p/software/juwelsbooster/stages/2020/software/magma/2.5.4-gcccoremkl-9.3.0-2020.2.254/include/magma_v2.h"

int magma_dgemm_async_gpu(double *A, double *B, double *C, int m, int n, int k, int lda, int ldb, int ldc)
{
 magma_init();
 // ... setup matrices on GPU:
 // m-by-k matrix dA,
 // k-by-n matrix dB,
 // m-by-n matrix dC.
 //int n = 100, m = 100, k = 100;
 //int lda = n, ldb = n, ldc = n;
 int ldda = magma_roundup( lda, 32 );
 int lddb = magma_roundup( ldb, 32 ); 
 int lddc = magma_roundup( ldc, 32 ); 
 int* ipiv = new int[ n ];

 double *dA, *dB, *dC;
 //double *A = new double[ lda*n ];
 //double *B = new double[ ldb*n ];
 //double *C = new double[ ldc*n ];
 
 magma_dmalloc_pinned( &dA, ldda*k ); 
 magma_dmalloc_pinned( &dB, lddb*n); 
 magma_dmalloc_pinned( &dC, lddc*n); 
 assert( dA != nullptr );
 assert( dB != nullptr );
 assert( dC != nullptr );

 int device;
 magma_queue_t queue;
 magma_getdevice( &device );
 magma_queue_create( device, &queue );

 // copy A, B to dA, dB
 magma_dsetmatrix( m, k,
 A, lda,
 dA, ldda, queue );

 magma_dsetmatrix( k, n,
 B, ldb,
 dB, lddb, queue );
 magma_queue_sync(queue);

 magma_dsetmatrix( m, n,
 C, ldc,
 dC, lddc, queue );
 magma_queue_sync(queue);

 // C = A B + C
 magma_dgemm( MagmaNoTrans,
 MagmaNoTrans, m, n, k,
 1.0, dA, ldda,
 dB, lddb,
 1.0, dC, lddc, queue );

 // ... use result in dC
 // copy result dC to C
 magma_dgetmatrix( m, n,
 dC, lddc,
 C, ldc, queue );

 // ... do concurrent work on CPU
 // wait for gemm to finish
 magma_queue_sync( queue );

 magma_queue_destroy( queue );
 // ... cleanup
 magma_free_pinned( dA ); 
 magma_free_pinned( dB );
 magma_free_pinned( dC );
 delete[] ipiv;

 magma_finalize();
 //printf("TEST5: Success DGEMM ASYNC CPU !\n");
 return 0;
}

extern "C" {
  void run_magma_dgemm_async_gpu_c(double *A, double *B, double *C, int m, int n, int k, int lda, int ldb, int ldc){
    //printf("Calling magma DGEMM -> (%d,%d,%d) (%d,%d,%d)\n",n,m,k,lda,ldb,ldc);
   int res = magma_dgemm_async_gpu(A,B,C,m,n,k,lda,ldb,ldc);
  }
}
