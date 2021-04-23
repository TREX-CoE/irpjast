/* Generated from qmckl_dgemm.org */

#include <cblas.h>

struct dgemm_args {
  double alpha; 
  double beta;
  double* A;
  double* B;
  double* C;
  int m;
  int n;
  int k;
  int lda;
  int ldb;
  int ldc;
  CBLAS_LAYOUT transa;
  CBLAS_LAYOUT transb;
}; 
  

static void qmckl_dgemm_rec(struct dgemm_args args) {

  cblas_dgemm(CblasColMajor, args.transa, args.transb, 
              args.m, args.n, args.k, args.alpha,
              args.A, args.lda, args.B, args.ldb,
              args.beta, args.C, args.ldc);

}

void qmckl_dgemm(char transa, char transb,
                 int m, int n, int k,
                 double alpha, 
                 double* A, int lda,
                 double* B, int ldb,
                 double beta,
                 double* C, int ldc)
{
  struct dgemm_args args; 

  args.alpha = alpha;
  args.beta  = beta ;
  args.A = A;
  args.B = B;
  args.C = C;
  args.m = m;
  args.n = n;
  args.k = k;
  args.lda = lda;
  args.ldb = ldb;
  args.ldc = ldc;

  if (transa == 'T' || transa == 't') {
    args.transa = CblasTrans; 
  } else {
    args.transa = CblasNoTrans; 
  }

  CBLAS_LAYOUT tb;
  if (transa == 'T' || transa == 't') {
    args.transb = CblasTrans; 
  } else {
    args.transb = CblasNoTrans; 
  }

  qmckl_dgemm_rec(args);
}
