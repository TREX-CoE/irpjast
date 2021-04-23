/* Generated from qmckl_dgemm.org */

#include <cblas.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>



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


#define MIN_SIZE 512
#include<stdio.h>

static void qmckl_dgemm_rec(struct dgemm_args args) {

//  printf("%5d %5d\n", args.m, args.n);

  if ( (args.m <= MIN_SIZE) || (args.n <= MIN_SIZE)) {
    #pragma omp task
    {
      printf("BLAS %5d %5d %5d\n", args.m, args.n, args.k);
      cblas_dgemm(CblasColMajor, args.transa, args.transb,
                  args.m, args.n, args.k, args.alpha,
                  args.A, args.lda, args.B, args.ldb,
                  args.beta, args.C, args.ldc);
    }
  } else {

    int m1 = args.m / 2;
    int m2 = args.m - m1;
    int n1 = args.n / 2;
    int n2 = args.n - n1;
      
    #pragma omp task
    {
      struct dgemm_args args_1 = args;
      args_1.m = m1;
      args_1.n = n1;
      qmckl_dgemm_rec(args_1);
    }
    
    #pragma omp task
    {
      // TODO: assuming 'N', 'N' here
      struct dgemm_args args_2 = args;
      args_2.B = args.B + args.ldb*n1;
      args_2.C = args.C + args.ldc*n1;
      args_2.m = m1;
      args_2.n = n2;
      qmckl_dgemm_rec(args_2);
    }
    
    #pragma omp task
    {
      struct dgemm_args args_3 = args;
      args_3.A = args.A + m1;
      args_3.C = args.C + m1;
      args_3.m = m2;
      args_3.n = n1;
      qmckl_dgemm_rec(args_3);
    }
    
    #pragma omp task
    {
      struct dgemm_args args_4 = args;
      args_4.A = args.A + m1;
      args_4.B = args.B + args.ldb*n1;
      args_4.C = args.C + m1 + args.ldc*n1;
      args_4.m = m2;
      args_4.n = n2;
      qmckl_dgemm_rec(args_4);
    }
  }

}

void qmckl_dgemm(char transa, char transb,
                 int m, int n, int k,
                 double alpha,
                 double* A, int lda,
                 double* B, int ldb,
                 double beta,
                 double* C, int ldc,
                 int64_t* result)
{
  struct dgemm_args* args = (struct dgemm_args*) malloc (sizeof(struct dgemm_args));
  assert (args != NULL);
  *result = (int64_t) args;

  args->alpha = alpha;
  args->beta  = beta ;
  args->A = A;
  args->B = B;
  args->C = C;
  args->m = m;
  args->n = n;
  args->k = k;
  args->lda = lda;
  args->ldb = ldb;
  args->ldc = ldc;

  if (transa == 'T' || transa == 't') {
    args->transa = CblasTrans;
  } else {
    args->transa = CblasNoTrans;
  }

  if (transa == 'T' || transa == 't') {
    args->transb = CblasTrans;
  } else {
    args->transb = CblasNoTrans;
  }

}


void qmckl_tasks_run(struct dgemm_args** gemms, int ngemms)
{
  #pragma omp parallel
  {
    #pragma omp single
    {
      for (int i=0 ; i<ngemms ; ++i)
      {
        qmckl_dgemm_rec(*(gemms[i]));
      }
    }
    #pragma omp taskwait
  }
}
