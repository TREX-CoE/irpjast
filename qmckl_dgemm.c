/* Generated from qmckl_dgemm.org */

#include <starpu.h>
#include <chameleon.h>

#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#define GPU_ENABLED 1

#ifdef GPU_ENABLED

#include <cuda.h>
#include <starpu_cublas_v2.h>

#endif

void f_dgemm(const char transa, const char transb, const int m, const int n, const int k,
             const double alpha, const double* A, const int lda, const double* B,
	     const int ldb, const double beta, double* C, const int ldc);



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
  char transa;
  char transb;
};


void dgemm_codelet_cpu(void *buffers[], void* cl_arg)
{
  struct dgemm_args *args = cl_arg;
  double* A = (double*) STARPU_MATRIX_GET_PTR(buffers[0]);
  double* B = (double*) STARPU_MATRIX_GET_PTR(buffers[1]);
  double* C = (double*) STARPU_MATRIX_GET_PTR(buffers[2]);

  int lda = STARPU_MATRIX_GET_LD(buffers[0]);
  int ldb = STARPU_MATRIX_GET_LD(buffers[1]);
  int ldc = STARPU_MATRIX_GET_LD(buffers[2]);

  int m = STARPU_MATRIX_GET_NX(buffers[2]);
  int n = STARPU_MATRIX_GET_NY(buffers[2]);
  int k = STARPU_MATRIX_GET_NY(buffers[0]);

  f_dgemm(args->transa, args->transb,
          m, n, k, args->alpha,
          A, lda, B, ldb, args->beta, C, ldc);

  free(args);
}

#ifdef GPU_ENABLED
void dgemm_codelet_gpu(void *buffers[], void* cl_arg)
{
  struct dgemm_args *args = cl_arg;
  double* A = (double*) STARPU_MATRIX_GET_PTR(buffers[0]);
  double* B = (double*) STARPU_MATRIX_GET_PTR(buffers[1]);
  double* C = (double*) STARPU_MATRIX_GET_PTR(buffers[2]);

  int lda = STARPU_MATRIX_GET_LD(buffers[0]);
  int ldb = STARPU_MATRIX_GET_LD(buffers[1]);
  int ldc = STARPU_MATRIX_GET_LD(buffers[2]);

  int m = STARPU_MATRIX_GET_NX(buffers[2]);
  int n = STARPU_MATRIX_GET_NY(buffers[2]);
  int k = STARPU_MATRIX_GET_NY(buffers[0]);

  char transa = (args->transa == 'T' || args->transa == 'T')? CUBLAS_OP_T : CUBLAS_OP_N;
  char transb = (args->transb == 'T' || args->transb == 'T')? CUBLAS_OP_T : CUBLAS_OP_N;

  cublasStatus_t status = cublasDgemm(starpu_cublas_get_local_handle(),
          transa, transb, m, n, k, &(args->alpha),
          A, lda, B, ldb, &(args->beta), C, ldc);

  if (status != CUBLAS_STATUS_SUCCESS)
                STARPU_CUBLAS_REPORT_ERROR(status);

  free(args);
}
#endif

static struct starpu_perfmodel perf_model =
{
    .type = STARPU_HISTORY_BASED,
    .symbol = "my_perfmodel",
};

struct starpu_codelet dgemm_cl =
  {
   .where = STARPU_CPU | STARPU_CUDA,
   .cpu_funcs = { dgemm_codelet_cpu },
   .cpu_funcs_name = { "dgemm_codelet_cpu" },
#ifdef GPU_ENABLED
   .cuda_funcs = { dgemm_codelet_gpu },
   .cuda_flags = {STARPU_CUDA_ASYNC},
#endif
   .nbuffers = 3,
   .max_parallelism = 1,
   .modes = {STARPU_R, STARPU_R, STARPU_RW},
   .model = &perf_model,
  };

#include<stdio.h>

static struct dgemm_args* qmckl_dgemm_to_struct(char transa, char transb,
                 int m, int n, int k,
                 double alpha,
                 double* A, int lda,
                 double* B, int ldb,
                 double beta,
                 double* C, int ldc)
{
  struct dgemm_args* args = (struct dgemm_args*) malloc (sizeof(struct dgemm_args));
  assert (args != NULL);


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
  args->transa = transa;
  args->transb = transb;
  return args;

}

#define MIN_SIZE 4096*4096
static void qmckl_dgemm_rec(struct dgemm_args args, int64_t* tasks, int64_t* ntasks)
{


  if ( args.m * args.n <= MIN_SIZE) {

      struct dgemm_args* args_new = (struct dgemm_args*) malloc (sizeof(struct dgemm_args));
      memcpy(args_new, &args, sizeof(args));
      tasks[*ntasks] = (int64_t) args_new;
      *ntasks += 1L;

  } else {

    int m1 = args.m / 2;
    int m2 = args.m - m1;
    int n1 = args.n / 2;
    int n2 = args.n - n1;

    {
      struct dgemm_args args_1 = args;
      args_1.m = m1;
      args_1.n = n1;
      qmckl_dgemm_rec(args_1, tasks, ntasks);
    }

    {
      // TODO: assuming 'N', 'N' here
      struct dgemm_args args_2 = args;
      args_2.B = args.B + args.ldb*n1;
      args_2.C = args.C + args.ldc*n1;
      args_2.m = m1;
      args_2.n = n2;
      qmckl_dgemm_rec(args_2, tasks, ntasks);
    }

    {
      struct dgemm_args args_3 = args;
      args_3.A = args.A + m1;
      args_3.C = args.C + m1;
      args_3.m = m2;
      args_3.n = n1;
      qmckl_dgemm_rec(args_3, tasks, ntasks);
    }

    {
      struct dgemm_args args_4 = args;
      args_4.A = args.A + m1;
      args_4.B = args.B + args.ldb*n1;
      args_4.C = args.C + m1 + args.ldc*n1;
      args_4.m = m2;
      args_4.n = n2;
      qmckl_dgemm_rec(args_4, tasks, ntasks);
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
                 int64_t* tasks, int64_t* ntasks)
{

  struct dgemm_args* args = qmckl_dgemm_to_struct (transa, transb,
      m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

  starpu_memory_pin(A, lda*k*sizeof(double));
  starpu_memory_pin(B, ldb*n*sizeof(double));
  starpu_memory_pin(C, ldc*n*sizeof(double));
  qmckl_dgemm_rec(*args, tasks, ntasks);
  free(args);
}

void qmckl_tasks_run(struct dgemm_args** gemms, int ngemms)
{
  int NCPU, NGPU;
  sscanf( getenv( "STARPU_NCPU" ), "%d", &NCPU );
  sscanf( getenv( "STARPU_NCUDA" ), "%d", &NGPU );

  int rc = CHAMELEON_Init(NCPU, NGPU);

  for (int i=0 ; i<ngemms ; ++i)
    {
      CHAMELEON_dgemm(ChamNoTrans, ChamNoTrans,
		      gemms[i]->m,
		      gemms[i]->n,
		      gemms[i]->k,
		      gemms[i]->alpha,
		      gemms[i]->A,
		      gemms[i]->lda,
		      gemms[i]->B,
		      gemms[i]->ldb,
		      gemms[i]->beta,
		      gemms[i]->C,
		      gemms[i]->ldc);
    }
  CHAMELEON_Finalize();

}

