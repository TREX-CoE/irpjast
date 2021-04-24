/* Generated from qmckl_dgemm.org */

#include <starpu.h>

#include <cblas.h>
#include <stdint.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>



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


void qmckl_dgemm_cl(struct dgemm_args args, double* A, double* B, double* C);

void dgemm_codelet(void *buffers[], void* cl_arg)
{
  struct dgemm_args *args = cl_arg;
  double* A = (double*) STARPU_MATRIX_GET_PTR(buffers[0]);
  double* B = (double*) STARPU_MATRIX_GET_PTR(buffers[1]);
  double* C = (double*) STARPU_MATRIX_GET_PTR(buffers[2]);
  qmckl_dgemm_cl(*args, A, B, C);
  free(args);
}

struct starpu_codelet dgemm_cl =
  {
   .where = STARPU_CPU,
   .cpu_funcs = { dgemm_codelet },
   .cpu_funcs_name = { "dgemm_codelet" },
   .nbuffers = 3,
   .max_parallelism = 1,
   .modes = {STARPU_R, STARPU_R, STARPU_RW},

  };

#include<stdio.h>

void qmckl_dgemm_cl(struct dgemm_args args, double* A, double* B, double* C) {
    cblas_dgemm(CblasColMajor, args.transa, args.transb,
                args.m, args.n, args.k, args.alpha,
                A, args.lda, B, args.ldb,
                args.beta, C, args.ldc);
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
  starpu_init(NULL);
  
  starpu_data_handle_t matrix_handle[ngemms][3];
  for (int i=0 ; i<ngemms ; ++i)
    {
      starpu_matrix_data_register(&(matrix_handle[i][0]),
                                    STARPU_MAIN_RAM,
                                    (uintptr_t) gemms[i]->A,
                                    gemms[i]->lda,
                                    gemms[i]->k,
                                    gemms[i]->m,
                                    sizeof(double));

      starpu_matrix_data_register(&(matrix_handle[i][1]),
                                    STARPU_MAIN_RAM,
                                    (uintptr_t) gemms[i]->B,
                                    gemms[i]->ldb,
                                    gemms[i]->n,
                                    gemms[i]->k,
                                    sizeof(double));

      starpu_matrix_data_register(&(matrix_handle[i][2]),
                                    STARPU_MAIN_RAM,
                                    (uintptr_t) gemms[i]->C,
                                    gemms[i]->ldc,
                                    gemms[i]->n,
                                    gemms[i]->m,
                                    sizeof(double));

      struct starpu_task *task = starpu_task_create();

      task->cl = &dgemm_cl;
      task->cl_arg = gemms[i];
      task->cl_arg_size = sizeof(*gemms[0]);
      task->handles[0] = matrix_handle[i][0];
      task->handles[1] = matrix_handle[i][1];
      task->handles[2] = matrix_handle[i][2];
      starpu_task_submit(task);
    }
  starpu_task_wait_for_all();

  for (int i=0 ; i<ngemms ; ++i)
    {
      starpu_data_unregister(matrix_handle[i][0]);
      starpu_data_unregister(matrix_handle[i][1]);
      starpu_data_unregister(matrix_handle[i][2]);
    }
  starpu_shutdown();
}
