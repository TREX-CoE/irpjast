/* Generated from qmckl_dgemm.org */

#include <starpu.h>

#include <mkl_cblas.h>
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
  return args;

}

#define MIN_SIZE 512
static void qmckl_dgemm_rec(struct dgemm_args args, int64_t* tasks, int64_t* ntasks)
{


  if ( args.m * args.n <= MIN_SIZE*MIN_SIZE) {

//    printf("%5d %5d\n", args.m, args.n);
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

  qmckl_dgemm_rec(*args, tasks, ntasks);
  free(args);
}

void qmckl_tasks_run(struct dgemm_args** gemms, int ngemms)
{
  int rc = starpu_init(NULL);
  assert (rc == 0);


  starpu_data_handle_t matrix_handle[ngemms][3];
  for (int i=0 ; i<ngemms ; ++i)
    {
      starpu_matrix_data_register(&(matrix_handle[i][0]),
                                    STARPU_MAIN_RAM,
                                    (uintptr_t) gemms[i]->A,
                                    gemms[i]->lda,
                                    gemms[i]->m,
                                    gemms[i]->k,
                                    sizeof(double));

      starpu_matrix_data_register(&(matrix_handle[i][1]),
                                    STARPU_MAIN_RAM,
                                    (uintptr_t) gemms[i]->B,
                                    gemms[i]->ldb,
                                    gemms[i]->k,
                                    gemms[i]->n,
                                    sizeof(double));

      starpu_matrix_data_register(&(matrix_handle[i][2]),
                                    STARPU_MAIN_RAM,
                                    (uintptr_t) gemms[i]->C,
                                    gemms[i]->ldc,
                                    gemms[i]->m,
                                    gemms[i]->n,
                                    sizeof(double));

      struct starpu_task *task = starpu_task_create();

      task->cl = &dgemm_cl;
      task->cl_arg = gemms[i];
      task->cl_arg_size = sizeof(*gemms[0]);
      task->handles[0] = matrix_handle[i][0];
      task->handles[1] = matrix_handle[i][1];
      task->handles[2] = matrix_handle[i][2];
      rc = starpu_task_submit(task);
      assert (rc == 0);
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


