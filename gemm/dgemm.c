/* StarPU --- Runtime system for heterogeneous multicore architectures.
 *
 * Copyright (C) 2009, 2010  Université de Bordeaux 1
 * Copyright (C) 2010  Centre National de la Recherche Scientifique
 *
 * StarPU is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * StarPU is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * See the GNU Lesser General Public License in COPYING.LGPL for more details.
 */

#include "double.h"
#include "xgemm.c" 

//int main(int argc, char **argv)
//{
//	double start, end;
//	int ret;
//  int nsteps, ns;
//
//  // Nsteps - Number of DGEMM tasks
//  nsteps = 40;
//
//	parse_args(argc, argv);
//
//#ifdef STARPU_QUICK_CHECK
//	niter /= 10;
//#endif
//
//	ret = starpu_init(NULL);
//	if (ret == -ENODEV)
//		return 77;
//	STARPU_CHECK_RETURN_VALUE(ret, "starpu_init");
//
//	starpu_cublas_init();
//
//	init_problem_data();
//	partition_mult_data();
//
//	start = starpu_timing_now();
//
//	unsigned x, y, iter;
//	for (iter = 0; iter < niter; iter++)
//	{
//		for (ns = 0; ns < nsteps; ++ns)
//		for (x = 0; x < nslicesx; ++x)
//		for (y = 0; y < nslicesy; ++y)
//		{
//			struct starpu_task *task = starpu_task_create();
//
//			task->cl = &cl;
//
//			task->handles[0] = starpu_data_get_sub_data(A_handle, 1, y);
//			task->handles[1] = starpu_data_get_sub_data(B_handle, 1, x);
//			task->handles[2] = starpu_data_get_sub_data(C_handle, 2, x, y);
//
//			ret = starpu_task_submit(task);
//			if (ret == -ENODEV)
//			{
//			     ret = 77;
//			     goto enodev;
//			}
//			STARPU_CHECK_RETURN_VALUE(ret, "starpu_task_submit");
//		}
//
//		starpu_task_wait_for_all();
//	}
//
//
//	end = starpu_timing_now();
//
//	double timing = end - start;
//
//	FPRINTF(stderr, "Time: %2.2f ms\n", timing/1000.0);
//
//	double flops = 2.0*((unsigned long)niter)*((unsigned long)xdim)
//				*((unsigned long)ydim)*((unsigned long)zdim*(unsigned long)nsteps);
//	FPRINTF(stderr, "GFlop/s: %.2f\n", flops/timing/1000.0);
//
//enodev:
//	starpu_data_unpartition(C_handle, 0);
//	starpu_data_unpartition(B_handle, 0);
//	starpu_data_unpartition(A_handle, 0);
//
//	starpu_data_unregister(A_handle);
//	starpu_data_unregister(B_handle);
//	starpu_data_unregister(C_handle);
//
//#ifndef STARPU_SIMGRID
//	if (check)
//		check_output();
//#endif
//
//	starpu_free(A);
//	starpu_free(B);
//	starpu_free(C);
//
//	starpu_cublas_shutdown();
//	starpu_shutdown();
//
//	return ret;
//}
