STARPU=~/StarPU/starpu-1.3.7/starpu-1.3.pc
IRPF90 = irpf90/bin/irpf90 --codelet=factor_een:2 --align=64 # -s nelec_8:504 -s nnuc:100 -s ncord:5 #-a -d
#FC     = ifort -mkl=sequential  -xCORE-AVX512 -g -qopt-zmm-usage=high 
FC     = ifort -mkl=sequential -march=core-avx2 -align array64byte -fma -ftz -fomit-frame-pointer
CC     = gcc $(shell pkg-config --cflags $(STARPU) chameleon) 
FCFLAGS= -O3 -I .
NINJA  = ninja
ARCHIVE = ar crs
RANLIB = ranlib

SRC= qmckl_blas_f.f90 qmckl_dgemm.c
OBJ= IRPF90_temp/qmckl_blas_f.o   IRPF90_temp/qmckl_dgemm.o
LIB= -mkl=sequential $(shell pkg-config --libs $(STARPU) chameleon magma)

-include irpf90.make
export

irpf90.make: IRPF90_temp/qmckl_blas_f.o
irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile
	$(IRPF90)

IRPF90_temp/%.f90: %.f90
IRPF90_temp/%.c: %.c

IRPF90_temp/%.o: %.f90
	$(FC) $(FCFLAGS) -g -c $< -o $@

IRPF90_temp/%.o: %.c
	$(CC) -g -c $< -o $@
