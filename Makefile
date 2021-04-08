FC     = ifort -O3 -ip -g  -march=core-avx2 -align array64byte -fma -ftz -fomit-frame-pointer -mkl=sequential #-O3 -march=native -mkl=sequential -g -I$(PWD) #-xHost -check all 
CPP    = g++ -O2 -march=native # -Wall
FCFLAGS= #-O2 -ffree-line-length-none -I .
NINJA  = ninja
AR = ar
ARCHIVE = ar crs
RANLIB = ranlib

SRC=
OBJ=
LIB= $(MAGMA_F90FLAGS) $(LDFLAGS) $(MAGMA_LIBS) -L${GLIB} -lstdc++  #magma_dgemm_async_gpu.o
MAGMA         = /p/software/juwelsbooster/stages/2020/software/magma/2.5.4-gcccoremkl-9.3.0-2020.2.254
MAGMADIR      = /p/software/juwelsbooster/stages/2020/software/magma/2.5.4-gcccoremkl-9.3.0-2020.2.254
FORTRAN       = /p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/lib64
GLIB          = /p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/lib64
CUDADIR      ?= /p/software/juwelsbooster/stages/2020/software/CUDA/11.0
OPENBLASDIR  ?= /p/software/juwelsbooster/stages/2020/software/GCC/
MAGMA_CFLAGS   := -DADD_ -I$(MAGMADIR)/include -I$(CUDADIR)/include
MAGMA_F90FLAGS := -I$(MAGMADIR)/include -Dmagma_devptr_t="integer(kind=8)"

MAGMA_LIBS   := #-L$(MAGMADIR)/lib -L$(CUDADIR)/lib64 -L$(OPENBLASDIR)/lib \
                -lmagma -lcublas -lcudart -lmkl

IRPF90 = irpf90 --codelet=elec_dist:2 -s tile_size:32
-include irpf90.make
export

irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile
	$(IRPF90)

#magma_dgemm_async_gpu.o: 
#	${CPP} $(CFLAGS) $(MAGMA_CFLAGS) -DCUBLAS_GFORTRAN -c magma_dgemm_async_gpu.cc -o magma_dgemm_async_gpu.o
#
#fortran.o: $(CUDADIR)/src/fortran.c
#	$(CC) $(CFLAGS) $(MAGMA_CFLAGS) -DCUBLAS_GFORTRAN -c -o $@ $<
#
#tiling_interface.o: tiling_interface.f90
#	$(FC) $(FFLAGS) -c -o $@ $<
