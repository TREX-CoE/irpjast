FC     = ifort -O3 -ip -g  -march=native -align array64byte -fma -ftz -fomit-frame-pointer -mkl=parallel #-O3 -march=native -mkl=sequential -g -I$(PWD) #-xHost -check all 
CC     = gcc
CPP    = g++ -O2 -march=native # -Wall
FCFLAGS= -I .#-O2 -ffree-line-length-none -I .
NINJA  = ninja
AR = ar
ARCHIVE = ar crs
RANLIB = ranlib

SRC=  tiling_interface.f90
OBJ=  tiling_interface.o
LIB= $(MAGMA_F90FLAGS) $(LDFLAGS) $(MAGMA_LIBS) $(STARPU_CFLAGS) $(STARPU_LIBS) -L${GLIB} -lstdc++  magma_dgemm_async_gpu.o dgemm.o gemm/common/blas.o
MAGMA         = /p/software/juwelsbooster/stages/2020/software/magma/2.5.4-gcccoremkl-9.3.0-2020.2.254
MAGMADIR      = /p/software/juwelsbooster/stages/2020/software/magma/2.5.4-gcccoremkl-9.3.0-2020.2.254
FORTRAN       = /p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/lib64
GLIB          = /p/software/juwelsbooster/stages/2020/software/GCCcore/9.3.0/lib64
CUDADIR      ?= /p/software/juwelsbooster/stages/2020/software/CUDA/11.0
OPENBLASDIR  ?= /p/software/juwelsbooster/stages/2020/software/GCC/
MAGMA_CFLAGS   := -DADD_ -I$(MAGMADIR)/include -I$(CUDADIR)/include
MAGMA_F90FLAGS := -I$(MAGMADIR)/include -Dmagma_devptr_t="integer(kind=8)"

MAGMA_LIBS   := -L$(MAGMADIR)/lib -L$(CUDADIR)/lib64 -L$(OPENBLASDIR)/lib \
                -lmagma -lcublas -lcudart -lmkl
## STAR PU ###
STARPU_VERSION=1.3
CPPFLAGS += $(shell pkg-config --cflags starpu-$(STARPU_VERSION))
LDLIBS += $(shell pkg-config --libs starpu-$(STARPU_VERSION))
STARPU_LIBS += $(shell pkg-config --libs starpu-$(STARPU_VERSION))
STARPU_CFLAGS += $(shell pkg-config --cflags starpu-$(STARPU_VERSION))

CFLAGS += -O3 -Wall -Wextra

# to avoid having to use LD_LIBRARY_PATH
LDLIBS += -Wl,-rpath -Wl,$(shell pkg-config --variable=libdir starpu-$(STARPU_VERSION))

# Automatically enable CUDA / OpenCL
STARPU_CONFIG=$(shell pkg-config --variable=includedir starpu-$(STARPU_VERSION))/starpu/$(STARPU_VERSION)/starpu_config.h
STARPU_CFLAGS=-I$(shell pkg-config --variable=includedir starpu-$(STARPU_VERSION))/starpu/$(STARPU_VERSION)
ifneq ($(shell grep "USE_CUDA 1" $(STARPU_CONFIG)),)
USE_CUDA=1
endif
ifneq ($(shell grep "USE_OPENCL 1" $(STARPU_CONFIG)),)
USE_OPENCL=1
endif
ifneq ($(shell grep "RELEASE_VERSION 99" $(STARPU_CONFIG)),)
USE_ENERGY=1
endif


IRPF90 = irpf90 --codelet=elec_dist:2 -s tile_size:32
-include irpf90.make
export

irpf90.make: fortran.o tiling_interface.o magma_dgemm_async_gpu.o dgemm.o $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile
	$(IRPF90)

magma_dgemm_async_gpu.o: 
	${CPP} $(CFLAGS) $(MAGMA_CFLAGS) -DCUBLAS_GFORTRAN -c magma_dgemm_async_gpu.cc -o magma_dgemm_async_gpu.o

fortran.o: $(CUDADIR)/src/fortran.c
	$(CC) $(CFLAGS) $(MAGMA_CFLAGS) -DCUBLAS_GFORTRAN -c -o $@ $<

tiling_interface.o: tiling_interface.f90
	$(FC) $(FFLAGS) -c -o $@ $<


#%.o: %.cu
#	nvcc $(CPPFLAGS) $< -c -o $@

CFLAGS+=-DSTARPU_OPENBLAS=0

dgemm.o:
	$(CC) $(CFLAGS) $(STARPU_CFLAGS) -c gemm/dgemm.c 

gemm/dgemm: gemm/dgemm.o gemm/common/blas.o

gemm/dgemm: LDLIBS+=-lmkl_intel_lp64
ifeq ($(USE_CUDA),1)
gemm/dgemm: LDLIBS+=-L$(CUDA_PATH)/lib64 -lcublas -lcudart
endif

