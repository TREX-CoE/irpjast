IRPF90 = irpf90  #-a -d
FC     = gfortran
FCFLAGS= -O2 -ffree-line-length-none -I .
NINJA  = ninja
ARCHIVE= ar crs
RANLIB = ranlib

SRC=jastrow_transirp.f
OBJ=jastrow_transirp.o
LIB=

-include irpf90.make
export

irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile 
	$(IRPF90)
