IRPF90 = irpf90  #-a -d
FC     = gfortran -g
FCFLAGS= -O2 -ffree-line-length-none -I .
FC     = ifort -C -traceback -g
#FCFLAGS= -O2 -ffree-line-length-none -I .
#NINJA  = ninja
AR     = ar
ARCHIVE= ar crs
RANLIB = ranlib

SRC=
OBJ=
LIB=

-include irpf90.make
export

irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile 
	$(IRPF90)
