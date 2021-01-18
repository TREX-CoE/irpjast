IRPF90 = irpf90 --codelet=factor_een:100000 #-a -d
FC     = ifort -xHost -g -mkl=sequential
FCFLAGS= -O2 -ffree-line-length-none -I .
NINJA  = ninja
AR = ar 
ARCHIVE = ar crs
RANLIB = ranlib

SRC=
OBJ=
LIB=

-include irpf90.make
export

irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile 
	$(IRPF90)
