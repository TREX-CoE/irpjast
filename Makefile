IRPF90 = irpf90  --codelet=factor_een_2:100000 #--codelet=factor_een:10000 
FC     = gfortran
FCFLAGS= -O2 -march=native -ffree-line-length-none -I .
NINJA  = ninja
ARCHIVE= ar crs
RANLIB = ranlib

SRC=
OBJ=
LIB=

-include irpf90.make
export

irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile 
	$(IRPF90)
