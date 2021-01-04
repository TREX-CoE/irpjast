IRPF90 = irpf90  #-a -d
FC     = ifort
FCFLAGS= -O2 -I .
NINJA  = ninja
ARCHIVE= ar crs
RANLIB = ranlib

SRC=inputs.f90
OBJ=
LIB=

-include irpf90.make
export

irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile 
	$(IRPF90)
