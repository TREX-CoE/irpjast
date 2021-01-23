IRPF90 = irpf90 --codelet=jastrow_full:1000 #-s nelec:10 -s nnuc:2 -s ncord:5 #-a -d
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
