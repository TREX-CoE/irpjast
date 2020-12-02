IRPF90 = ~/irpf90/bin/irpf90 --codelet factor_een:200
#FC = gfortran
#FCFLAGS= -g -msse4.2  -fcheck=all -Waliasing -Wampersand -Wconversion -Wsurprising -Wintrinsics-std -Wno-tabs -Wintrinsic-shadow -Wline-truncation -Wreal-q-constant -Wuninitialized  -fbacktrace -ffpe-trap=zero,overflow,underflow -finit-real=nan
FC     = ifort -g
FCFLAGS= -O2 -xHost -I .
NINJA  = ninja
AR     = ar crs
RANLIB = ranlib

SRC=
OBJ=
LIB=

-include irpf90.make
export

irpf90.make: $(filter-out IRPF90_temp/%, $(wildcard */*.irp.f)) $(wildcard *.irp.f) $(wildcard *.inc.f) Makefile 
	$(IRPF90)
