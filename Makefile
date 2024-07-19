#Unix makefile for fortran-file	

# Put the name of the target program here
TARGET =  nucleosome-mer
# the list of source files
SRC =  mpivars.f90  precision.f90  mathconst.f90 physconst.f90 quaternion.f90 eigen.f90 Eq_triangle.f90 globals.f90 myutils.f90  molecule.f90  dielectfcn.f90  loop.f90 chains.f90 rands.f90 volume.f90 L2norm.f90 parameter.f90 GBpotential.f90  poissonEq.f90 field.f90 VdW.f90 surface.f90  fcnMg-expl.f90 confEntropy.f90  fenergy.f90 initcha.f90  myio.f90 rota.f90 cadenas.f90 cadenas-sequence.f90  fcn.Mg.f90 fcn.brush.f90  init.f90 VdWpotential.f90 chaingenerator.f90 kinsolsolver.f90  solver.f90  main.f90
# the object file
OBJ := $(patsubst %.f90, %.o, $(SRC))

# some definitions
SHELL = /bin/bash

# get git version
GIT_VERSION := $(shell git describe --abbrev=6 --dirty --always --tags)

ifeq ($(shell hostname),gadol)

FFLAGS= -cpp -DVERSION=\"$(GIT_VERSION)\" -fbounds-check -Warray-bounds -g -fbacktrace #-O3

LDFLAGS=-L/opt/local/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)

FF= mpif90


else ifeq ($(shell hostname),alfajor)

FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\" #-fbounds-check -Warray-bounds #-O3

LDFLAGS=  -O3 -lm /usr/lib/x86_64-linux-gnu/librt.so -L/opt/local/sundials-2.6.1-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/opt/local/sundials-2.6.1-openmpi/lib


LFFLAGS=$(LDFLAGS)


FF= mpif90




else ifeq ($(shell hostname),empanada)

FFLAGS=  -cpp -DVERSION=\"$(GIT_VERSION)\"  -fbounds-check -Warray-bounds #-O3

LDFLAGS=-L/opt/local/sundials-2.5.0-openmpi-2.0.2/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/5 -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/5/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/5/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)


FF= mpif90


else ifeq ($(shell hostname),ensalada)

FFLAGS=  -cpp -DVERSION=\"$(GIT_VERSION)\"  -fbounds-check -Warray-bounds #-O3                                                                 


LDFLAGS=  -lm /usr/lib/x86_64-linux-gnu/librt.so -L/opt/local/sundials-2.6.1-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/opt/local/sundials-2.6.1-openmpi/lib


LFFLAGS=$(LDFLAGS)

FF= mpif90

else ifeq ($(shell hostname),master)


FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS=-L/shared/software/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)
FF= /shared/software/openmpi-1.6.1/bin/mpif90

else ifeq ($(shell hostname),chiquita)


FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"

#LDFLAGS=-L/opt/local/sundials-2.5.0-openmpi-atlas/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1 -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../.. -lgfortran -lm -lgcc_s -lquadmath

LDFLAGS=  -lm /usr/lib/x86_64-linux-gnu/librt.so -L/opt/local/sundials-2.6.1-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/opt/local/sundials-2.6.1-openmpi/lib



LFFLAGS=$(LDFLAGS)
FF= mpif90

else ifeq ($(shell hostname),orange)

FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"

LDFLAGS=-L/opt/local/sundials-2.5.0-openmpi-atlas/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1 -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6.1/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)
FF= gfortran


else ifeq ($(shell hostname),pear)

FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"
LDFLAGS=-L/opt/local/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-linux-gnu/4.6 -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.6/../../.. -lgfortran -lm -lgcc_s -lquadmath

LFFLAGS=$(LDFLAGS)
FF= gfortran

else ifeq ($(shell hostname),master.bw01.bme.northwestern.edu)

FFLAGS= -O3 -cpp -DVERSION=\"$(GIT_VERSION)\"

LDFLAGS= -L/export/apps/sundials-2.5.0-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial -lm -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../../../lib64 -L/lib/../lib64 -L/usr/lib/../lib64 -L/usr/lib/gcc/x86_64-redhat-linux/4.4.6/../../.. -lgfortranbegin -lgfortran -lm -lgcc_s

LFFLAGS=$(LDFLAGS)

FF= gfortran


else ifeq ($(shell hostname),quser34)

        is_quest = yes

else ifeq ($(shell hostname),quser33)

	is_quest = yes

else ifeq ($(shell hostname),quser32)

        is_quest = yes

else ifeq ($(shell hostname),quser31)

	is_quest = yes

else 


FFLAGS=  -std=f2008 -cpp -DVERSION=\"$(GIT_VERSION)\" -fcheck=all -fbounds-check -Warray-bounds -g -fbacktrace  -Wpedantic #-Wall

#FFLAGS=  -std=f2008 -cpp -DVERSION=\"$(GIT_VERSION)\" -O3

#LDFLAGS=-lm -L/opt/local/kinsol-2.8.2-stat/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/opt/local/kinsol-2.8.2-stat/lib

LDFLAGS= -lm -L/opt/local/sundials-2.6.1-openmpi/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/opt/local/sundials-2.6.1-openmpi/lib -framework Accelerate


LFFLAGS=$(LDFLAGS)

FF= mpif90


endif


ifdef is_quest 


FFLAGS=  -std=f2008  -cpp -DVERSION=\"$(GIT_VERSION)\"  -O3 # -fcheck=all -fbounds-check -Warray-bounds -g -fbacktrace # -Wargument-mismatch -Wpedantic #-Wall


LDFLAGS= -lm /usr/lib64/librt.so -L/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran/lib -lsundials_fkinsol -lsundials_kinsol -lsundials_fnvecserial -lsundials_nvecserial     -Wl,-rpath,/projects/p31445/sundials/sundials-2.6.1-openmpi-gfortran/lib



LFFLAGS=$(LDFLAGS)

FF= mpif90

endif



all:	$(TARGET)

$(TARGET): $(SRC:.f90=.o)
	$(FF) -o $(TARGET) $(SRC:.f90=.o)  $(LDFLAGS) $(LFLAGS)


%.o: %.f90
	$(FF) -c $(FFLAGS) $< -o $@


install: all
ifdef is_cooley        
	cp $(TARGET) ~/bincooley
else 
	cp $(TARGET) ~/bin
endif



clean:	
	@rm -f $(SRC:.f90=.o) $(SRC:.f90=.d) $(TARGET) *~ *.mod

realclean: clean
	@rm -f .depend

depend dep:
	@$(FF)  $(CFLAGS) -MM $(SRC) > .depend 

ifeq (.depend, $(wildcard .depend))
include .depend
endif













