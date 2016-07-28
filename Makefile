LIBS = -L/TEMP/LIB -L/opt/intel/mkl/10.1.2.024/lib/em64t -larpack64new  -llapack64new -lblas64new 
F90 =ifort

PROG = ./hydrogen_fedvr

F90FLAGS = -O3 -static -openmp -parallel -fpp -u -132

LDFLAGS=$(F90FLAGS)



SRCS =	hydrogen_fedvr.f90                \
	accuracy_real.f90                 \
	globalmod_constants.f90           \
	globalmod_discrete.f90            \
	readin.f90 hamiltonian_matrix.f90 \
	lagrange_dvr_derivative.f90       \
	lgngr.f90 gaussq.f90		  \
	expokit.f90 mataid.f90   	  \
	localincludes.f90 prop.f90        \
	laser.f90



OBJS =	$(SRCS:.f90=.o)

all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

clean:
	rm -f $(PROG) $(OBJS) *.mod *~

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<


hydrogen_fedvr.o          : hydrogen_fedvr.f90 accuracy_real.o globalmod_constants.o globalmod_discrete.o localincludes.o laser.o \
                            readin.o hamiltonian_matrix.o
accuracy_real.o           : accuracy_real.f90
globalmod_constants.o     : globalmod_constants.f90 accuracy_real.o 
globalmod_discrete.o      : globalmod_discrete.f90 accuracy_real.o
readin.o                  : readin.f90 accuracy_real.o globalmod_constants.o globalmod_discrete.o
hamiltonian_matrix.o      : hamiltonian_matrix.f90 accuracy_real.o globalmod_constants.o globalmod_discrete.o
lagrange_dvr_derivative.o : lagrange_dvr_derivative.f90 accuracy_real.o lgngr.o gaussq.o
lgngr.o                   : lgngr.f90
gaussq.o                  : gaussq.f90
localincludes.o		  : localincludes.f90
laser.o		      	  : laser.f90
