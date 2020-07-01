.PHONY: clean default debug release extrae run run-trace

FILE_NAMES 	?= matrix_example.cpp
PROG 		?= main
MATRIX_SIZE ?= 300
RUN_SETTINGS=OMP_PLACES=cores OMP_PROC_BIND=close I_MPI_PIN=1 I_MPI_PIN_DOMAIN=auto

APP_COMPILE_FLAGS=-std=c++11 -qopenmp
APP_LINKER_FLAGS=-lm -lstdc++

# default is release
FLAGS_OPTIMIZATION=-g -O3

default: release

release: FLAGS_OPTIMIZATION=-g -O3
release: generic

extrae: FLAGS_OPTIMIZATION=-g -O3
extrae: APP_COMPILE_FLAGS+= -I${EXTRAE_HOME}/include -DTRACE_EXTRAE
extrae: APP_LINKER_FLAGS+= -L${EXTRAE_HOME}/lib -lompitrace -lpapi
extrae: generic

debug: FLAGS_OPTIMIZATION=-g -O0
debug: generic

generic: clean
	${MPICXX}  $(FLAGS_OPTIMIZATION) -o $(PROG) $(APP_COMPILE_FLAGS) $(FILE_NAMES) $(APP_LINKER_FLAGS)

run:
	$(RUN_SETTINGS) mpiexec.hydra -np 4 -genvall $(PROG) $(MATRIX_SIZE) 700 500 300 200

run-trace:
	$(RUN_SETTINGS) mpiexec.hydra -np 4 -genvall ./trace.sh $(PROG) $(MATRIX_SIZE) 700 500 300 200

clean:
	rm -f $(PROG) *.o 
