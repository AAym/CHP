CC=mpif90 -fbounds-check
PROG_MPI = run_mpi
PROG = run
SRC_MPI = src_MPI/Data.f90 src_MPI/Num.f90 src_MPI/Func.f90 src_MPI/System.f90 src_MPI/GC.f90 src/Main.f90
SRC = src/Data.f90 src/Num.f90 src/Func.f90 src/System.f90 src/GC.f90 src/Main.f90
OBJ = $(SRC:.f90=.o)
OBJ_MPI = $(SRC_MPI:.f90=.o)

vpath %.f90 src_MPI src

mpi : $(OBJ_MPI) 
	$(CC) $^ -o $(PROG_MPI)

seq: CC=gfortran
seq: $(OBJ) 
	$(CC) $^ -o $(PROG)


%.o: %.f90
	$(CC) -c $< -o $@

.PHONY: clean
clean :
	rm -rf *.o *.mod *~ $(PROG) $(PROG_MPI)
