CC=mpif90 -fbounds-check
PROG = run
SRC = Data.f90 Num.f90 Func.f90 System.f90 GC.f90
OBJ = $(SRC:.f90=.o)

vpath %.f90 src_MPI

all : $(OBJ) Main.o
	$(CC) $^ -o $(PROG)

%.o: %.f90
	$(CC) -c $<

.PHONY: clean
clean :
	rm -f *.o *.mod *~ $(PROG)
