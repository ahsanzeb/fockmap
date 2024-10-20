##################################
# makefile: test
##################################

FC 	= gfortran
FFLAGS= -O3 -g

SRC = main.f

OBJ = $(SRC:.f=.o)

.SUFFIXES: .o .f

test: main.o fockmap.a
	$(FC) -o test $(FFLAGS) main.o fockmap.a

main.o: main.f
	 $(FC) -c $(FFLAGS) main.f	

clean:
	rm -rf *.o *.mod test
  
