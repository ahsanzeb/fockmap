##################################
# makefile: fockmap
##################################

FC 	= gfortran
FFLAGS= -O3 -g

SRC = factorial.f  binomial.f mapping.f mapn2.f  triangles.f \
 sizes.f basisall.f basis.f ratiosall.f ratios.f  

OBJ = $(SRC:.f=.o)

.SUFFIXES: .o .f

fockmap.a: $(OBJ)
	ar cr fockmap.a $(OBJ)
	cp fockmap.a ../test/.

.f.o: 
	$(FC) $(FFLAGS) -c $< -o $@

clean:
	rm -rf *.o fockmap.a
  
