GF=gfortran
SRC=module.f90 main.f90
OBJ=${SRC:.f90=.o}

%.o: %.f90
	$(GF) -o $@ -c $<

program.out: $(OBJ)
	$(GF) -o $@ $(OBJ)

run: program.out
	./program.out

clean:
	@rm -f *.mod *.o program.out