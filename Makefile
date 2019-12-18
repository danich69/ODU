rk: prog
	./prog rk
ae: prog
	time ./prog ae
ai: prog
	time ./prog ai
prog: main.o odu.o non_linear_system.o diffur.o linear_system.o
	gfortran $^ -o $@ -g
main.o: main.f90 odu.mod
	gfortran main.f90 -c -g
odu.mod odu.o: odu.f90 non_linear_system.mod diffur.mod
	gfortran $^ -c -g
non_linear_system.mod non_linear_system.o: non_linear_system.f90 linear_system.mod
	gfortran $^ -c -g
linear_system.mod linear_system.o: linear_system.f90
	gfortran $^ -c -g
diffur.mod diffur.o: diffur.f90
	gfortran $^ -c -g
clean: 
	rm -f *.o *mod
