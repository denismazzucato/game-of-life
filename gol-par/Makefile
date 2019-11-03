gol-par: gol-par.c
	mpicc -Wall -O3 -o gol-par.o gol-par.c -lm

run: gol-par
	mpirun gol-par.o a b c

clean:
	rm -f *.o gol-par.o *~ *core