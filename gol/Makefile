gol-seq: gol-seq.c
	gcc -Wall -O3 -o gol-seq gol-seq.c -lm

gol-par: gol-par.c
	mpicc -Wall -O3 -o gol-par gol-par.c -lm

clean:
	rm -f *.o *~ *core

run: gol-par
	mpirun gol-par a b c
