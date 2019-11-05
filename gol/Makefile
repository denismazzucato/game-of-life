gol-seq: gol-seq.c
	gcc -Wall -O3 -o gol-seq gol-seq.c -lm

gol-par: gol-par.c
	mpicc -Wall -O3 -o gol-par gol-par.c -lm

goal-par: gol-par
	prun -v -1 -np 16 -script $PRUN_ETC/prun-openmpi ./gol-par

clean:
	rm -f *.o *~ *core
