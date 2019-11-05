gol-seq: gol-seq.c
	gcc -Wall -O3 -o gol-seq gol-seq.c -lm

gol-par: gol-par.c
	mpicc -Wall -O3 -o gol-par gol-par.c -lm

gol-par-opt: gol-par-opt.c
	mpicc -Wall -O3 -o gol-par-opt gol-par-opt.c -lm

goal-par-das: gol-par
	prun -v -1 -np 16 -script $PRUN_ETC/prun-openmpi ./gol-par 3000 3000 1000 0

goal-par-loc: gol-par
	mpirun -np 6 gol-par 3000 3000 1000 0

goal-par-opt-das: gol-par-opt
	prun -v -1 -np 16 -script $PRUN_ETC/prun-openmpi ./gol-par-opt 3000 3000 1000 0

goal-par-opt-loc: gol-par-opt
	mpirun -np 6 gol-par-opt 3000 3000 1000 0


clean:
	rm -f *.o *~ *core gol-par-opt gol-par gol-seq
