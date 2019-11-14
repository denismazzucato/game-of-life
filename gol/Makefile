gol-seq: gol-seq.c
	gcc -Wall -O3 -o gol-seq gol-seq.c -lm

gol-par: gol-par.c
	mpicc -Wall -O3 -o gol-par gol-par.c -lm

gol-par-opt: gol-par-opt.c
	mpicc -Wall -O3 -o gol-par-opt gol-par-opt.c -lm

gol-par-nbl: gol-par-nbl.c
	mpicc -Wall -O3 -o gol-par-nbl gol-par-nbl.c -lm

clean:
	rm -f *.o *~ *core gol-par-opt gol-par gol-seq
