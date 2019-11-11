clean-out:
	rm -r ./out
	mkdir out

all:
	cd gol
	make gol-seq
	make gol-par
	make gol-par-opt
