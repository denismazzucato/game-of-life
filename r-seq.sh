#!/bin/sh

EXEC='gol-seq'
DIR="./$EXEC/"

cd $DIR
make
prun -v -1 -np 1 -script $PRUN_ETC/prun-openmpi ./$EXEC.o
make clean
