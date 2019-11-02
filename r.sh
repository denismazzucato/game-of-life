#!/bin/sh

EXEC='gol-par'
DIR="./$EXEC/"

proc_number=1
if [ -n "$1" ]; then
  proc_number=$1
fi

cd $DIR
make
prun -v -1 -np $proc_number -script $PRUN_ETC/prun-openmpi ./$EXEC.o
make clean
