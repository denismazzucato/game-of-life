#!/bin/sh

EXEC='gol-par'
DIR="./gol/"

proc_number=1
if [ -n "$1" ]; then
  proc_number=$1
fi

cd $DIR
make gol-par
prun -v -1 -np $proc_number -script $PRUN_ETC/prun-openmpi ./$EXEC
make clean
