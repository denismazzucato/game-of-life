/***********************

Conway's Game of Life

************************/

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

int
main(int argc, char *argv[])
{
  if (argc > 0)
    for (int i = 0; i < argc; i++)
      printf("%d-nth arg: '%s'\n", i, argv[i]);

  return 0;
}
