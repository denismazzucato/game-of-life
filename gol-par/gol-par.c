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

  MPI_Init(NULL, NULL); // argc, argv

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  // Get the name of the processor
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);

  // Print off a hello world message
  if (world_rank % 2 == 0 && world_rank < world_size) {
    int data = world_rank;
    int next_node = (world_rank + 1); // safe by the condition
    printf("%d", next_node);
    MPI_Send(&data, 1, MPI_INT, next_node, 0, MPI_COMM_WORLD);
  }

  if (world_rank % 2 == 1) {
    int buf;
    MPI_Recv(&buf, 1, MPI_INT, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    printf(
      "Hello world from processor %s, data %d, rank %d out of %d processors\n",
      processor_name,
      buf,
      world_rank,
      world_size
    );
    if (argc > 0) {
      for (int i = 0; i < argc; i++) {
        printf("%d-nth arg: '%s'\n", i, argv[i]);
      }
    }
  }

  MPI_Finalize();

  return 0;
}
