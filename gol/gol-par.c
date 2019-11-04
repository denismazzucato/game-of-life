#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

static int rank, numberOfNodes;

static int bwidth, // board width
					 bheight, // board height, in the seq algorithm bheight == nrows
					 nrows, // entire board height, not only the part of this node
					 nsteps; // time limit

// matrices (bwidth+2)x(bheight+2), to wrap-around the world structure
static int **old, **new;

static int print_world = 0;

// use fixed world or random world?
// static int random_world = 0;
static int random_world = 1;

// 22x42
char *start_world[] = {
    /* Gosper glider gun */
    /* example from https://bitstorm.org/gameoflife/ */
    "..........................................",
    "..........................................",
    "..........................................",
    "..........................................",
    "..........................................",
    "..........................................",
    "........................OO.........OO.....",
    ".......................O.O.........OO.....",
    ".OO.......OO...........OO.................",
    ".OO......O.O..............................",
    ".........OO......OO.......................",
    ".................O.O......................",
    ".................O........................",
    "....................................OO....",
    "....................................O.O...",
    "....................................O.....",
    "..........................................",
    "..........................................",
    ".........................OOO..............",
    ".........................O................",
    "..........................O...............",
    "..........................................",
};

int isMaster() {
	return rank == 0;
}

int isSlave() {
	return rank;
}


// I've already check that argv contains exactly 5 args
void initParameters(char *argv[]) {
	bwidth = atoi(argv[1]);
	nrows = atoi(argv[2]);
	nsteps = atoi(argv[3]);
	print_world = atoi(argv[4]);
}

void initMPI() {
	MPI_Init(NULL, NULL); // TODO: shall I use something instead of NULL, NULL?

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfNodes);
}

/* the array allocation can be done without communication, each processors
 * already has all the information to instantiate the old and new matrices
 */
// void allocateArrays() {
// 	int rowsPerNode = nrows / numberOfNodes; // this is in avarage

// 	int otherNodes = nrows % numberOfNodes; // nodes with N/p +1 rows
// 	int fullNodes = numberOfNodes - otherNodes; // nodes with N/p rows

// 	// each node in the first fullNodes holds N/p rows
// 	// each node in the last otherNodes holds (N/p)+1 rows
// 	int realRowsPerFullNodes = rowsPerNode + 2;
// 	int realRowsPerOtherNodes = realRowsPerFullNodes + 1; // rowsPerNode + 2 + 1

// 	int numberOfRealRows = (rank < fullNodes
// 		? realRowsPerFullNodes
// 		: realRowsPerOtherNodes);
// 	int numberOfRealColumns = bwidth + 2;

// 	old = malloc(numberOfRealRows * sizeof(int*));
// 	new = malloc(numberOfRealRows * sizeof(int*));

// 	for (int i = 0; i < numberOfRealRows; i++) {
// 		old[i] = malloc(numberOfRealColumns * sizeof(int*));
// 		new[i] = malloc(numberOfRealColumns * sizeof(int*));
// 	}

// 	// now we have to set the correct value for bheight
// 	// because then I want to iterate over bheight
// 	bheight = numberOfRealRows - 2;
// }

float randomNumber() {
	return rand() / ((float)RAND_MAX + 1);
}

int isInsideStartWorld(int i, int j) {
	return (i < sizeof(start_world) / sizeof(char *)) &&
				 (j < strlen(start_world[i]));
}

/* this function has to return a linear representation of the matrix
 *
 * TODO: send data as first as possible,
 * I don't need the entire matrix on the master node
 */
int* buildBoard() {
	int *tmpBoard = malloc((nrows*bwidth) * sizeof(int));
	for (int i = 0; i < bwidth; i++)
		for (int j = 0; j < nrows; j++)
			if (random_world)
				tmpBoard[i*nrows + j] = (randomNumber() < 0.5 ? 0 : 1);
			else if (isInsideStartWorld(i, j))
				tmpBoard[i*nrows + j] = (start_world[i][j] != '.');
			else
				tmpBoard[i*nrows + j] = 0;

	return tmpBoard;
}

void receiveRows() {
}

void scatterRows(int* data) {
	if (isMaster())
		for (int i = 0; i < bwidth; i++) {
			for (int j = 0; j < nrows; j++)
				printf("%d ", data[i*nrows + j]);
			printf("\n");
		}
}

/* the initialization of old/new can be done without communication with this
 * configuration, but this is due to the fact that the static strt_world is
 * very small
 * In the other hand, the big problem that we want to compute is random
 * generated, so each node can do this procedure itself
 * In a general way, the master will manage the very big input and than he
 * scatter the data to all the other nodes, so we need to store the whole big
 * matrix just in a single node
 */
void initializeBoard() {
	int *tempBoard;


	if (isMaster())
		tempBoard = buildBoard();

	scatterRows(tempBoard);

	if (isMaster())
		free(tempBoard);
}

void freeMatrices() {
	if(new && old) {

		for (int i = 0; i < bwidth + 2; i++) {
			free(old[i]);
			free(new[i]);
		}

		free(old);
		free(new);
	}
}

int main (int argc, char *argv[]) {

	initMPI();

	if (argc == 5)
		initParameters(argv);
	else if (isMaster()) {
		fprintf(stderr,
			"Usage: %s board_width board_height steps_count print_iter\n",
			argv[0]);

		MPI_Abort(MPI_COMM_WORLD, 1);
		MPI_Finalize();
		return 0;
	}
	// bwidth, nrows, nsteps, print_world are now initialized

	// allocateArrays(rank, numberOfNodes);
	// bheight is now initialized

	initializeBoard();

	MPI_Finalize();

	freeMatrices();
	return 0;
}