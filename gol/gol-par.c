#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

static int PRINT_TAG = 28301321;

static int rank, numberOfNodes;

static int bwidth, // board width
					 bheight, // board height, in the seq algorithm bheight == nrows
					 nrows, // entire board height, not only the part of this node
					 nsteps; // time limit

// matrices (bwidth+2)x(bheight+2), to wrap-around the world structure
static int **old = 0, **new = 0;

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

void debug(void) {
	MPI_Barrier(MPI_COMM_WORLD);
	fprintf(stderr,"\n--------------------- debug [%d] -------------------------", rank);
	MPI_Barrier(MPI_COMM_WORLD);
}

void abort(void) {
	MPI_Abort(MPI_COMM_WORLD, 1);
	MPI_Finalize();
	exit(1);
}

void printCellUnsafe(void) {
	for (int i = 1; i <= bheight && old; i++) {
		fprintf(stdout,"\n(Node: %3d; Line: %3d): ", rank, i);
		fflush(stdout);
		for (int j = 1; j <= bwidth && old[i]; j++) {
			fprintf(stdout,"%s", (old[i][j]) ? "0" : "_");
			fflush(stdout);
		}
	}
}

void printCells(void) {
	int m = 0; // fake message to provide a correct int to send and receive
	if (rank == 0) {
		printCellUnsafe();
		if (rank + 1 != numberOfNodes)
    	MPI_Send(&m, 1, MPI_INT, rank+1, PRINT_TAG, MPI_COMM_WORLD);
	} else {
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, PRINT_TAG, MPI_COMM_WORLD, &status);
		MPI_Recv(&m, 1, MPI_INT, MPI_ANY_SOURCE, PRINT_TAG, MPI_COMM_WORLD, &status);
		printCellUnsafe();
		if (rank + 1 != numberOfNodes)
			MPI_Send(&m, 1, MPI_INT, rank+1, PRINT_TAG, MPI_COMM_WORLD);
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

int isMaster(void) {
	return rank == 0;
}

int isSlave(void) {
	return rank;
}


// I've already check that argv contains exactly 5 args
void initParameters(char *argv[]) {
	bwidth = atoi(argv[1]);
	nrows = atoi(argv[2]);
	nsteps = atoi(argv[3]);
	print_world = atoi(argv[4]);
}

void initMPI(void) {
	MPI_Init(NULL, NULL); // TODO: shall I use something instead of NULL, NULL?

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numberOfNodes);
}

/* the array allocation can be done without communication, each processors
 * already has all the information to instantiate the old and new matrices
 *
 * the rows that exceed from splitting in perfect halves the data will be
 * distribuited in the last nodes
 */
void allocateArrays(void) {
	int rowsPerNode = nrows / numberOfNodes; // this is in avarage

	int otherNodes = nrows % numberOfNodes; // nodes with N/p +1 rows
	int fullNodes = numberOfNodes - otherNodes; // nodes with N/p rows

	// each node in the first fullNodes holds N/p rows
	// each node in the last otherNodes holds (N/p)+1 rows
	int realRowsPerFullNodes = rowsPerNode + 2;
	int realRowsPerOtherNodes = realRowsPerFullNodes + 1; // rowsPerNode + 2 + 1

	int numberOfRealRows = (rank < fullNodes
		? realRowsPerFullNodes
		: realRowsPerOtherNodes);
	int numberOfRealColumns = bwidth + 2;

	old = malloc(numberOfRealRows * sizeof(int*));
	new = malloc(numberOfRealRows * sizeof(int*));

	for (int i = 0; i < numberOfRealRows; i++) {
		old[i] = malloc(numberOfRealColumns * sizeof(int*));
		new[i] = malloc(numberOfRealColumns * sizeof(int*));
	}

	// now we have to set the correct value for bheight
	// because then I want to iterate over bheight
	bheight = numberOfRealRows - 2;

	// fprintf(stderr,
	// 	"Node %d: {\n\t rowsPerFullNode: %d\n\t otherNodes: %d\n\t fullNodes: %d\n\t numberOfRealRows: %d\n\t numberOfRealColumns: %d\n\t bheight: %d\n\t bwidth: %d}\n",
	// 	rank,
	// 	rowsPerNode,
	// 	otherNodes,
	// 	fullNodes,
	// 	numberOfRealRows,
	// 	numberOfRealColumns,
	// 	bheight,
	// 	bwidth
	// 	);
}

float randomNumber(void) {
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
int* buildBoard(void) {
	int *tmpBoard = malloc((nrows*bwidth) * sizeof(int));

	fprintf(stderr, "\n[Initial Matrix:\n");
	for (int i = 0; i < bwidth; i++) {
		fprintf(stderr, " > ");
		for (int j = 0; j < nrows; j++) {
			if (random_world)
				tmpBoard[i*nrows + j] = (randomNumber() < 0.5 ? 0 : 1);
			else if (isInsideStartWorld(i, j))
				tmpBoard[i*nrows + j] = (start_world[i][j] != '.');
			else
				tmpBoard[i*nrows + j] = 0;
			fprintf(stderr, "%s", (tmpBoard[i*nrows + j]) ? "0" : "_");
		}
		fprintf(stderr, "\n");
	}
	fprintf(stderr, "]\n");


	return tmpBoard;
}

void buildMatrixFromBuffer(int* recbuf) {
	int k = 0;
	for (int i = 1; i <= bheight; i++)
		for (int j = 1; j <= bwidth; j++)
			old[i][j] = recbuf[k++]; // k == (i-1)*bwidth + (j-1)
}

void scatterRows(int* data) {
	int *recbuf = malloc(bheight*bwidth * sizeof(int));
  int *sendcounts = malloc(numberOfNodes * sizeof(int));
  int *displs = malloc(numberOfNodes * sizeof(int));

	int sum = 0;
	int iter, otherNodes = nrows % numberOfNodes; // nodes with N/p +1 rows

  for (iter = 0; iter < numberOfNodes-otherNodes; iter++)
    sendcounts[iter] = bwidth*(nrows/numberOfNodes);
  for (;iter < numberOfNodes; iter++)
    sendcounts[iter] = bwidth*(nrows/numberOfNodes + 1);
  for (int i = 0; i < numberOfNodes; i++) {
    displs[i] = sum;
    sum += sendcounts[i];
  }

	// divide the data among processes as described by sendcounts and displs
  MPI_Scatterv(data, sendcounts, displs, MPI_INT, recbuf, bheight*bwidth, MPI_INT, 0, MPI_COMM_WORLD);

	buildMatrixFromBuffer(recbuf);

	free(recbuf);
	free(sendcounts);
	free(displs);
	if (data && isMaster()) free(data);
}

/* this function built all the line,
 * when (~)N/P line are created are immediately sent to the corresponding node
 *
 * this is because I don't want to store all the data in a single node even
 * only for the initialization because this table may be too large to be stored
 */
void initializeBoard(void) {
	int *tempBoard = 0;

	if (isMaster()) tempBoard = buildBoard();

	scatterRows(tempBoard);

	printCells();
}

void freeMatrices(void) {
	if(new && old) {

		for (int i = 0; i < bheight + 2; i++) {
			free(old[i]);
			free(new[i]);
		}

		free(old);
		free(new);
	}

	fprintf(stdout, "\n");
	fflush(stdout);
}

void argsNotProvided(void) {
	if (isMaster()) {
		fprintf(
			stderr,
			"Usage: board_width board_height steps_count print_iter\n");
		abort();
	}
}

void checkArgs(void) {
	if (nrows <= numberOfNodes && isMaster()) {
		fprintf(
			stderr,
			"Too few rows provided, provide at least %d rows\n",
			(numberOfNodes + 1));
		abort();
	}
}

int main (int argc, char *argv[]) {

	initMPI();

	if (argc == 5) {
		initParameters(argv);
		checkArgs();
	} else
		argsNotProvided();
	// bwidth, nrows, nsteps, print_world are now initialized and correct

	allocateArrays();
	// bheight is now initialized

	initializeBoard();

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	freeMatrices();
	return 0;
}