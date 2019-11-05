#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

static int
	PRINT_TAG = 0,
	FIRST_ROW_TAG = 1,
	LAST_ROW_TAG = 2,
	ROOT = 0;

static int rank, numberOfNodes;

static int bwidth, // board width
					 bheight, // board height, in the seq algorithm bheight == nrows
					 nrows, // entire board height, not only the part of this node
					 nsteps; // time limit

// matrices (bwidth+2)x(bheight+2), to wrap-around the world structure
static int **old = 0, **new = 0;

static int print_world = 0;

// use fixed world or random world?
static int random_world = 0; // fixed
// static int random_world = 1; // random

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
		// fprintf(stdout,"(Node: %3d; Line: %3d): ", rank, i);
		fflush(stdout);
		for (int j = 1; j <= bwidth && old[i]; j++) {
			fprintf(stdout,"%s", (old[i][j]) ? "0" : "_");
			fflush(stdout);
		}
		fprintf(stdout, "\n");
	}
}

// TODO: this may be optimized using ghater instruction
void printCells(int timestep) {
	int m = 0; // fake message to provide a correct int to send and receive

	MPI_Barrier(MPI_COMM_WORLD);
	if (rank == 0) {
		if (timestep > 0) printf("\nafter time step %d:\n\n", timestep);
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
	return rank == ROOT;
}

int isLast(void) {
	return rank == (numberOfNodes - 1);
}


// I've already check that argv contains exactly 5 args
void initParameters(char *argv[]) {
	nrows = atoi(argv[1]);
	bwidth = atoi(argv[2]);
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
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < bwidth; j++) {
			if (random_world)
				tmpBoard[i*bwidth + j] = (randomNumber() < 0.5 ? 0 : 1);
			else if (isInsideStartWorld(i, j))
				tmpBoard[i*bwidth + j] = (start_world[i][j] != '.');
			else
				tmpBoard[i*bwidth + j] = 0;
		}
	}

	return tmpBoard;
}

void buildMatrixFromBuffer(int* recbuf) {
	int k = 0;
	for (int i = 1; i <= bheight; i++)
		for (int j = 1; j <= bwidth; j++)
			old[i][j] = recbuf[k++]; // k == (i-1)*bwidth + (j-1)

}

void scatterRows(int* data) {
	int	*recbuf = malloc(bheight*bwidth * sizeof(int)),
		  *sendcounts = malloc(numberOfNodes * sizeof(int)),
		  *displs = malloc(numberOfNodes * sizeof(int));
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
  MPI_Scatterv(data, sendcounts, displs, MPI_INT, recbuf, bheight*bwidth, MPI_INT, ROOT, MPI_COMM_WORLD);

	buildMatrixFromBuffer(recbuf);

	free(recbuf);
	free(sendcounts);
	free(displs);
}

/* this function built all the line,
 * when (~)N/P line are created are immediately sent to the corresponding node
 *
 * this is because I don't want to store all the data in a single node even
 * only for the initialization because this table may be too large to be stored
 */
void initializeBoard(void) {
	int *tempBoard = 0;

	if (isMaster()) {
		tempBoard = buildBoard();

		// fprintf(stderr, "\n[Initial Matrix:\n");
		// for (int i = 0; i < nrows; i++) {
		// 	fprintf(stderr, " %2d > ", i);
		// 	for (int j = 0; j < bwidth; j++)
		// 		fprintf(stderr, "%s", (tempBoard[i*nrows + j]) ? "0" : "_");
		// 	fprintf(stderr, "\n");
		// }
		// fprintf(stderr, "]\n");
	}

	scatterRows(tempBoard);

	if (tempBoard && isMaster()) free(tempBoard);
}

void freeMatrix(int **data) {
	if(data) {
		for (int i = 0; i < bheight + 2; i++)
			free(data[i]);

		free(data);
	}
}

void argsNotProvided(void) {
	if (isMaster()) {
		fprintf(
			stderr,
			"Usage: board_height board_width steps_count print_iter\n");
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

// TODO: maybe with caotic search we can overlap communication e computation
void exchangeColumn(
		// MPI_Request *reqSendFirstRow,
		// MPI_Request *reqSendLastRow,
		// MPI_Request *reqRecvLastRow,
		// MPI_Request *reqRecvLastRow,
		MPI_Status *s) {
	int precNode = (rank == 0) ? numberOfNodes - 1 : rank - 1;
	int nextNode = (rank == (numberOfNodes -1)) ? 0 : rank + 1;

	MPI_Send(
		&old[1][1],
		bwidth, MPI_INT,
		precNode,
		FIRST_ROW_TAG,
		MPI_COMM_WORLD);
	MPI_Send(
		&old[bheight][1],
		bwidth, MPI_INT,
		nextNode,
		LAST_ROW_TAG,
		MPI_COMM_WORLD);

	MPI_Recv(
		&old[0][1],
		bwidth,	MPI_INT,
		precNode,
		LAST_ROW_TAG,
		MPI_COMM_WORLD, s
		);

	MPI_Recv(
		&old[bheight+1][1],
		bwidth,	MPI_INT,
		nextNode,
		FIRST_ROW_TAG,
		MPI_COMM_WORLD, s
		);
}

/* Take world wrap-around into account: */
void boundaryConditions(void) {
	// if (isMaster()) { // first node
	// } else if (isLast()) { // last node
	// } else // middle node

	for (int i = 0; i < bheight + 2; i++) {
		old[i][0] = old[i][bwidth];
		old[i][bwidth + 1] = old[i][1];
	}
}

void updateBoard(void) {
	int im, ip, jm, jp, nsum;

	// TODO: first and last row for last to wait the receive
	for (int i = 1; i <= bheight; i++)
		for (int j = 1; j <= bwidth; j++) {

			// sum surrounding cells
			im = i - 1;
			ip = i + 1;
			jm = j - 1;
			jp = j + 1;
			nsum = old[im][jp] + old[i][jp] + old[ip][jp]
					 + old[im][j]               + old[ip][j]
					 + old[im][jm] + old[i][jm] + old[ip][jm];

			switch (nsum) {
			case 3:
				// a new cell is born
				new[i][j] = 1;
				break;
			case 2:
				// nothing happens
				new[i][j] = old[i][j];
				break;
			default:
				// the cell, if any, dies
				new[i][j] = 0;
			}
		}

}

void updateState(void) {
	// wait here

	for (int i = 1; i <= bheight; i++)
		for (int j = 1; j <= bwidth; j++)
			old[i][j] = new[i][j];
}

void doTimeStep(int timestep) {
	// MPI_Request reqSendFirstRow, reqSendLastRow, reqRecvFirstRow, reqRecvLastRow;
	MPI_Status s;
	exchangeColumn(
		// &reqSendFirstRow,
		// &reqSendLastRow,
		// &reqRecvFirstRow,
		// &reqRecvLastRow,
		&s);

	boundaryConditions();

	updateBoard();

	updateState();

	if (print_world > 0 && (timestep % print_world) == (print_world - 1)) {
		printCells(timestep);
	}
}

void core(void) {
	struct timeval start, end;

	MPI_Barrier(MPI_COMM_WORLD);
	if (gettimeofday(&start, 0) != 0) {
		fprintf(stderr, "could not do timing\n");
		abort();
	}

	/*  time steps */
	for (int n = 0; n < nsteps; n++) {
		doTimeStep(n);
	}

	if (gettimeofday(&end, 0) != 0) {
		fprintf(stderr, "could not do timing\n");
		abort();
	}
	MPI_Barrier(MPI_COMM_WORLD);

	if (isMaster()) {
		// compute running time
		double rtime = (end.tv_sec + (end.tv_usec / 1000000.0)) -
								(start.tv_sec + (start.tv_usec / 1000000.0));

		/*  Iterations are done; sum the number of live cells */
		int isum = 0;
		for (int i = 1; i <= bheight; i++) {
			for (int j = 1; j <= bwidth; j++) {
				isum = isum + new[i][j];
			}
		}


		printf("\nNumber of live cells = %d\n", isum);
		fprintf(stderr, "Game of Life took %10.3f seconds\n", rtime);
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

	if (print_world > 0) {
		if (isMaster()) printf("\ninitial world:\n\n");
		printCells(-1);
	}

	core();

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	freeMatrix(old);
	freeMatrix(new);
	return 0;
}