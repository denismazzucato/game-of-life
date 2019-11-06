#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

static int
	FIRST_ROW_TAG = 1, // tag for communication over first row
	LAST_ROW_TAG = 2, // tag for communication over last row
	ROOT = 0, // master rank
	NO_TIMESTEP = -1; // to print cells without initial line

static int rank, number_nodes; // initializated by mpi_init

static int *counts, *displs, sum; // used in gather and scatter operations

static int bwidth, // board width
					 bheight, // board height, in the seq algorithm bheight == nrows
					 nrows, // entire board height, not only the part of this node
					 nsteps; // time limit

// matrices (bwidth+2)x(bheight+2), to wrap-around the world structure
static int **old = 0, **new = 0;

static int print_world = 0;

// use fixed world or random world?
// static int random_world = 0; // fixed
static int random_world = 1; // random

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

// close the MPI communication and exit
void abort(void) {
	MPI_Abort(MPI_COMM_WORLD, 1);
	MPI_Finalize();
	exit(1);
}

// these two functions are required to scatter and gather data
void build_matrix_from_buffer(int* recbuf) {
	for (int i = 1; i <= bheight; i++)
		for (int j = 1; j <= bwidth; j++)
			old[i][j] = recbuf[(i-1)*bwidth + (j-1)];
}
int* build_buffer_from_old() {
	int *buffer = malloc(bheight * bwidth * sizeof(int));

	for (int i = 1; i <= bheight; i++)
		for (int j = 1; j <= bwidth; j++)
			buffer[(i-1)*bwidth + (j-1)] = old[i][j];

	return buffer;
}

// function that avoid the explicit use of rank
int is_master(void) {
	return rank == ROOT;
}
int is_last(void) {
	return rank == (number_nodes - 1);
}
int next_node(void) {
	return is_last() ? 0 : rank + 1;
}
int prec_node(void) {
	return is_master() ? number_nodes - 1 : rank - 1;
}

// gather to print cells
void print_cells(int timestep) {
	int *global = malloc(bwidth * nrows * sizeof(int)), // global data in master
			*data = build_buffer_from_old(); // local data

	MPI_Gatherv(
		data, bheight*bwidth, MPI_INT,
		global, counts, displs, MPI_INT,
		ROOT, MPI_COMM_WORLD);

	if (is_master()) {
		if (timestep != NO_TIMESTEP)
			fprintf(stdout, "\nafter time step %d:\n\n", timestep);

		for (int i = 0; i < nrows; i++) {
			for (int j = 0; j < bwidth; j++)
				fprintf(stdout,"%s", (global[i * bwidth + j]) ? "0" : " ");
			fprintf(stdout, "\n");
		}
		free(global);
	}
	free(data);
}

// I've already check that argv contains exactly 5 args (name + args)
void init_parameters(char *argv[]) {
	nrows = atoi(argv[1]);
	bwidth = atoi(argv[2]);
	nsteps = atoi(argv[3]);
	print_world = atoi(argv[4]);
}

void init_MPI(void) {
	MPI_Init(NULL, NULL);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &number_nodes);
}

/* the array allocation can be done without communication, each processors
 * already has all the information to instantiate the old and new matrices
 *
 * the rows that exceed from splitting in perfect halves the data will be
 * distribuited in the last nodes
 */
void allocateArrays(void) {
	int other_nodes = nrows % number_nodes; // nodes with N/p +1 rows
	int full_nodes = number_nodes - other_nodes; // nodes with N/p rows

	// this is different on every node
	int number_of_real_rows = (rank < full_nodes
		? nrows / number_nodes + 2
			// each node in the first full_nodes holds N/p rows
		: nrows / number_nodes + 3);
			// each node in the last other_nodes holds (N/p)+1 rows
	int numberOfRealColumns = bwidth + 2;

	old = malloc(number_of_real_rows * sizeof(int*));
	new = malloc(number_of_real_rows * sizeof(int*));

	for (int i = 0; i < number_of_real_rows; i++) {
		old[i] = malloc(numberOfRealColumns * sizeof(int*));
		new[i] = malloc(numberOfRealColumns * sizeof(int*));
	}

	// now we have to set the correct value for bheight
	// because then I want to iterate over bheight
	bheight = number_of_real_rows - 2;

	// now allcoate counts and displs for gather and scatter operations
	counts = malloc(number_nodes * sizeof(int)),
	displs = malloc(number_nodes * sizeof(int));
	sum = 0;
	int iter;

  for (iter = 0; iter < number_nodes-other_nodes; iter++)
    counts[iter] = bwidth*(nrows/number_nodes);
  for (;iter < number_nodes; iter++)
    counts[iter] = bwidth*(nrows/number_nodes + 1);
  for (int i = 0; i < number_nodes; i++) {
    displs[i] = sum;
    sum += counts[i];
  }
}

float randomNumber(void) {
	return rand() / ((float)RAND_MAX + 1);
}

int isInsideStartWorld(int i, int j) {
	return (i < sizeof(start_world) / sizeof(char *)) &&
				 (j < strlen(start_world[i]));
}

/* this function has to return a linear representation of the initial matrix
 */
int* build_board(void) {
	int *tmp_board = malloc((nrows*bwidth) * sizeof(int));
	srand(1); // to fix rand seme
	for (int i = 0; i < nrows; i++) {
		for (int j = 0; j < bwidth; j++) {
			if (random_world) {
				float x = randomNumber();
				tmp_board[i*bwidth + j] = (x < 0.5 ? 0 : 1);
			} else if (isInsideStartWorld(i, j))
				tmp_board[i*bwidth + j] = (start_world[i][j] != '.');
			else
				tmp_board[i*bwidth + j] = 0;
		}
	}

	return tmp_board;
}

void scatter_rows(int* data) {
	int	*recbuf = malloc(bheight*bwidth * sizeof(int));

	// divide the data among processes as described by counts and displs
  MPI_Scatterv(
		data, counts, displs, MPI_INT,
		recbuf, bheight*bwidth, MPI_INT,
		ROOT, MPI_COMM_WORLD);

	build_matrix_from_buffer(recbuf);

	free(recbuf);
}

/* this function built all the line,
 * when (~)N/P line are created are immediately sent to the corresponding node
 *
 * this is because I don't want to store all the data in a single node even
 * only for the initialization because this table may be too large to be stored
 */
void initialize_board(void) {
	int *temp_board = 0;

	if (is_master()) temp_board = build_board();

	scatter_rows(temp_board);

	if (temp_board) free(temp_board);
}

void free_matrices(void) {
	if(old) {
		for (int i = 0; i < bheight + 2; i++)
			free(old[i]);
	}
	free(old);
	if(new) {
		for (int i = 0; i < bheight + 2; i++)
			free(new[i]);
	}
	free(new);
	if(counts) free(counts);
	if(displs) free(displs);
}

void args_not_provided(void) {
	if (is_master()) {
		fprintf(
			stderr,
			"Usage: board_height board_width steps_count print_iter\n");
		abort();
	}
}

// our goal is to solve large problem fast,
// so requiring at least number_nodes rows isn't a problem
void check_args(void) {
	if (nrows <= number_nodes && is_master()) {
		fprintf(
			stderr,
			"Too few rows provided, provide at least %d rows\n",
			(number_nodes + 1));
		abort();
	}
}

void exchange_column(void) {
	MPI_Status s;

	// this provide a buffer to MPI
	int bsize, *buffer, length;
	MPI_Pack_size(bwidth, MPI_INT, MPI_COMM_WORLD, &length);
	bsize = 2 * (MPI_BSEND_OVERHEAD + length);
	buffer = malloc(bsize);
	MPI_Buffer_attach(buffer, bsize);

	// buffered send to avoid unsafe program
	MPI_Bsend(
		&old[1][1],	bwidth, MPI_INT,
		prec_node(), FIRST_ROW_TAG,
		MPI_COMM_WORLD);
	MPI_Bsend(
		&old[bheight][1],	bwidth, MPI_INT,
		next_node(), LAST_ROW_TAG,
		MPI_COMM_WORLD);

	MPI_Recv(
		&old[0][1],	bwidth,	MPI_INT,
		prec_node(), LAST_ROW_TAG,
		MPI_COMM_WORLD, &s);
	MPI_Recv(
		&old[bheight+1][1],	bwidth,	MPI_INT,
		next_node(), FIRST_ROW_TAG,
		MPI_COMM_WORLD, &s);

	MPI_Buffer_detach( &buffer, &bsize );
	free(buffer);
}

/* Take world wrap-around into account: */
void boundary_conditions(void) {
	for (int i = 0; i < bheight + 2; i++) {
		old[i][0] = old[i][bwidth];
		old[i][bwidth + 1] = old[i][1];
	}
}

void update_board(void) {
	int im, ip, jm, jp, nsum;

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

void update_state(void) {
	for (int i = 1; i <= bheight; i++)
		for (int j = 1; j <= bwidth; j++)
			old[i][j] = new[i][j];
}

void do_timestep(int timestep) {
	exchange_column();
	boundary_conditions();
	update_board();
	update_state();

	if (print_world > 0 && (timestep % print_world) == (print_world - 1))
		print_cells(timestep);
}

void core(void) {
	struct timeval start, end;

	if (gettimeofday(&start, 0) != 0) {
		fprintf(stderr, "could not do timing\n");
		abort();
	}

	// compute all timesteps
	for (int n = 0; n < nsteps; n++)
		do_timestep(n);

	if (gettimeofday(&end, 0) != 0) {
		fprintf(stderr, "could not do timing\n");
		abort();
	}

	// compute running time
	double
		rtime = (end.tv_sec + (end.tv_usec / 1000000.0)) -
						(start.tv_sec + (start.tv_usec / 1000000.0)),
		maxrtime;
	MPI_Reduce(&rtime, &maxrtime, 1, MPI_DOUBLE, MPI_MAX, ROOT, MPI_COMM_WORLD);

	/*  Iterations are done; sum the number of live cells */
	int isum = 0, totalisum;
	for (int i = 1; i <= bheight; i++)
		for (int j = 1; j <= bwidth; j++)
			isum = isum + old[i][j];
	MPI_Reduce(&isum, &totalisum, 1, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);

	if (is_master()) {
		printf("Number of live cells = %d\n", totalisum);
		fprintf(stderr, "Game of Life took %10.3f seconds\n", maxrtime);
	}
}

int main (int argc, char *argv[]) {

	init_MPI();

	if (argc == 5) {
		init_parameters(argv);
		check_args();
	} else
		args_not_provided();
	// bwidth, nrows, nsteps, print_world are now initialized and correct

	allocateArrays();
	// bheight is now initialized

	initialize_board();

	if (print_world > 0) {
		if (is_master()) printf("\ninitial world:\n\n");
		print_cells(NO_TIMESTEP);
	}

	// computation core
	core();

	MPI_Finalize();
	free_matrices();
	return 0;
}