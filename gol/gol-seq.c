/***********************

Conway's Game of Life

Based on https://web.cs.dal.ca/~arc/teaching/CS4125/2014winter/Assignment2/Assignment2.html

************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

static int bwidth, bheight;
static int **old, **new;

static int print_world = 0;

// use fixed world or random world?
// static int random_world = 0;
static int random_world = 1;

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

void
printCells(void)
{
    int i, j;

    for (i = 1; i <= bwidth; i++) {
        for (j = 1; j <= bheight; j++) {
            if (old[i][j]) {
                printf("O");
            } else {
                printf(" ");
            }
        }
        printf("\n");
    }
}

// update board for step n
void
doTimeStep(int n)
{
    int i, j;

    /* Take world wrap-around into account: */

    /* corner boundary conditions */
    old[0][0] = old[bwidth][bheight];
    old[0][bheight + 1] = old[bwidth][1];
    old[bwidth + 1][bheight + 1] = old[1][1];
    old[bwidth + 1][0] = old[1][bheight];

    /* left-right boundary conditions */
    for (i = 1; i <= bwidth; i++) {
        old[i][0] = old[i][bheight];
        old[i][bheight + 1] = old[i][1];
    }

    /* top-bottom boundary conditions */
    for (j = 1; j <= bheight; j++) {
        old[0][j] = old[bwidth][j];
        old[bwidth + 1][j] = old[1][j];
    }

    // update board
    for (i = 1; i <= bwidth; i++) {
        for (j = 1; j <= bheight; j++) {
            int im, ip, jm, jp, nsum;

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

    /* copy new state into old state */
    for (i = 1; i <= bwidth; i++) {
        for (j = 1; j <= bheight; j++) {
            old[i][j] = new[i][j];
        }
    }

    if (print_world > 0 && (n % print_world) == (print_world - 1)) {
        printf("\nafter time step %d:\n\n", n);
        printCells();
    }
}

int
main(int argc, char *argv[])
{
    int i, j, n;
    int ni, nj;
    int isum;
    int nsteps;
    double rtime;
    struct timeval start, end;

    /* Get Parameters */
    if (argc != 5) {
        fprintf(stderr,
            "Usage: %s board_width board_height steps_count print_iter\n",
            argv[0]);
        exit(1);
    }
    bwidth = atoi(argv[1]);
    bheight = atoi(argv[2]);
    nsteps = atoi(argv[3]);
    print_world = atoi(argv[4]);

    /* allocate arrays */
    ni = bwidth + 2;  /* add 2 for left and right ghost cells */
    nj = bheight + 2;
    old = malloc(ni * sizeof(int*));
    new = malloc(ni * sizeof(int*));

    for (i = 0; i < ni; i++) {
        old[i] = malloc(nj * sizeof(int));
        new[i] = malloc(nj * sizeof(int));
    }

    /*  initialize board */
    if (random_world) {
        srand(1);
        for (i = 1; i <= bwidth; i++) {
            for (j = 1; j <= bheight; j++) {
                float x = rand() / ((float)RAND_MAX + 1);
                if (x < 0.5) {
                    old[i][j] = 0;
                } else {
                    old[i][j] = 1;
                }
            }
        }
    } else {
        /* predefined start_world */
        for (i = 1; i <= bwidth; i++) {
            for (j = 1; j <= bheight; j++) {
                if ((i <= sizeof(start_world) / sizeof(char *)) &&
                    (j <= strlen(start_world[i - 1]))) {
                    old[i][j] = (start_world[i - 1][j - 1] != '.');
                } else {
                    old[i][j] = 0;
                }
            }
        }
    }

    if (print_world > 0) {
        printf("\ninitial world:\n\n");
        printCells();
    }

    if (gettimeofday(&start, 0) != 0) {
        fprintf(stderr, "could not do timing\n");
        exit(1);
    }

    /*  time steps */
    for (n = 0; n < nsteps; n++) {
        doTimeStep(n);
    }

    if (gettimeofday(&end, 0) != 0) {
        fprintf(stderr, "could not do timing\n");
        exit(1);
    }

    // compute running time
    rtime = (end.tv_sec + (end.tv_usec / 1000000.0)) -
        (start.tv_sec + (start.tv_usec / 1000000.0));

    /*  Iterations are done; sum the number of live cells */
    isum = 0;
    for (i = 1; i <= bwidth; i++) {
        for (j = 1; j <= bheight; j++) {
            isum = isum + new[i][j];
        }
    }

    printf("Number of live cells = %d\n", isum);
    fprintf(stderr, "Game of Life took %10.3f seconds\n", rtime);

    return 0;
}
