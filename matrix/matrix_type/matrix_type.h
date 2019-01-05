#ifndef MATRIX_TYPE_H
#define MATRIX_TYPE_H

/* general */
#define ROW 0
#define COLUMN 1
#define NUMBER_DIMENSIONS 2

/* configure matrix types here */
#define N 2
#define NxN_MAX 16
#define Nx1_MAX 16

/* since struct _matrix_entry is convenient from an implementation point */
/* of view but cumbersome when using the api, ME stands for matrix entry */
#define ME(i, j) ((struct _matrix_entry) {(i), (j)})

/* structures related to matrices */
struct _matrix {
	unsigned dimensions[NUMBER_DIMENSIONS];
	double * elements;
};

struct _matrix_entry {
	unsigned row;
	unsigned col;
};

#endif
