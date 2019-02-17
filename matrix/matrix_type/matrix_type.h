#ifndef MATRIX_TYPE_H
#define MATRIX_TYPE_H

/* configure matrix types here */
#define N 48U
#define NxN_MAX 5
#define Nx1_MAX 10

/* since struct _matrix_entry is convenient from an implementation point */
/* of view but cumbersome when using the api, ME stands for matrix entry */
#define ME(i, j) ((struct _matrix_entry) {(i), (j)})

#define MATRIX_GET_ROW(m) (GETM32(16, 31, m->dimensions) >> 16)
#define MATRIX_GET_COL(m) GETM32(0, 15, m->dimensions)
#define MATRIX_SET_ROW(m, n) SETM32(16, 31, m->dimensions, n)
#define MATRIX_SET_COL(m, n) SETM32(0, 15, m->dimensions, n)

/* structures related to matrices */
struct _matrix {
	unsigned dimensions;
	double * elements;
};

struct _matrix_entry {
	unsigned row;
	unsigned col;
};

#endif
