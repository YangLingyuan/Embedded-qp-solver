#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"

#define TEST_PRECISION 1e-6
#define NUM_TEST_RUNS 128
#define RAND_ENTRY_MIN -1e6
#define RAND_ENTRY_MAX 1e6

#define ABS(x) (0 > (x) ? -(x) : (x))

static unsigned almost_identity(struct _matrix m, double precision)
{
	double x = 0;
	for (unsigned i = 0; m.nrows > i; i++)
		for (unsigned j = 0; m.nrows > j; j++) {
			x = matrix_get_entry(m, i, j);
			if (isnan(x)) {
				fprintf(stderr,
					"entry (%u, %u) is nan\n", i, j);
				return 0;
			} else if (i == j && precision < ABS(1 - x)) {
				fprintf(stderr, "entry (%u, %u) is %.12lf "
						"but should be 1\n", i, j, x);
				return 0;
			} else if (i != j && precision < ABS(x)) {
				fprintf(stderr, "entry (%u, %u) is %.12lf "
						"but should be 0\n", i, j, x);
				return 0;
			}
		}
	return 1U;
}

int main(void)
{
	unsigned n = MATRIX_ROWS;
	unsigned m = MATRIX_ROWS;

	double a_entries[MATRIX_ROWS * MATRIX_ROWS];
	struct _matrix a;
	matrix_set_up(&a, n, m, a_entries);

	double b_entries[MATRIX_ROWS * MATRIX_ROWS];
	struct _matrix b;
	matrix_set_up(&b, n, m, b_entries);

	double c_entries[MATRIX_ROWS * MATRIX_ROWS];
	struct _matrix c;
	matrix_set_up(&c, n, m, c_entries);

	double d_entries[MATRIX_ROWS * MATRIX_ROWS];
	struct _matrix d;
	matrix_set_up(&d, n, m, d_entries);

	for (unsigned i = 0; NUM_TEST_RUNS > i; i++) {
		matirx_random_pos_def(a, b, c, RAND_ENTRY_MIN, RAND_ENTRY_MAX);
		matrix_copy(c, a);

		matrix_invert(a, b);

		matrix_mult(d, a, c);
		if (!almost_identity(d, TEST_PRECISION))
			printf("error\n");
	}

	return EXIT_SUCCESS;
}
