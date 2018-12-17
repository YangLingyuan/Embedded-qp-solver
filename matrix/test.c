#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix.h"

#define N_ROWS 500

#define TEST_PRECISION 1e-6
#define NUM_TEST_RUNS 12
#define RAND_ENTRY_MIN -1e3
#define RAND_ENTRY_MAX 1e3

#define ABS(x) (x < 0 ? -x : x)

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
	unsigned n = N_ROWS;
	unsigned m = N_ROWS;

	double a_entries[N_ROWS * N_ROWS];
	unsigned a_permutation[N_ROWS];
	struct _matrix a;
	matrix_set_up(&a, n, m, a_entries, a_permutation);

	double b_entries[N_ROWS * N_ROWS];
	unsigned b_permutation[N_ROWS];
	struct _matrix b;
	matrix_set_up(&b, n, m, b_entries, b_permutation);

	double c_entries[N_ROWS * N_ROWS];
	unsigned c_permutation[N_ROWS];
	struct _matrix c;
	matrix_set_up(&c, n, m, c_entries, c_permutation);

	double d_entries[N_ROWS * N_ROWS];
	unsigned d_permutation[N_ROWS];
	struct _matrix d;
	matrix_set_up(&d, n, m, d_entries, d_permutation);

	for (unsigned i = 0; NUM_TEST_RUNS > i; i++) {
		matirx_random_pos_def(a, b, c, RAND_ENTRY_MIN, RAND_ENTRY_MAX);
		matrix_copy(c, a);

		double tmp1[N_ROWS];
		double tmp2[N_ROWS];
		matrix_invert(a, b, tmp1, tmp2);

		matrix_mult(d, a, c);
		if (!almost_identity(d, TEST_PRECISION))
			printf("error\n");
	}

	return EXIT_SUCCESS;
}
