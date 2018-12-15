#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

#define N_ROWS 8

void set_test(struct _matrix m)
{
	double tmp = 0;
	matrix_zero_up(m);
	for (unsigned i = 0; N_ROWS > i; i++)
		for (unsigned j = 0; N_ROWS > j; j++) {
			tmp = rand();
			if ((unsigned)tmp % 2)
				tmp *= -1;
			tmp /= (double)RAND_MAX;
			matrix_set_entry(m, i, j, tmp);
		}
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

	set_test(a);
	printf("\ntest matrix:\n");
	matrix_print(a);
	matrix_copy(c, a);

	matrix_lup_decompose(a);
	printf("\ntest matrix lu decompose (l and u in one, ones of l omitted on diagonal):\n");
	matrix_print(a);

	double tmp1[N_ROWS];
	double tmp2[N_ROWS];
	matrix_invert(a, b, tmp1, tmp2);
	printf("inverse of the test matrix\n");
	matrix_print(a);

	matrix_mult(d, a, c);
	printf("this should give the identity\n");
	matrix_print(d);

	return EXIT_SUCCESS;
}
