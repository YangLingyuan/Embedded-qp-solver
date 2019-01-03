#include <stdio.h>
#include <stdlib.h>

#include "matrix.h"

#define ADD 1U
#define SUB 0U

#define ABS(x) (0 > (x) ? -(x) : (x))

#define NERVOUS

struct _pivot {
	double pivot;
	unsigned pivot_idx;
};

double __inline
matrix_get_entry(struct _matrix m, unsigned row, unsigned col)
{
#ifdef NERVOUS
	if (!m.elements) {
		fprintf(stderr, "invalid matrix in matrix_get_entry\n");
		return 0;
	}
#endif
	return *(m.elements + row * m.ncols + col);
}

void __inline
matrix_set_entry(struct _matrix m, unsigned row, unsigned col, double val)
{
#ifdef NERVOUS
	if (!m.elements) {
		fprintf(stderr, "invalid matrix in matrix_set_entry\n");
		return;
	}
#endif
	*(m.elements + row * m.ncols + col) = val;
}

void __inline
matrix_swap_entries(struct _matrix m,
		unsigned row1, unsigned col1, unsigned row2, unsigned col2) 
{
	double a1 = matrix_get_entry(m, row1, col1);
	double a2 = matrix_get_entry(m, row2, col2);
	matrix_set_entry(m, row1, col1, a2);
	matrix_set_entry(m, row2, col2, a1);
}

void __inline
matrix_negate_entry(struct _matrix m, unsigned row, unsigned col)
{
	double tmp = -matrix_get_entry(m, row, col);
	matrix_set_entry(m, row, col, tmp);
}

void
matrix_set_up(struct _matrix * m, unsigned nrows,
                                  unsigned ncols, double * elements)
{
#ifdef NERVOUS
	if (!m) {
		fprintf(stderr, "invalid matrix "
				"pointer in matrix_set_up\n");
		return;
	}
	if (!elements) {
		fprintf(stderr, "invalid matrix element "
				"array pointer in matrix_set_up\n");
		return;
	}
#endif

	m->nrows = nrows;
	m->ncols = ncols;
	m->elements = elements;
}

static unsigned __inline
matrix_mult_dimension_compatible(struct _matrix prod,
				struct _matrix a, struct _matrix b)
{
	if (prod.nrows != a.nrows)
		return 0;
	if (a.ncols != b.nrows)
		return 0;
	if (prod.ncols != b.ncols)
		return 0;
	return 1;
}

static unsigned __inline
matrix_sum_dimension_compatible(struct _matrix sum,
				struct _matrix a, struct _matrix b)
{
	unsigned rows_ok = sum.nrows == a.nrows && a.nrows == b.nrows;
	unsigned cols_ok = sum.ncols == a.ncols && a.ncols == b.ncols;
	if (rows_ok && cols_ok)
		return 1;
	return 0;
}

void matrix_zero_up(struct _matrix m)
{
	unsigned nrows = m.nrows;
	unsigned ncols = m.ncols;
	for (unsigned i = 0; nrows > i; i++)
		for (unsigned j = 0; ncols > j; j++)
			matrix_set_entry(m, i, j, 0);
}

void matrix_copy(struct _matrix a, struct _matrix b)
{
#ifdef NERVOUS
	if (!matrix_sum_dimension_compatible(a, a, b)) {
		fprintf(stderr, "incompatible matrix "
				"dimensions in matrix_copy\n");
		return;
	}
#endif

	double tmp = 0;
	unsigned nrows = a.nrows;
	unsigned ncols = a.ncols;
	for (unsigned i = 0; nrows > i; i++)
		for (unsigned j = 0; ncols > j; j++) {
			tmp = matrix_get_entry(b, i, j);
			matrix_set_entry(a, i, j, tmp);
		}
}

void matrix_mult(struct _matrix prod,
		 struct _matrix a, struct _matrix b)
{
#ifdef NERVOUS
	if (!matrix_mult_dimension_compatible(prod, a, b)) {
		fprintf(stderr, "incompatible matrix "
				"dimensions in matrix_mult\n");
		return;
	}
#endif

	unsigned n = prod.nrows;
	unsigned m = prod.ncols;
	for (unsigned i = 0; n > i; i++)
		for (unsigned j = 0; m > j; j++) {
			double acc = 0;
			for (unsigned k = 0; m > k; k++)
				acc += matrix_get_entry(a, i, k)
				       * matrix_get_entry(b, k, j);
			matrix_set_entry(prod, i, j, acc);
		}
}

static void __inline
matrix_add_sub(struct _matrix sum,
		struct _matrix a, struct _matrix b, unsigned add_sub)
{
#ifdef NERVOUS
	if (!matrix_sum_dimension_compatible(sum, a, b)) {
		fprintf(stderr, "incompatible matrix dimensions in "
				"matrix_%s\n", ADD == add_sub ? "add" : "sub");
		return;
	}
#endif

	unsigned n = sum.nrows;
	unsigned m = sum.ncols;
	double aa = 0;
	double bb = 0;
	for (unsigned i = 0; n > i; i++)
		for (unsigned j = 0; m > j; j++) {
			aa = matrix_get_entry(a, i, j);
			bb = matrix_get_entry(b, i, j);
			matrix_set_entry(sum, i, j,
					ADD == add_sub ? aa + bb : aa - bb);
		}
}

void matrix_add(struct _matrix sum,
		struct _matrix a, struct _matrix b)
{
	matrix_add_sub(sum, a, b, ADD);
}

void matrix_sub(struct _matrix sum,
		struct _matrix a, struct _matrix b)
{
	matrix_add_sub(sum, a, b, SUB);
}

void matrix_trans(struct _matrix m)
{
	unsigned nrows = m.nrows;
	unsigned ncols = m.ncols;
	for (unsigned i = 0; nrows > i; i++)
		for (unsigned j = i + 1; ncols > j; j++)
			matrix_swap_entries(m, i, j, j, i);
}

void matrix_neg(struct _matrix m)
{
	unsigned nrows = m.nrows;
	unsigned ncols = m.ncols;
	for (unsigned i = 0; nrows > i; i++)
		for (unsigned j = 0; ncols > j; j++)
			matrix_negate_entry(m, i, j);
}

static void __inline
permutation_swap(unsigned * permutation, unsigned i, unsigned j)
{
#ifdef NERVOUS
	if (!permutation) {
		fprintf(stderr, "invalid permutation in permutation_swap");
		return;
	}
#endif

	double tmp = permutation[i];
	permutation[i] = permutation[j];
	permutation[j] = tmp;
}

static struct _pivot
matrix_lup_pivot(struct _matrix m, unsigned k)
{
	double pivot = 0;
	unsigned pivot_idx = 0;

	for (unsigned i = k; m.nrows > i; i++) {
		double tmp = matrix_get_entry(m, i, k);
		tmp = ABS(tmp);
		if(tmp > pivot) {
			pivot = tmp;
			pivot_idx = i;
		}
	}

	return (struct _pivot){ .pivot = pivot, .pivot_idx = pivot_idx, };
}

static void matrix_lup_decompose(struct _matrix m, unsigned * permutation)
{
#ifdef NERVOUS
	if (!permutation) {
		fprintf(stderr, "invalid permutation in matrix_lup_decompose");
		return;
	}
#endif

	for (unsigned k = 0; m.nrows > k + 1; k++) {

		struct _pivot pivot = matrix_lup_pivot(m, k);
		if(!pivot.pivot) {
			fprintf(stderr, "singular matrix "
					"in matrix_lup_decompose\n");
			return;
		}

		permutation_swap(permutation, pivot.pivot_idx, k);
		
		for (unsigned i = 0; m.nrows > i; i++)
			matrix_swap_entries(m, pivot.pivot_idx, i, k, i);

		for (unsigned i = k + 1; m.nrows > i; i++) {
			double tmp = matrix_get_entry(m, i, k);
			tmp /= matrix_get_entry(m, k, k);
			matrix_set_entry(m, i, k, tmp);
			for (unsigned j = k + 1; m.nrows > j; j++) {
				tmp = matrix_get_entry(m, i, j);
				tmp -= matrix_get_entry(m, i, k)
				       * matrix_get_entry(m, k, j);
				matrix_set_entry(m, i, j, tmp);
			}
		}
	}
	return;
}

void matrix_invert(struct _matrix m, struct _matrix tmp)
{
#ifdef NERVOUS
	if (!(m.nrows <= tmp.nrows) || !(m.ncols <= tmp.ncols)) {
		fprintf(stderr, "invalid matrix tmporary "
				"matrix pair in matrix_lup_decompose");
		return;
	}
#endif

	unsigned nrows = m.nrows;
	double v[MATRIX_ROWS];
	double w[MATRIX_ROWS];
	for (unsigned n = 0; nrows > n; n++)
		v[n] = w[n] = 0;
	unsigned permutation[MATRIX_ROWS];
	for (unsigned i = 0; m.nrows > i; i++)
		permutation[i] = i;

	matrix_lup_decompose(m, permutation);

	for (unsigned i = 0; nrows > i; i++) {
		for(unsigned j = 0; nrows > j; j++)
			matrix_set_entry(tmp, i, j, 0);
		matrix_set_entry(tmp, i, i, 1);

		for (unsigned n = 0; nrows > n; n++) {
			double t = 0;
			for (unsigned k = 0; n >= k + 1; k++)
				t += matrix_get_entry(m, n, k) * w[k];
			w[n] = matrix_get_entry(tmp, i, permutation[n]) - t;
		}

		for (unsigned n = nrows - 1; 1 <= n + 1; n--) {
			double t = 0;
			for(unsigned k = n+1; nrows > k; k++)
				t += matrix_get_entry(m, n, k) * v[k];
			v[n] = (w[n] - t) / matrix_get_entry(m, n, n);
		}

		for (unsigned j = 0; nrows > j; j++)
			matrix_set_entry(tmp, i, j, v[j]);
	}

	matrix_trans(tmp);
	matrix_copy(m, tmp);

	return;
}

void matrix_print(struct _matrix m)
{
	unsigned nrows = m.nrows;
	unsigned ncols = m.ncols;
	for (unsigned i = 0; nrows > i; i++) {
		for (unsigned j = 0; ncols > j; j++)
			printf("%lf ", matrix_get_entry(m, i, j));
		printf("\n");
	}
	printf("\n");
}

static double random_number(double min, double max)
{
	double r = rand();
	return min + r * (max - min) / (double)RAND_MAX;
}

void matrix_random(struct _matrix m, double min, double max)
{
	unsigned nrows = m.nrows;
	unsigned ncols = m.ncols;
	for (unsigned i = 0; nrows > i; i++)
		for (unsigned j = 0; ncols > j; j++)
			matrix_set_entry(m, i, j, random_number(min, max));
}

void matirx_random_pos_def(struct _matrix a,
		struct _matrix b, struct _matrix c, double min, double max)
{
	/* TODO check rank */
	matrix_random(b, min, max);
	matrix_copy(c, b);
	matrix_trans(b);
	matrix_mult(a, b, c);
}
