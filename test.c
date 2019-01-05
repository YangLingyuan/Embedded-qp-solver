#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix_ops.h"
#include "qp.h"

//#define DK_ALPHA_INSPECT
#define LHS_RHS_IN_ARMIJO

#define TEST_PRECISION 1e-9
#define NUM_INVERSION_TEST_RUNS 128
#define RAND_ENTRY_MIN -1e1
#define RAND_ENTRY_MAX 1e1

#define ITERATIONS 1e4

#define ABS(x) (x < 0 ? -x : x)

static unsigned almost_identity(struct _matrix * m, double precision)
{
	double x = 0;
	unsigned nrows = m->dimensions[ROW];
	for (unsigned i = 0; nrows > i; i++)
		for (unsigned j = 0; nrows > j; j++) {
			x = matrix_get_entry(m, ME(i, j));
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

unsigned armijo(struct _quadratic_form * qf, struct _matrix * x,
  		      struct _matrix * dk, double alpha, double c2)
{
	double lhs = 0;
	double rhs = 0;
	struct _matrix * tmp = matrix_alloc(Nx1);
	struct _matrix * dk_alpha = matrix_alloc(Nx1);
	matrix_copy(dk_alpha, dk);
	matrix_scalar_mult(dk_alpha, alpha);
#ifdef DK_ALPHA_INSPECT
	matrix_print(dk_alpha);
#endif
	matrix_add(tmp, x, dk_alpha);

	lhs = qf->eval(qf, tmp);

	struct _matrix * qf_grad = qf->eval_grad(qf, x);
	rhs += c2 * alpha * matrix_scalar_prod(qf_grad, dk);
	rhs += qf->eval(qf, x);

	matrix_free(qf_grad);
	matrix_free(dk_alpha);
	matrix_free(tmp);

	return lhs <= rhs;
}

double
line_search(struct _quadratic_form * qf,
            struct _matrix * x, struct _matrix * dk, double c1, double c2)
{
	double alpha = 1;
	while (!armijo(qf, x, dk, alpha, c1))
		alpha *= c2;
	return alpha;
}

void inversion_test(void)
{
	for (unsigned i = 0; NUM_INVERSION_TEST_RUNS > i; i++) {
		struct _matrix * a = matrix_alloc(NxN);
		struct _matrix * b = matrix_alloc(NxN);
		struct _matrix * c = matrix_alloc(NxN);

		matirx_random_pos_def(a, RAND_ENTRY_MIN, RAND_ENTRY_MAX);
		matrix_copy(b, a);
		matrix_invert(a);
		matrix_mult(c, b, a);
		if (!almost_identity(c, TEST_PRECISION))
			printf("error\n");
		matrix_free(a);
		matrix_free(b);
		matrix_free(c);
	}
}

void mult_test(void)
{
	struct _matrix * a = matrix_alloc(NxN);
	for (unsigned i = 0; N > i; i++)
		matrix_set_entry(a, ME(i, 0), i);

	struct _matrix * b = matrix_alloc(Nx1);
	matrix_set_entry(b, ME(0, 0), 1);

	struct _matrix * c = matrix_alloc(Nx1);
	matrix_mult(c, a, b);

	struct _matrix * d = matrix_alloc(Nx1);
	for (unsigned i = 0; N > i; i++)
		matrix_set_entry(d, ME(i, 0), i);
	matrix_trans(c);
	matrix_trans(d);

	if (14 != matrix_scalar_prod(d, c))
		fprintf(stderr, "error in mult test\n");

	matrix_free(a);
	matrix_free(b);
	matrix_free(c);
	matrix_free(d);
}

int main(void)
{
	/* important */
	kmalloc_init();

	inversion_test();
	mult_test();

	/* newton method */
	/* prepare quadratic cost function */
	struct _matrix * p = matrix_alloc(NxN);
	matirx_random_pos_def(p, RAND_ENTRY_MIN, RAND_ENTRY_MAX);

	struct _matrix * q = matrix_alloc(Nx1);
	matrix_zero_up(q);

	double r = 0;

	struct _quadratic_form * qf = quadratic_form_alloc(p, q, r);

	struct _matrix * tmp = 0;
	struct _matrix * dk = matrix_alloc(Nx1);
	struct _matrix * h2_inv = matrix_alloc(NxN);
	matrix_mult(h2_inv, qf->p, qf->p);
	matrix_invert(h2_inv);
	struct _matrix * x = matrix_alloc(Nx1);
	for (unsigned i = 0; N > i; i++)
		matrix_set_entry(x, ME(i, 0), 1);
	double alpha = 1;
	double c1 = 1e-4;
	double c2 = 0.9;

	for (unsigned i = 0; ITERATIONS > i; i++) {
		tmp = qf->eval_grad(qf, x);
		matrix_mult(dk, h2_inv, tmp);
		matrix_free(tmp);
		alpha = line_search(qf, x, dk, c1, c2);
		matrix_scalar_mult(dk, alpha);
		matrix_sub(x, x, dk);
	}
	matrix_print(x);

	matrix_free(h2_inv);
	matrix_free(dk);
	matrix_free(p);
	matrix_free(q);
	matrix_free(x);
	quadratic_form_free(qf);

	return EXIT_SUCCESS;
}
