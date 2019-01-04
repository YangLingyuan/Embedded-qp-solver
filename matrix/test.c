#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "matrix_ops.h"

#define TEST_PRECISION 1e-6
#define NUM_TEST_RUNS 128
#define RAND_ENTRY_MIN -1e1
#define RAND_ENTRY_MAX 1e1

#define ITERATIONS 1e4

#define ABS(x) (x < 0 ? -x : x)

struct _quadratic_form {
	struct _matrix * p;
	struct _matrix * q;
	double r;
	double (*eval)(struct _quadratic_form * qf, struct _matrix * x);
	struct _matrix * (*eval_grad)(struct _quadratic_form * qf,
			                           struct _matrix * x);
};

static double
quadratic_form_eval(struct _quadratic_form * qf, struct _matrix * x)
{
	struct _matrix * tmp = matrix_alloc(Nx1);
	if (!tmp) {
		fprintf(stderr, "no matrix available "
				"to alloc in quadratic_form_eval\n");
		return 0;
	}
	matrix_mult(tmp, qf->p, x);
	double a = 0.5 * matrix_scalar_prod(x , tmp);
	kfree(tmp, Nx1);

	a += matrix_scalar_prod(qf->q, x);
	a += qf->r;

	return a;
}

static struct _matrix *
quadratic_form_eval_grad(struct _quadratic_form * qf, struct _matrix * x)
{
	struct _matrix * ret = matrix_alloc(Nx1);
	if (!ret) {
		fprintf(stderr, "no matrix available "
				"to alloc in quadratic_form_eval_grad\n");
		return 0;
	}
	matrix_mult(ret, qf->p, x);
	matrix_add(ret, ret, qf->q);

	return ret;
}

void quadratic_form_init(struct _quadratic_form * qf,
		         struct _matrix * p, struct _matrix * q, double r)
{
	qf->p = p;
	qf->q = q;
	qf->r = r;
	qf->eval = quadratic_form_eval;
	qf->eval_grad = quadratic_form_eval_grad;
}

static unsigned almost_identity(struct _matrix * m, double precision)
{
	double x = 0;
	unsigned nrows = m->dimensions[ROW];
	for (unsigned i = 0; nrows > i; i++)
		for (unsigned j = 0; nrows > j; j++) {
			x = matrix_get_entry(m, (struct _matrix_entry) {i, j});
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
	matrix_add(tmp, x, dk_alpha);

	lhs = qf->eval(qf, tmp);

	struct _matrix * qf_grad = qf->eval_grad(qf, x);
	rhs += c2 * alpha * matrix_scalar_prod(qf_grad, dk);
	rhs += qf->eval(qf, x);

	kfree(qf_grad, Nx1);
	kfree(dk_alpha, Nx1);
	kfree(tmp, Nx1);

	return lhs <= rhs;
}

double
line_search(struct _quadratic_form * qf,
            struct _matrix * x, struct _matrix * dk, double c1, double c2)
{
	double alpha = 1;
	while (!armijo(qf, x, dk, alpha, c1))
		alpha *= c2;
	if (1e-9 > alpha)
		return 1e-2;
	fprintf(stderr, "no zero alpha\n");
	return alpha;
}

int main(void)
{
	/* important */
	kmalloc_init();

	/* inversion test */
	/*
	for (unsigned i = 0; NUM_TEST_RUNS > i; i++) {
		struct _matrix * a = matrix_alloc(NxN);
		struct _matrix * b = matrix_alloc(NxN);
		struct _matrix * c = matrix_alloc(NxN);

		matirx_random_pos_def(a, RAND_ENTRY_MIN, RAND_ENTRY_MAX);
		matrix_copy(b, a);
		matrix_invert(a);
		matrix_mult(c, b, a);
		if (!almost_identity(c, TEST_PRECISION))
			printf("error\n");
		kfree(a, NxN);
		kfree(b, NxN);
		kfree(c, NxN);
	}
	*/

	/* mult test */
	/*
	struct _matrix * a = matrix_alloc(NxN);
	for (unsigned i = 0; N > i; i++)
		matrix_set_entry(a, (struct _matrix_entry) {i, 0}, i);

	struct _matrix * b = matrix_alloc(Nx1);
	matrix_set_entry(b, (struct _matrix_entry) {0, 0}, 1);

	struct _matrix * c = matrix_alloc(Nx1);
	matrix_mult(c, a, b);

	struct _matrix * d = matrix_alloc(Nx1);
	for (unsigned i = 0; N > i; i++)
		matrix_set_entry(d, (struct _matrix_entry) {i, 0}, i);
	matrix_trans(c);
	matrix_trans(d);

	if (690880 != matrix_scalar_prod(d, c))
		fprintf(stderr, "error in mult test\n");

	kfree(a, NxN);
	kfree(b, NxN);
	kfree(c, NxN);
	kfree(d, NxN);
	*/

	/* newton method */
	/* prepare cost function */
	struct _matrix * p = matrix_alloc(NxN);
	matirx_random_pos_def(p, RAND_ENTRY_MIN, RAND_ENTRY_MAX);

	struct _matrix * q = matrix_alloc(Nx1);
	matrix_zero_up(q);

	double r = 0;

	struct _quadratic_form qf;
	quadratic_form_init(&qf , p, q, r);

	struct _matrix * tmp = 0;
	struct _matrix * dk = matrix_alloc(Nx1);
	struct _matrix * h2_inv = matrix_alloc(NxN);
	matrix_mult(h2_inv, qf.p, qf.p);
	matrix_invert(h2_inv);
	struct _matrix * x = matrix_alloc(Nx1);
	for (unsigned i = 0; N > i; i++)
		matrix_set_entry(x, (struct _matrix_entry) {i, 0}, 1);
	double alpha = 1;
	double c1 = 1e-4;
	double c2 = 0.99;

	for (unsigned i = 0; ITERATIONS > i; i++) {
		tmp = qf.eval_grad(&qf, x);
		matrix_mult(dk, h2_inv, tmp);
		kfree(tmp, Nx1);
		alpha = line_search(&qf, x, dk, c1, c2);
		matrix_scalar_mult(dk, alpha);
		matrix_sub(x, x, dk);
	}
	matrix_print(x);

	kfree(h2_inv, Nx1);
	kfree(dk, Nx1);
	kfree(p, Nx1);
	kfree(q, Nx1);
	kfree(x, Nx1);

	return EXIT_SUCCESS;
}
