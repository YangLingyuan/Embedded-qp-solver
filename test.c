#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

#include "matrix_ops.h"
#include "qp.h"

#define TEST_PRECISION 1e-6
#define NUM_INVERSION_TEST_RUNS 128
#define RAND_ENTRY_MIN -1e4
#define RAND_ENTRY_MAX 1e4

#define C1 0.9
#define C2 1e-4
#define ALPHA0_GRAD 1e1
#define ALPHA_WARM_INCREASE_GRAD 1e7
#define ALPHA0_NEWTON 1e18
#define ALPHA_WARM_INCREASE_NEWTON 1e7
#define NUM_OPT_TEST_RUNS 8

#define MIN_GRAD 8e-1
#define ITERATIONS 1e6

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

unsigned armijo(struct _matrix * xk,
		struct _quadratic_form * qf,
		double fx, struct _matrix * grad_fx,
  		struct _matrix * d_alpha, double c2)
{
	struct _matrix * lhs_arg = matrix_alloc(Nx1);
	matrix_add(lhs_arg, xk, d_alpha);

	double lhs = quadratic_form_eval(qf, lhs_arg);
	matrix_free(lhs_arg);

	double rhs = fx + c2 * matrix_scalar_prod(grad_fx, d_alpha);

	return lhs <= rhs;
}

double
line_search(struct _quadratic_form * qf, struct _matrix * x,
	    struct _matrix * d, double alpha0, double c1, double c2)
{
	double fx = quadratic_form_eval(qf, x);
	struct _matrix * grad_fx = quadratic_form_eval_grad(qf, x);

	struct _matrix * xk = matrix_alloc(Nx1);
	struct _matrix * d_alpha = matrix_alloc(Nx1);
	double alpha = alpha0;

	matrix_copy(d_alpha, d);
	matrix_scalar_mult(d_alpha, alpha);
	matrix_add(xk, x, d_alpha);
	while (!armijo(xk, qf, fx, grad_fx, d_alpha, c2)) {
		alpha *= c1;
		matrix_copy(d_alpha, d);
		matrix_scalar_mult(d_alpha, alpha);
		matrix_add(xk, x, d_alpha);
	}

	matrix_free(grad_fx);
	matrix_free(xk);
	matrix_free(d_alpha);

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
			fprintf(stderr, "error\n");
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

	double should_be = (N*(N+1)*(2*N+1)/6) - N*N;
	if (should_be != matrix_scalar_prod(d, c))
		fprintf(stderr, "error in mult test\n");

	matrix_free(a);
	matrix_free(b);
	matrix_free(c);
	matrix_free(d);
}

void
gradient_descent_with_line_search_test(struct _matrix * x0,
		                       struct _quadratic_form * qf)
{
	time_t start = time(0);
	unsigned end_by_grad_min = 0;

	struct _matrix * x = matrix_alloc(Nx1);
	matrix_copy(x, x0);

	double alpha = ALPHA0_GRAD;
	double c1 = C1;
	double c2 = C2;

	struct _matrix * grad_qf = 0;
	struct _matrix * d = matrix_alloc(Nx1);
	for (unsigned i = 0; ITERATIONS > i; i++) {
		grad_qf = quadratic_form_eval_grad(qf, x);
		if (!grad_qf) {
			fprintf(stderr, "got invalid grad_qf in test\n");
			exit(EXIT_FAILURE);
		}
		if (MIN_GRAD > matrix_norm(grad_qf)) {
			end_by_grad_min = 1U;
			matrix_free(grad_qf);
			break;
		}
		matrix_copy(d, grad_qf);
		matrix_free(grad_qf);

		matrix_scalar_mult(d, -1);

		alpha = line_search(qf, x, d,
				alpha * ALPHA_WARM_INCREASE_GRAD, c1, c2);

		matrix_scalar_mult(d, alpha);
		matrix_add(x, x, d);

	}
	matrix_free(d);

	time_t end = time(0);

	printf("gradient_descent_with_line_search_test gave:\n");
	printf("f(x_min) = %e (took around %li s) %s\n",
	       quadratic_form_eval(qf, x), end - start,
	       end_by_grad_min ? "(end by grad min)" : "");
	fflush(stdout);

	matrix_free(x);
}

void
newton_method_with_line_search_test(struct _matrix * x0,
		                    struct _quadratic_form * qf)
{
	time_t start = time(0);
	unsigned end_by_grad_min = 0;

	struct _matrix * x = matrix_alloc(Nx1);
	matrix_copy(x, x0);

	double alpha = ALPHA0_NEWTON;
	double c1 = C1;
	double c2 = C2;

	struct _matrix * grad_qf = 0;
	struct _matrix * d = matrix_alloc(Nx1);
	struct _matrix * hessian_2_inv = matrix_alloc(NxN);
	matrix_mult(hessian_2_inv, qf->p, qf->p);
	matrix_invert(hessian_2_inv);

	for (unsigned i = 0; ITERATIONS > i; i++) {
		grad_qf = quadratic_form_eval_grad(qf, x);
		if (!grad_qf) {
			fprintf(stderr, "got invalid grad_qf in test\n");
			exit(EXIT_FAILURE);
		}
		if (MIN_GRAD > matrix_norm(grad_qf)) {
			end_by_grad_min = 1U;
			matrix_free(grad_qf);
			break;
		}
		matrix_mult(d, hessian_2_inv, grad_qf);
		matrix_free(grad_qf);

		matrix_scalar_mult(d, -1);

		alpha = line_search(qf, x, d,
				alpha * ALPHA_WARM_INCREASE_NEWTON, c1, c2);

		matrix_scalar_mult(d, alpha);
		matrix_add(x, x, d);

	}
	matrix_free(hessian_2_inv);
	matrix_free(d);

	time_t end = time(0);

	printf("newton_method_with_line_search_test gave:\n");
	printf("f(x_min) = %e (took around %li s) %s\n",
	       quadratic_form_eval(qf, x), end - start,
	       end_by_grad_min ? "(end by grad min)" : "");
	fflush(stdout);

	matrix_free(x);
}

int main(void)
{
	/* important */
	kmalloc_init();

	inversion_test();
	mult_test();

	/* prepare quadratic cost function */
	struct _matrix * p = matrix_alloc(NxN);
	struct _matrix * q = matrix_alloc(Nx1);
	double r = 0;
	struct _quadratic_form * qf = quadratic_form_alloc(p, q, r);
	if (!qf) {
		fprintf(stderr, "invalid quadratic "
				"from quadratic_form_alloc in test\n");
		return EXIT_FAILURE;
	}

	/* to hold initial state */
	struct _matrix * x0 = matrix_alloc(Nx1);

	for (unsigned i = 0; NUM_OPT_TEST_RUNS > i; i++) {
		matirx_random_pos_def(p, RAND_ENTRY_MIN, RAND_ENTRY_MAX);
		matrix_random(q, RAND_ENTRY_MIN, RAND_ENTRY_MAX);
		r = rand();
		matrix_random(x0, RAND_ENTRY_MIN, RAND_ENTRY_MAX);
		gradient_descent_with_line_search_test(x0, qf);
		newton_method_with_line_search_test(x0, qf);
		printf("\v");
	}

	/* clean up */
	matrix_free(x0);
	matrix_free(p);
	matrix_free(q);
	quadratic_form_free(qf);

	return EXIT_SUCCESS;
}
