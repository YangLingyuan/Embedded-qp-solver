#include <stdio.h>
#include <stdlib.h>

#include "qp.h"
#include "qp_solvers.h"

#define C1 0.9
#define C2 1e-4

#define ALPHA0_GRAD 1
#define ALPHA0_NEWTON 1

#define MIN_GRAD_GRAD 1e-1
#define MIN_GRAD_NEWTON 1e-1

#define ABSTOL = 1e-4
#define RELTOL = 1e-2
#define ALPHA_ADMM 1.4

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

struct _matrix *
gradient_descent_with_line_search(struct _matrix * x0,
		                  unsigned iterations,
		                  struct _quadratic_form * qf)
{
	struct _matrix * x = matrix_alloc(Nx1);
	matrix_copy(x, x0);

	double alpha = ALPHA0_GRAD;

	struct _matrix * grad_qf = 0;
	struct _matrix * d = matrix_alloc(Nx1);
	for (unsigned i = 0; iterations > i; i++) {
		grad_qf = quadratic_form_eval_grad(qf, x);
		if (!grad_qf) {
			fprintf(stderr, "got invalid grad_qf in test\n");
			exit(EXIT_FAILURE);
		}
		if (MIN_GRAD_GRAD > matrix_norm(grad_qf)) {
			matrix_free(grad_qf);
			break;
		}
		matrix_copy(d, grad_qf);
		matrix_free(grad_qf);

		matrix_scalar_mult(d, -1);

		alpha = line_search(qf, x, d, ALPHA0_GRAD, C1, C2);

		matrix_scalar_mult(d, alpha);
		matrix_add(x, x, d);

	}
	matrix_free(d);

	return x;
}

struct _matrix *
newton_method_with_line_search(struct _matrix * x0,
		               unsigned iterations,
		               struct _quadratic_form * qf)
{
	struct _matrix * x = matrix_alloc(Nx1);
	matrix_copy(x, x0);

	double alpha = ALPHA0_NEWTON;

	struct _matrix * grad_qf = 0;
	struct _matrix * d = matrix_alloc(Nx1);
	struct _matrix * hessian_inv = matrix_alloc(NxN);
	matrix_copy(hessian_inv, qf->p);
	matrix_invert(hessian_inv);

	for (unsigned i = 0; iterations > i; i++) {
		grad_qf = quadratic_form_eval_grad(qf, x);
		if (!grad_qf) {
			fprintf(stderr, "got invalid grad_qf in test\n");
			exit(EXIT_FAILURE);
		}
		if (MIN_GRAD_NEWTON > matrix_norm(grad_qf)) {
			matrix_free(grad_qf);
			break;
		}
		matrix_mult(d, hessian_inv, grad_qf);
		matrix_free(grad_qf);

		matrix_scalar_mult(d, -1);

		alpha = line_search(qf, x, d, ALPHA0_NEWTON, C1, C2);

		matrix_scalar_mult(d, alpha);
		matrix_add(x, x, d);

	}
	matrix_free(hessian_inv);
	matrix_free(d);

	return x;
}

struct _matrix *
admm(struct _matrix * x0, unsigned iterations,
		          struct _quadratic_form * qf)
{
	struct _matrix * x = matrix_alloc(Nx1);
	matrix_copy(x, x0);

	for (unsigned i = 0; iterations > i; i++) {
	}

	return x;
}
