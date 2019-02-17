#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "qp.h"
#include "qp_solvers.h"

#define C1 0.9
#define C2 1e-4

#define ALPHA0_GRAD 1
#define ALPHA0_NEWTON 1

#define MIN_GRAD_GRAD 1e-1
#define MIN_GRAD_NEWTON 1e-1

#define ABSTOL 1e-4
#define RELTOL 1e-2
#define ALPHA_ADMM 1.0

static unsigned armijo(struct _matrix * xk,
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

static double
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

static
void admm_update_x(struct _matrix * R, double rho,
		   struct _matrix * x, struct _matrix * z,
		   struct _matrix * u, struct _matrix * q)
{
	struct _matrix * y1 = matrix_alloc(Nx1);

	matrix_sub(y1, z, u);
	matrix_scalar_mult(y1, rho);
	matrix_sub(y1, y1, q);
	matrix_mult(x, R, y1);

	matrix_free(y1);
}

static
void admm_update_x_hat(struct _matrix * x, double alpha,
		       struct _matrix * z, struct _matrix * x_hat)
{
	struct _matrix * tmp = matrix_alloc(Nx1);

	matrix_copy(x_hat, x);
	matrix_scalar_mult(x_hat, alpha);
	matrix_copy(tmp, z);
	matrix_scalar_mult(tmp, 1 - alpha);
	matrix_add(x_hat, x_hat, tmp);

	matrix_free(tmp);
}

static
void admm_update_z(struct _matrix * ub,
		   struct _matrix * lb,
                   struct _matrix * u,
                   struct _matrix * z,
		   struct _matrix * x_hat)
{
	struct _matrix * tmp = matrix_alloc(Nx1);

	matrix_add(tmp, x_hat, u);
	matrix_max(tmp, tmp, lb);
	matrix_min(z, tmp, ub);

	matrix_free(tmp);
}

static
void admm_update_u(struct _matrix * u,
		   struct _matrix * z,
		   struct _matrix * x_hat)
{
	struct _matrix * tmp = matrix_alloc(Nx1);

	matrix_sub(tmp, x_hat, z);
	matrix_add(u, u, tmp);

	matrix_free(tmp);
}

static
double admm_r_norm(struct _matrix * x, struct _matrix * z)
{
	struct _matrix * tmp = matrix_alloc(Nx1);

	matrix_sub(tmp, x, z);
	double ret = matrix_norm(tmp);

	matrix_free(tmp);

	return ret;
}

static
double admm_s_norm(double rho,
		   struct _matrix * z, struct _matrix * z_old)
{
	struct _matrix * tmp = matrix_alloc(Nx1);

	matrix_sub(tmp, z, z_old);
	matrix_scalar_mult(tmp, -rho);
	double ret = matrix_norm(tmp);

	matrix_free(tmp);

	return ret;
}

static
double admm_eps_pri(unsigned nrows,
		    double abstol, double restol,
		    struct _matrix * x, struct _matrix * z)
{
	double norm_x = matrix_norm(x);
	double norm_z = matrix_norm(z);
	double norm_max = norm_x > norm_z ? norm_x : norm_z;

	return sqrt(nrows) * abstol + restol * norm_max;
}

static
double admm_eps_dual(unsigned nrows,
		     double abstol, double restol,
		     double rho, struct _matrix * u)
{
	double norm_u = matrix_norm(u);

	return sqrt(nrows) * abstol + restol * rho * norm_u;
}

struct _matrix *
admm(struct _matrix * x0, unsigned iterations,
		          struct _quadratic_form * qf)
{
	struct _matrix * x = matrix_alloc(Nx1);
	struct _matrix * x_hat = matrix_alloc(Nx1);
	struct _matrix * z = matrix_alloc(Nx1);
	struct _matrix * u = matrix_alloc(Nx1);
	matrix_zero_up(z);
	matrix_zero_up(u);
	struct _matrix * R = matrix_alloc(NxN);
	struct _matrix * P = qf->p; 
	unsigned nrows = MATRIX_GET_ROW(P);
	struct _matrix * q = qf->q; 
	double r_nrom = 0;
	double s_nrom = 0;
	double eps_pri = 0;
	double eps_dual = 0;
	double abstol = ABSTOL;
	double restol = RELTOL;
	struct _matrix * ub = matrix_alloc(Nx1);
	struct _matrix * lb = matrix_alloc(Nx1);
	for (unsigned i = 0; nrows > i; i++) {
		matrix_set_entry(ub, ME(i, 0), ADMM_BOX_CONSTRAINT_MAX);
		matrix_set_entry(lb, ME(i, 0), ADMM_BOX_CONSTRAINT_MIN);
	}
	double rho = 1;
	double alpha = 1;

	for (unsigned i = 0; iterations > i; i++) {
		if (!i) {
			matrix_copy(R, P);
			for (unsigned i = 0; nrows > i; i++) {
				double s = matrix_get_entry(R, ME(i, i));
				matrix_set_entry(R, ME(i, i), s + rho);
			}
			matrix_invert(R);
		}
		struct _matrix * z_old = matrix_alloc(Nx1);
		matrix_copy(z_old, z);

		admm_update_x(R, rho, x, z, u, q);
		admm_update_x_hat(x, alpha, z, x_hat);
		admm_update_z(ub, lb, u, z, x_hat);
		admm_update_u(u, z, x_hat);

		r_nrom = admm_r_norm(x, z);
		s_nrom = admm_s_norm(rho, z, z_old);
		eps_pri = admm_eps_pri(nrows, abstol, restol, x, z);
		eps_dual = admm_eps_dual(nrows, abstol, restol, rho, u);

		matrix_free(z_old);

		if (r_nrom < eps_pri && s_nrom < eps_dual)
			break;
	}
	matrix_free(R);
	matrix_free(u);
	matrix_free(z);
	matrix_free(lb);
	matrix_free(ub);
	matrix_free(x_hat);

	return x;
}
