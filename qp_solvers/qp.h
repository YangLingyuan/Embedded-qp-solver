#ifndef QP_H
#define QP_H

#include "matrix_ops.h"

#define QUADRATIC_FORM_MAX 1

struct _quadratic_form {
	/* represents the quadratic form: (1/2) x^T.p.x + q^T.x + r */
	struct _matrix * p;
	struct _matrix * q;
	double r;
};

struct _quadratic_form * 
quadratic_form_alloc(struct _matrix * p, struct _matrix * q, double r);
void quadratic_form_free(struct _quadratic_form * qf);
/* gives f(x) (scalar) */
double
quadratic_form_eval(struct _quadratic_form * qf, struct _matrix * x);
/* gives (grad f)(x) (Nx1 matrix (vector), remember to free !) */
struct _matrix *
quadratic_form_eval_grad(struct _quadratic_form * qf, struct _matrix * x);

#endif
