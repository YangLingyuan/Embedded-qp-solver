#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define ABS(x) (x < 0 ? -x : x)

#include "qp.h"
#include "kmalloc.h"

double
quadratic_form_eval(struct _quadratic_form * qf, struct _matrix * x)
{
	/* gives (1/2) x^T.p.x + q^T.x + r */
	struct _matrix * tmp = matrix_alloc(Nx1);
	if (!tmp) {
		fprintf(stderr, "no matrix available "
				"to alloc in quadratic_form_eval\n");
		return 0;
	}
	matrix_mult(tmp, qf->p, x);
	double a = 0.5 * matrix_scalar_prod(x , tmp);
	matrix_free(tmp);

	a += matrix_scalar_prod(qf->q, x);
	a += qf->r;

	return a;
}

struct _matrix *
quadratic_form_eval_grad(struct _quadratic_form * qf, struct _matrix * x)
{
	/* gives p.x + q, which is a Nx1 matrix (vector), */
	/* hence u should not forget to matrix_free it afterwards */
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

static void quadratic_form_init(struct _quadratic_form * qf,
		         struct _matrix * p, struct _matrix * q, double r)
{
	qf->p = p;
	qf->q = q;
	qf->r = r;
}

struct _quadratic_form * 
quadratic_form_alloc(struct _matrix * p, struct _matrix * q, double r)
{
	/* alloc and init in one */
	struct _quadratic_form * qf = kmalloc(QUADRATIC_FORM, 0);
	if (!qf) {
		fprintf(stderr, "no quadratic form available to alloc\n");
		return 0;
	}
	quadratic_form_init(qf, p, q, r);

	return qf;
}

void quadratic_form_free(struct _quadratic_form * qf)
{
	kfree(qf, QUADRATIC_FORM);
}
