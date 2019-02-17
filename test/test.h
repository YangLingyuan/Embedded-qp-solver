#ifndef TEST_H
#define TEST_H

#include "matrix_ops.h"
#include "qp.h"
#include "qp_solvers.h"

#define TEST_PRECISION 1e-6
#define NUM_INVERSION_TEST_RUNS 8

#define P_RAND_ENTRY_MIN -1e3
#define P_RAND_ENTRY_MAX 1e3
#define Q_RAND_ENTRY_MIN -1e3
#define Q_RAND_ENTRY_MAX 1e3
#define X_RAND_ENTRY_MIN -1e3
#define X_RAND_ENTRY_MAX 1e3

#define NUM_OPT_TEST_RUNS 16

#define GRAD_ITERATIONS 1e4
#define HESS_ITERATIONS 1e1
#define ADMM_ITERATIONS 1e4

#define GRAD_TEST
#define HESS_TEST
#define ADMM_TEST
#define REF_TEST

void inversion_test(void);
void scalar_prod_test(void);
void
test_optimizer(struct _quadratic_form * qf, unsigned iterations,
	       struct _matrix * x0, const char * opt_name,
	       struct _matrix * (* opt)(struct _matrix *,
		                        unsigned,
		                        struct _quadratic_form *));
void test_reference(struct _quadratic_form * qf);

#endif
