#ifndef TEST_H
#define TEST_H

#include "matrix_ops.h"
#include "qp.h"
#include "qp_solvers.h"

#ifndef NO_BOARD
#include "config_board.h"
#else
#include "config.h"
#endif

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
